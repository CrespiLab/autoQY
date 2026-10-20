from pathlib import Path
import base64
import json
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
try:
    from dash import dcc, html
    from dash._callback_context import context_value
    from dash._utils import AttributeDict
except ImportError:
    dcc = html = None

from autoqy_core.smoother import SpectralDataset
from autoqy_core.tools.spectral_recovery import (
    SESSION_MAGIC, SpectralRecovery, decode_session, encode_session,
)
from autoqy_core.tools.smoother_gui import (
    _append_packed, _combine_loaded, _compact_packed, _loaded_spectrum_rows,
    _pack, _remove_packed, _unpack, create_app,
)


def components(root):
    if isinstance(root, (list, tuple)):
        for item in root:
            yield from components(item)
    else:
        yield root
        yield from components(getattr(root, "children", []) or [])


def spectrum(value=1):
    return SpectralDataset(np.array([400., 410.]), np.array([0.]),
                           np.array([[value], [value + 1.]]), "csv")


@unittest.skipUnless(html is not None, "Dash GUI dependencies are not installed")
class SpectralRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "session.json"
        self.data = _compact_packed(_pack(spectrum(), ["sample"], ["sample.csv"]))

    def test_atomic_failure_preserves_previous_snapshot(self):
        store = SpectralRecovery(self.path)
        snapshot = {"version": 1, "dataset": self.data}
        store.write(snapshot)
        with patch("autoqy_core.tools.spectral_recovery.os.replace", side_effect=OSError):
            with self.assertRaises(OSError):
                store.write({"version": 1, "dataset": None})
        self.assertEqual(store.read(), snapshot)
        self.assertEqual(list(self.path.parent.glob("*.tmp")), [])

    def test_binary_roundtrip_preserves_full_precision_and_compresses_arrays(self):
        matrix = np.random.default_rng(4).random((1000, 100))
        snapshot = {"version": 1, "dataset": {"absorbance": matrix.tolist()},
                    "settings": {"baseline-enabled": ["on"]},
                    "isomer_request": {"reference": 0, "mixture": 99, "percent": 42}}
        payload = encode_session(snapshot)
        self.assertTrue(payload.startswith(SESSION_MAGIC))
        self.assertLess(len(payload), len(json.dumps(snapshot).encode()) / 2)
        self.assertEqual(decode_session(payload), snapshot)
        np.testing.assert_array_equal(decode_session(payload)["dataset"]["absorbance"], matrix)
        with self.assertRaises(ValueError):
            decode_session(payload[:50])
        with self.assertRaises(ValueError):
            encode_session({"version": 2})

    def test_legacy_autosave_is_read_and_new_writes_are_binary(self):
        legacy = self.path.parent / "autosave.json"
        snapshot = {"version": 1, "dataset": self.data}
        legacy.write_text(json.dumps(snapshot), encoding="utf-8")
        store = SpectralRecovery(legacy.with_suffix(".aqys"))
        self.assertEqual(store.read(), snapshot)
        store.write(snapshot)
        self.assertTrue(store.path.read_bytes().startswith(SESSION_MAGIC))
        self.assertTrue(legacy.exists())

    def test_manual_save_uses_output_and_preserves_dynamic_settings(self):
        app = create_app(self.path)
        by_id = {item.id: item for item in components(app.layout)
                 if isinstance(getattr(item, "id", None), str)}
        callback = next(value for key, value in app.callback_map.items()
                        if "session-save-message.children" in key)
        settings = {"baseline-enabled": ["on"]}
        values = [settings.get(item["id"], getattr(by_id[item["id"]], "value", None))
                  for item in callback["state"][4:-6]]
        args = (1, self.data, None, "source", None, *values,
                [0.002], [0.5], [[]], ["Custom legend"],
                [{"type": "legend-spectrum", "index": 0}],
                [{"type": "legend-name", "index": 0}])
        token = context_value.set(AttributeDict(triggered_inputs=[{"prop_id": "save-session.n_clicks"}]))
        try:
            save = callback["callback"].__wrapped__
            destination = self.path.parent / "portable.aqys"
            with patch("autoqy_core.tools.smoother_gui._choose_save_file", return_value=str(destination)):
                result = save(*args)
            self.assertIn("Saved session", result)
            snapshot = decode_session(destination.read_bytes())
            self.assertEqual(snapshot["settings"]["baseline-enabled"], ["on"])
            self.assertEqual(snapshot["dataset"]["legend_names"], ["Custom legend"])
            self.assertEqual(snapshot["dataset"]["legend_visibility"], [False])
            self.assertEqual(snapshot["dataset"]["concentrations"], [0.002])
            self.assertEqual(snapshot["dataset"]["path_lengths"], [0.5])
            previous = destination.read_bytes()
            with patch("autoqy_core.tools.smoother_gui._choose_save_file", return_value=None):
                self.assertIn("cancelled", save(*args))
            self.assertEqual(destination.read_bytes(), previous)
        finally:
            context_value.reset(token)

    def test_session_drop_replaces_data_restores_settings_and_rejects_stale_autosave(self):
        app = create_app(self.path)
        before = {item.id: item for item in components(app.layout)
                  if isinstance(getattr(item, "id", None), str)}
        old_generation = before["session-generation"].data
        snapshot = {"version": 1, "dataset": self.data, "source": "saved-folder",
                    "settings": {"save-folder": self.directory.name, "baseline-enabled": ["on"],
                                 "slice-time-multiplier": 30},
                    "isomer_request": {"reference": 0, "mixture": 0, "percent": 50}}
        key, callback = next((key, value) for key, value in app.callback_map.items()
                             if "dataset-store.data" in key and "load-message.children" in key)
        overrides = {("upload-spectra", "contents"): ["data:application/octet-stream;base64," +
                     base64.b64encode(encode_session(snapshot)).decode()],
                     ("upload-spectra", "filename"): ["saved.aqys"]}
        # Use the real Dash request context so the reload side-effect is exercised.
        def dependency(item):
            identity, prop = item["id"], item["property"]
            value = overrides.get((identity, prop))
            if (identity, prop) not in overrides:
                value = getattr(before.get(identity), prop, None) if not identity.startswith("{") else []
            return {**item, "value": value}
        output = callback["output"]
        client = app.server.test_client()
        response = client.post("/_dash-update-component", json={
            "output": key, "outputs": [{"id": item.component_id, "property": item.component_property}
                                         for item in output],
            "inputs": [dependency(item) for item in callback["inputs"]],
            "state": [dependency(item) for item in callback["state"]],
            "changedPropIds": ["upload-spectra.contents"],
        })
        self.assertEqual(response.status_code, 200, response.data)
        self.assertIn("session-reload", response.json["sideUpdate"])
        after = {item.id: item for item in components(app.layout)
                 if isinstance(getattr(item, "id", None), str)}
        self.assertEqual(after["dataset-store"].data, self.data)
        self.assertEqual(after["baseline-enabled"].value, ["on"])
        self.assertNotIn("save-folder", after)
        self.assertEqual(after["isomer-request-store"].data, snapshot["isomer_request"])
        stale = client.post("/_autoqy_spectral_autosave", json={
            "version": 1, "dataset": None, "generation": old_generation})
        self.assertEqual(stale.status_code, 409)
        self.assertEqual(SpectralRecovery(self.path).read(), snapshot)

    def test_saved_session_restores_original_data_and_settings_after_restart(self):
        app = create_app(self.path)
        snapshot = {"version": 1, "dataset": self.data, "source": "source-folder",
                    "settings": {"slice-time-multiplier": 30, "baseline-enabled": ["on"],
                                 "wavelength-low": 402, "wavelength-high": 408}}
        response = app.server.test_client().post("/_autoqy_spectral_autosave", json=snapshot)
        self.assertEqual(response.status_code, 204)
        restarted = create_app(self.path)
        controls = {item.id: item for item in components(restarted.layout)
                    if isinstance(getattr(item, "id", None), str)}
        self.assertEqual(controls["slice-time-multiplier"].value, 30)
        self.assertEqual(controls["wavelength-low"].value, 402)
        self.assertEqual(controls["baseline-enabled"].value, ["on"])
        np.testing.assert_array_equal(_unpack(controls["dataset-store"].data).absorbance,
                                      spectrum().absorbance)
        self.assertIn("Recovered", controls["autosave-message"].children)
        self.assertEqual(controls["edit-legend"].value, [])

    def test_cleared_session_does_not_resurrect_spectra(self):
        SpectralRecovery(self.path).write({"version": 1, "dataset": None})
        app = create_app(self.path)
        store = next(item for item in components(app.layout)
                     if getattr(item, "id", None) == "dataset-store")
        self.assertIsNone(store.data)

    def test_corrupt_file_does_not_prevent_startup_or_get_overwritten(self):
        self.path.write_text('{broken', encoding='utf-8')
        app = create_app(self.path)
        message = next(item for item in components(app.layout)
                       if getattr(item, "id", None) == "autosave-message")
        self.assertIn("could not be read", message.children)
        self.assertEqual(self.path.read_text(), '{broken')

    def test_failed_autosave_returns_visible_failure_status(self):
        app = create_app(self.path)
        with patch.object(SpectralRecovery, "write", side_effect=OSError):
            response = app.server.test_client().post(
                "/_autoqy_spectral_autosave", json={"version": 1, "dataset": self.data})
        self.assertEqual(response.status_code, 500)

    def test_legend_editors_are_not_created_until_requested_and_edits_survive_toggle(self):
        rows = _loaded_spectrum_rows(html, dcc, self.data)
        ids = [item.id for item in components(rows) if isinstance(getattr(item, "id", None), dict)]
        self.assertFalse(any(item['type'].startswith('legend-') for item in ids))
        app = create_app(self.path)
        toggle = next(value['callback'].__wrapped__ for key, value in app.callback_map.items()
                      if 'show-all-legends.disabled' in key)
        updated, disabled, _ = toggle([], self.data, [], [], [[]], ['Custom'],
                                     [{'type': 'legend-spectrum', 'index': 0}],
                                     [{'type': 'legend-name', 'index': 0}])
        self.assertTrue(disabled)
        self.assertEqual(updated['legend_names'], ['Custom'])
        rows = _loaded_spectrum_rows(html, dcc, updated, edit_legend=True)
        name = next(item for item in components(rows)
                    if isinstance(getattr(item, 'id', None), dict) and item.id['type'] == 'legend-name')
        self.assertEqual(name.value, 'Custom')

    def test_multiple_drop_then_single_and_batch_additions_keep_unit_spacing(self):
        dataset, labels, _ = _combine_loaded([(spectrum(i), f'sample{i}.csv') for i in range(4)])
        data = _pack(dataset, labels, labels)
        single = _pack(spectrum(4), ['extra'], ['extra.csv'])
        data = _append_packed(data, single, elapsed_seconds=183.71)
        self.assertEqual(data['coordinates'], [0., 1., 2., 3., 4.])
        # Removing intermediate samples must not change the next unit step.
        data = _remove_packed(data, [1, 2, 3])
        data = _append_packed(data, single, elapsed_seconds=901)
        self.assertEqual(data['coordinates'], [0., 4., 5.])
        data = _append_packed(data, _pack(dataset, labels, labels), elapsed_seconds=1000)
        self.assertEqual(data['coordinates'], [0., 4., 5., 6., 7., 8., 9.])
        np.testing.assert_array_equal(np.array(data['coordinates']) * 30,
                                      [0, 120, 150, 180, 210, 240, 270])

    def test_existing_explicit_time_axis_keeps_its_spacing(self):
        series = SpectralDataset(np.array([400., 410.]), np.array([0., 10., 20.]),
                                 np.ones((2, 3)), 'csv')
        data = _append_packed(_pack(series, ['a', 'b', 'c'], ['series.csv']), self.data,
                              elapsed_seconds=187)
        self.assertEqual(data['coordinates'], [0., 10., 20., 30.])
