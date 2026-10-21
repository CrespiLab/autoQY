import unittest
import base64
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
import shutil
import subprocess

from autoqy_core.epsilon import EpsilonResult, calculate_epsilon_statistics, load_epsilon_table
from autoqy_core.plot_style import COLORBLIND_TRACE_PALETTE
from autoqy_core.smoother import SpectralDataset, reconstruct_isomer_absorbance
from autoqy_core.tools.smoother_gui import (
    MAX_INTERACTIVE_SPECTRA,
    _append_packed,
    _absorbance_figure,
    _colors,
    _combine_loaded,
    _compact_packed,
    _csv_filename,
    _csv_exports,
    _epsilon_figure,
    _export_csv_payload,
    _fit_exponential_decay,
    _figure_sheet_page,
    _loaded_spectrum_rows,
    _pack,
    _pack_epsilon,
    _parameter_cards,
    _plotted_spectrum_indices,
    _prepare_processing,
    _remove_packed,
    _reorder_packed,
    _spectrum_colors,
    _save_path,
    _unpack,
    _with_spectrum_state,
    _wavelength_slice,
    _wavelength_slice_figure,
    create_app,
)

try:
    from dash import dcc, html
    from dash._callback_context import context_value
    from dash._utils import AttributeDict
except ImportError:
    dcc = html = None


def _components(root):
    if isinstance(root, (list, tuple)):
        for child in root:
            yield from _components(child)
        return
    yield root
    children = getattr(root, "children", None)
    if children is not None:
        yield from _components(children)


class IsomerReconstructionTests(unittest.TestCase):
    def test_known_mixture_recovers_pure_isomer_without_epsilon(self):
        stable = np.array([0.2, 0.8, 0.1])
        metastable = np.array([0.9, 0.3, -0.02])
        mixture = 0.25 * stable + 0.75 * metastable
        values = np.column_stack([stable, mixture])
        np.testing.assert_allclose(reconstruct_isomer_absorbance(values, 0, 1, 75),
                                   metastable)
        np.testing.assert_array_equal(values[:, 1], mixture)
        np.testing.assert_array_equal(reconstruct_isomer_absorbance(values, 0, 1, 100),
                                      mixture)

    def test_invalid_compositions_and_spectrum_selections_are_rejected(self):
        values = np.ones((3, 2))
        for percent in (None, 0, -1, 101, float("nan"), float("inf")):
            with self.subTest(percent=percent), self.assertRaises(ValueError):
                reconstruct_isomer_absorbance(values, 0, 1, percent)
        for reference, mixture in ((0, 0), (0, 2), (-1, 1), (None, 1)):
            with self.subTest(reference=reference), self.assertRaises(ValueError):
                reconstruct_isomer_absorbance(values, reference, mixture, 50)


class ProcessedAbsorbanceExportTests(unittest.TestCase):
    def setUp(self):
        self.dataset = SpectralDataset(
            np.array([400.0, 401.0, 402.0]),
            np.array([0.0]),
            np.array([[0.1], [0.2], [0.3]]),
            source_format="csv",
        )
        self.packed = _pack(self.dataset, ["irradiation"], ["irradiation.csv"])

    def test_processed_absorbance_has_a_csv_export_without_epsilon(self):
        filename, text = _export_csv_payload(None, self.packed)
        self.assertEqual(filename, "processed-absorbance.csv")
        self.assertTrue(text.startswith("Wavelength,irradiation\n"))
        self.assertIn("2.00000000e-01", text)

    def test_negative_processed_absorbance_can_be_converted_to_zero(self):
        self.packed["absorbance"][1][0] = -0.2
        _, text = _export_csv_payload(None, self.packed, nonnegative=True)
        frame = pd.read_csv(StringIO(text))
        np.testing.assert_array_equal(frame["irradiation"], [0.1, 0.0, 0.3])

    def test_csv_headers_use_edited_legends_and_quote_special_characters(self):
        self.packed["legend_names"] = ['Sample, "edited"']
        _, text = _export_csv_payload(None, self.packed)
        frame = pd.read_csv(StringIO(text))
        self.assertEqual(frame.columns.tolist(), ["Wavelength", 'Sample, "edited"'])
        np.testing.assert_allclose(frame.iloc[:, 1], self.dataset.absorbance[:, 0])
        result = calculate_epsilon_statistics(self.dataset.wavelengths, self.dataset.absorbance,
                                             np.array([0.001]), np.array([1.]))
        _, text = _export_csv_payload(_pack_epsilon(result, self.packed["legend_names"]), None)
        restored, labels = load_epsilon_table(text)
        self.assertEqual(labels, self.packed["legend_names"])
        np.testing.assert_allclose(restored.absorbance, result.absorbance)

    def test_sheet_pages_only_plotted_data_and_preserves_curve_names(self):
        x = np.arange(105, dtype=float)
        figure = {"data": [{"x": x, "y": x * 2, "name": "Edited legend"},
                           {"x": x, "y": x * 3, "name": "Hidden", "visible": "legendonly"},
                           {"x": x, "y": x * 4, "name": "Reconstructed pure isomer"}],
                  "layout": {"xaxis": {"title": {"text": "Wavelength (nm)"}}}}
        columns, rows, count, message = _figure_sheet_page(figure, 1, 25)
        self.assertEqual([item["name"] for item in columns],
                         ["Wavelength (nm)", "Edited legend", "Reconstructed pure isomer"])
        self.assertEqual(count, 5)
        self.assertEqual(len(rows), 25)
        self.assertEqual(rows[0], {"x": 25, "c0": 50, "c1": 100})
        self.assertIn("105 rows", message)
        self.assertEqual(len(_figure_sheet_page(figure, 4, 25)[1]), 5)
        self.assertEqual(len(_figure_sheet_page(figure, visibility={"0": "legendonly"})[0]), 2)
        _, all_rows, count, _ = _figure_sheet_page(figure, size=None)
        self.assertEqual(len(all_rows), 105)
        self.assertEqual(count, 1)
        self.assertEqual(all_rows[-1]["x"], 104)
        self.assertTrue(all(column["selectable"] for column in columns))

    def test_sheet_accepts_plotly_binary_arrays_and_empty_plots(self):
        array = np.array([400., 401.], dtype="float64")
        encoded = {"dtype": "f8", "bdata": base64.b64encode(array.tobytes()).decode()}
        columns, rows, _, _ = _figure_sheet_page({"data": [{"x": encoded, "y": encoded,
                                                           "name": "sample"}]})
        self.assertEqual(rows[1], {"x": 401., "c0": 401.})
        self.assertEqual(columns[1]["name"], "sample")
        self.assertEqual(_figure_sheet_page(None)[0:3], ([], [], 0))

    def test_negative_absorbance_and_epsilon_can_be_converted_to_zero(self):
        result = EpsilonResult(
            wavelengths=np.array([400.0, 401.0]),
            absorbance=np.array([[-0.1], [0.2]]),
            concentrations_m=np.array([1e-4]),
            path_lengths_cm=np.array([1.0]),
            individual=np.array([[-1000.0], [2000.0]]),
            mean=np.array([-1000.0, 2000.0]),
            standard_deviation=np.array([10.0, 20.0]),
            standard_error=np.array([10.0, 20.0]),
        )
        _, text = _export_csv_payload(
            _pack_epsilon(result, ["sample"]), None, nonnegative=True
        )
        frame = pd.read_csv(StringIO(text))
        for column in ("Absorbance__sample", "Epsilon_M-1_cm-1__sample",
                       "Epsilon_mean_M-1_cm-1"):
            self.assertGreaterEqual(frame[column].min(), 0.0)

    def test_custom_csv_name_accepts_an_optional_extension(self):
        self.assertEqual(_csv_filename("treated irradiation", "default.csv"),
                         "treated irradiation.csv")
        self.assertEqual(_csv_filename("treated.CSV", "default.csv"), "treated.CSV")
        self.assertEqual(_csv_filename("", "default.csv"), "default.csv")

    def test_custom_csv_name_rejects_paths_and_other_extensions(self):
        with self.assertRaisesRegex(ValueError, "without a folder path"):
            _csv_filename("subfolder/treated.csv", "default.csv")
        with self.assertRaisesRegex(ValueError, "must use the .csv extension"):
            _csv_filename("treated.tsv", "default.csv")

    def test_pasted_save_folder_accepts_copy_as_path_and_whitespace(self):
        with TemporaryDirectory(prefix="AutoQY output ") as directory:
            self.assertEqual(_save_path(f'  "{directory}"  ', "output.csv"),
                             Path(directory) / "output.csv")

    def test_large_plot_selection_is_even_and_includes_both_endpoints(self):
        selected = _plotted_spectrum_indices(571)
        self.assertEqual(len(selected), MAX_INTERACTIVE_SPECTRA)
        self.assertEqual(selected[0], 0)
        self.assertEqual(selected[-1], 570)
        self.assertEqual(selected, sorted(set(selected)))
        spacing = np.diff(selected)
        self.assertLessEqual(int(spacing.max() - spacing.min()), 1)

    def test_small_plot_selection_keeps_every_spectrum(self):
        self.assertEqual(_plotted_spectrum_indices(4), [0, 1, 2, 3])

    def test_compact_browser_store_round_trips_absorbance(self):
        compact = _compact_packed(self.packed)
        self.assertNotIn("absorbance", compact)
        self.assertEqual(compact["absorbance_encoding"], "zlib-base64-float64")
        restored = _unpack(compact)
        np.testing.assert_array_equal(restored.wavelengths, self.dataset.wavelengths)
        np.testing.assert_array_equal(restored.coordinates, self.dataset.coordinates)
        np.testing.assert_array_equal(restored.absorbance, self.dataset.absorbance)


@unittest.skipUnless(html is not None, "Dash GUI dependencies are not installed")
class SpectralGuiTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.app = create_app(Path(directory.name) / "autosave.json")

    def test_isomer_overlay_uses_preprocessing_and_visibility_via_dash_endpoint(self):
        dataset = SpectralDataset(np.array([400., 401., 402.]), np.array([0., 1.]),
                                  np.array([[0.2, 0.5], [0.8, 0.4], [0.1, 0.1]]), "csv")
        packed = _pack(dataset, ["stable", "PSS"], ["spectra.csv"])
        key, callback = next((key, callback) for key, callback in self.app.callback_map.items()
                             if "processed-store.data" in key)
        values = {"dataset-store": packed, "wavelength-low": 400., "wavelength-high": 402.,
                  "baseline-enabled": ["on"], "baseline-low": 401.5, "baseline-high": 402.5,
                  "smoothing-method": "off", "savgol-window": 5, "savgol-order": 2,
                  "svd-enabled": [], "svd-rank": 1, "minimal-spectrum-colors": [],
                  "show-original-spectra": [], "main-x-axis-label": "Wavelength (nm)",
                  "main-absorbance-axis-label": "Absorbance", "main-epsilon-axis-label": "ε",
                  "isomer-request-store": {"reference": 0, "mixture": 1, "percent": 50},
                  "show-reconstructed-isomer": ["on"]}
        inputs = [dict(item, value=values.get(item["id"], [])) for item in callback["inputs"]]
        state = [dict(item, value=[]) for item in callback["state"]]
        payload = {"output": key, "outputs": [item.to_dict() for item in callback["output"]],
                   "inputs": inputs, "state": state,
                   "changedPropIds": ["isomer-request-store.data"]}
        client = self.app.server.test_client()
        response = client.post("/_dash-update-component", json=payload)
        self.assertEqual(response.status_code, 200, response.data)
        self.assertEqual(response.json["response"]["preview-error"]["children"], "")
        figure = response.json["response"]["epsilon-plot"]["figure"]
        y = figure["data"][-1]["y"]
        if isinstance(y, dict):
            y = np.frombuffer(base64.b64decode(y["bdata"]), dtype=y["dtype"])
        np.testing.assert_allclose(y, [0.7, -0.1, 0.0], atol=1e-14)
        self.assertEqual(figure["data"][-1]["name"], "Reconstructed pure metastable isomer")
        self.assertEqual(figure["data"][-1]["line"]["dash"], "dot")
        next(item for item in inputs if item["id"] == "show-reconstructed-isomer")["value"] = []
        hidden = client.post("/_dash-update-component", json=payload)
        self.assertEqual(hidden.status_code, 200, hidden.data)
        self.assertEqual(len(hidden.json["response"]["epsilon-plot"]["figure"]["data"]), 2)

    def test_reconstruction_has_short_copy_and_nearby_accessible_help(self):
        components = list(_components(self.app.layout))
        by_id = {item.id: item for item in components if isinstance(getattr(item, "id", None), str)}
        self.assertEqual(by_id["upload-spectra"].children.children, "Drop spectra or a saved session")
        helpers = [item for item in components if getattr(item, "className", None) == "info-popup"]
        topics = {item.to_plotly_json()["props"]["aria-label"]: item for item in helpers}
        for topic in ("Spectral data", "Stable spectrum", "PSS spectrum", "Metastable percentage",
                      "Composition error", "ε reference", "Reference error", "Reconstruct", "Copy selection"):
            helper = topics[f"Help: {topic}"]
            self.assertIsInstance(helper, html.Button)
            self.assertEqual(helper.type, "button")
            tooltip = helper.children[1]
            self.assertEqual(tooltip.role, "tooltip")
            self.assertEqual(helper.to_plotly_json()["props"]["aria-describedby"], tooltip.id)
        self.assertEqual(len({item.children[1].id for item in helpers}), len(helpers))
        for removed in ("Load files", "Negative values in CSV", "Save CSV", "Save figures"):
            self.assertNotIn(f"Help: {removed}", topics)
        self.assertIn("One .aqys session", topics["Help: Spectral data"].children[1].children)
        self.assertIn("one dialog for each", topics["Help: Save data, sessions and figures"].children[1].children)
        help_text = " ".join(str(item.children[1].children) for item in helpers)
        self.assertIn("same total concentration and path length", help_text)
        self.assertIn("They do not include PSS measurement error", help_text)
        self.assertIn("1; enter 0", help_text)

        def visible_text(node):
            if isinstance(node, str):
                return node
            if isinstance(node, (list, tuple)):
                return " ".join(visible_text(child) for child in node)
            if getattr(node, "className", None) == "info-popup":
                return ""
            return visible_text(getattr(node, "children", "")) if node else ""

        reconstruction_copy = visible_text(by_id["reconstruction-panel"])
        self.assertNotIn("envelope", reconstruction_copy)
        self.assertNotIn("preprocessing", reconstruction_copy)
        self.assertNotIn("Drop a triplicate", reconstruction_copy)
        self.assertIn("Metastable (%)", reconstruction_copy)
        self.assertLess(len(reconstruction_copy.split()), 50)

    def test_isomer_status_is_short_keeps_warnings_and_defers_calibrated_status(self):
        status = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                      if "isomer-message.children" in key)
        packed = _pack(SpectralDataset(np.array([400., 401.]), np.array([0., 1.]),
                                       np.array([[.2, .5], [.8, .1]]), "csv"),
                       ["stable", "PSS"], ["spectra.csv"])
        request = {"reference": 0, "mixture": 1, "percent": 50, "error_percent": 1}
        message, disabled = status(packed, request)
        self.assertFalse(disabled)
        self.assertTrue(message.startswith("Metastable absorbance ready."))
        self.assertIn("1 negative value", message)
        self.assertIn("Check composition and baseline", message)
        self.assertNotIn("uncertainty", message)
        self.assertIn("below 5%", status(packed, request | {"percent": 2})[0])
        self.assertEqual(status(packed, request | {"mode": "calibrated"}), ("", True))

    def test_isomer_export_retains_negatives_and_uses_save_as(self):
        packed = _pack(SpectralDataset(np.array([400., 401.]), np.array([0., 1.]),
                                      np.array([[0.2, 0.5], [0.8, 0.1]]), "csv"),
                       ["stable", "PSS"], ["spectra.csv"])
        request = {"reference": 0, "mixture": 1, "percent": 50}
        save = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                    if "epsilon-save-message.children" in key)
        with TemporaryDirectory() as directory:
            token = context_value.set(AttributeDict(
                triggered_inputs=[{"prop_id": "export-epsilon.n_clicks"}]))
            try:
                path = Path(directory) / "reconstructed-metastable-isomer.csv"
                with patch("autoqy_core.tools.smoother_gui._choose_save_file", return_value=str(path)):
                    result = save(1, None, packed, [], "isomer", request)
                self.assertEqual(result[1], "", result)
                frame = pd.read_csv(path)
                self.assertEqual(frame.columns[1], "Reconstructed pure metastable isomer")
                np.testing.assert_allclose(frame.iloc[:, 1], [0.8, -0.6])
                with patch("autoqy_core.tools.smoother_gui._choose_save_file", return_value=None):
                    self.assertIn("cancelled", save(2, None, packed, [], "isomer", request)[0])
            finally:
                context_value.reset(token)

    def test_slice_csv_uses_output_name_clamping_and_fit_details(self):
        slice_data = {"wavelength": 405, "coordinates": [0, 30], "values": [0.4, -0.1],
                      "x_label": "Time (s)", "y_label": "Absorbance", "fit_values": [0.4, -0.05],
                      "lifetime": 12, "lifetime_error": 0.5, "time_unit": "s"}
        name, text = _csv_exports("slice", None, None, "chosen", True,
                                  slice_data=slice_data)[0]
        self.assertEqual(name, "chosen.csv")
        self.assertIn("30,0,0", text)
        self.assertIn("Lifetime,12,s", text)
        self.assertIn("Lifetime standard error,0.5,s", text)
        self.assertEqual(_csv_exports("slice", None, None, slice_data=slice_data)[0][0],
                         "wavelength-slice-405-nm.csv")

    def test_changing_isomer_selection_or_dataset_invalidates_request(self):
        callback = self.app.callback_map["isomer-request-store.data"]["callback"].__wrapped__
        for trigger in ("dataset-store.data", "isomer-percent.value", "isomer-reference.value",
                        "isomer-mixture.value"):
            token = context_value.set(AttributeDict(triggered_inputs=[{"prop_id": trigger}]))
            try:
                self.assertIsNone(callback(1, {"labels": ["stable", "PSS"]}, 0, 1, 50))
            finally:
                context_value.reset(token)

    def test_only_data_panel_is_open_initially(self):
        details = []

        def visit(component):
            if isinstance(component, (list, tuple)):
                for child in component:
                    visit(child)
                return
            if isinstance(component, html.Details):
                details.append(component)
            children = getattr(component, "children", None)
            if children is not None:
                visit(children)

        visit(self.app.layout)
        opened = [
            component for component in details
            if getattr(component, "open", False) is True
        ]
        self.assertEqual(len(opened), 1)
        self.assertIn("1 · Data", str(opened[0].to_plotly_json()))

    def test_large_data_preprocessing_inputs_are_debounced(self):
        components = {
            component.id: component
            for component in _components(self.app.layout)
            if isinstance(getattr(component, "id", None), str)
        }
        for component_id in (
            "wavelength-low", "wavelength-high", "baseline-low", "baseline-high",
            "savgol-window", "savgol-order",
        ):
            self.assertTrue(components[component_id].debounce, component_id)

    def test_native_save_dialog_has_extension_overwrite_prompt_and_remembers_folder(self):
        from autoqy_core.tools.smoother_gui import _choose_save_file
        with patch("autoqy_core.tools.smoother_gui._run_powershell_dialog", return_value=["C:/chosen.csv"]) as dialog:
            self.assertEqual(_choose_save_file("suggested.csv", "C:/source"), "C:/chosen.csv")
            script = dialog.call_args.args[0]
            self.assertIn("SaveFileDialog", script)
            self.assertIn("OverwritePrompt = $true", script)
            self.assertIn("DefaultExt = 'csv'", script)
        save = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                    if "epsilon-save-message.children" in key)
        with TemporaryDirectory() as directory:
            data = _pack(SpectralDataset(np.array([400., 401.]), np.array([0.]),
                                         np.array([[.2], [.3]]), "csv"), ["stable"], ["input.csv"])
            path = str(Path(directory) / "chosen.csv")
            with patch("autoqy_core.tools.smoother_gui._choose_save_file", return_value=path):
                self.assertIn("Saved", save(1, None, data, [], source="C:/input")[0])
            with patch("autoqy_core.tools.smoother_gui._choose_save_file", return_value=None) as dialog:
                self.assertIn("cancelled", save(2, None, data, [], source="C:/other")[0])
                self.assertEqual(dialog.call_args.args[1], directory)

    def test_svd_analysis_is_skipped_while_svd_is_off(self):
        callback = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "svd-rank.options" in key
        )
        dataset = SpectralDataset(
            np.array([400.0, 401.0, 402.0]), np.array([0.0]),
            np.array([[0.1], [0.2], [0.3]]), source_format="csv",
        )
        with patch("autoqy_core.tools.smoother_gui.analyze_svd") as analyze:
            result = callback(
                _compact_packed(_pack(dataset, ["sample"], ["sample.csv"])),
                400.0, 402.0, [], None, None, "off", 5, 2, [],
            )
        analyze.assert_not_called()
        self.assertEqual(result, ([], None, True, "SVD is off.", ""))

    def test_plot_controls_offer_legend_toggle_slice_and_exports(self):
        components = list(_components(self.app.layout))
        by_id = {
            component.id: component for component in components
            if getattr(component, "id", None)
        }
        self.assertNotIn("spectrum-legend-labels", by_id)
        self.assertNotIn("show-spectrum-legend", by_id)
        self.assertIn("loaded-spectrum-manager", by_id)
        self.assertIn("show-all-legends", by_id)
        self.assertIn("hide-all-legends", by_id)
        self.assertEqual(by_id["minimal-spectrum-colors"].value, [])
        self.assertEqual(by_id["show-original-spectra"].value, [])
        self.assertEqual(by_id["include-plot-title"].value, [])
        self.assertEqual(by_id["include-plot-legend"].value, ["on"])
        self.assertEqual(by_id["include-export-axis-labels"].value, ["on"])
        self.assertEqual(by_id["origin-epsilon-export"].value, ["on"])
        self.assertEqual(by_id["show-epsilon-export-grid"].value, [])
        self.assertEqual(by_id["zero-epsilon-export-y"].value, ["on"])
        self.assertNotIn("include-slice-title", by_id)
        self.assertNotIn("origin-slice-export", by_id)
        self.assertEqual(by_id["slice-time-multiplier"].value, 1)
        self.assertEqual(by_id["slice-time-multiplier"].type, "number")
        self.assertEqual(by_id["fit-slice-exponential"].value, [])
        image_export_menus = [
            component for component in components
            if isinstance(component, html.Details)
            and "image-export-options" in (getattr(component, "className", "") or "")
        ]
        self.assertEqual(len(image_export_menus), 1)
        menu_ids = [
            {
                component.id for component in _components(menu)
                if getattr(component, "id", None)
            }
            for menu in image_export_menus
        ]
        self.assertIn({
            "include-plot-title", "include-plot-legend",
            "origin-epsilon-export", "show-epsilon-export-grid",
            "zero-epsilon-export-y", "include-export-axis-labels",
        }, menu_ids)
        for message in ("epsilon-image-message",):
            state_ids = {state["id"] for state in self.app.callback_map[message + ".children"]["state"]}
            self.assertTrue(menu_ids[0].issubset(state_ids))
            self.assertIn("csv-export-kind", state_ids)
        self.assertFalse(any(
            {"minimal-spectrum-colors", "fit-slice-exponential"} & ids
            for ids in menu_ids
        ))
        self.assertTrue(all(not getattr(menu, "open", False) for menu in image_export_menus))
        menu_rows = [
            component for component in components
            if getattr(component, "className", None) == "plot-option-menus"
        ]
        self.assertEqual(len(menu_rows), 1)
        aligned_menu_labels = [
            str(component.children)
            for component in _components(menu_rows[0])
            if isinstance(component, html.Summary)
        ]
        self.assertEqual(
            aligned_menu_labels,
            ["Legend options", "Axis names"],
        )
        self.assertTrue(any(
            getattr(component, "className", None) == "info-popup"
            and "Existing timestamps are multiplied by this number" in str(
                component.to_plotly_json()
            )
            for component in components
        ))
        self.assertFalse(by_id["wavelength-slice-panel"].open)
        self.assertNotIn("responsive", by_id["epsilon-plot"].config)
        for component_id in (
            "save-epsilon-png", "save-epsilon-svg", "csv-export-kind", "save-session",
        ):
            self.assertIn(component_id, by_id)
        for removed in ("save-slice-csv", "export-isomer", "export-nmr", "save-folder", "choose-save-folder",
                        "save-filename", "session-filename", "nmr-reactant-filename", "nmr-product-filename",
                        "nmr-export-raw", "save-slice-png", "save-nmr-png"):
            self.assertNotIn(removed, by_id)
        self.assertIn("epsilon-image-message.children", self.app.callback_map)
        self.assertNotIn("nmr-image-message.children", self.app.callback_map)

    @unittest.skipUnless(shutil.which("node"), "Node is needed to exercise browser export code")
    def test_image_export_can_hide_title_legend_and_axis_labels_without_changing_preview(self):
        script = self.app.index_string.split("window.autoqyDownloadPlot =", 1)[1].split("window.autoqySaveText =", 1)[0]
        program = """
const assert = require('node:assert/strict');
let captured;
global.document = {
  querySelector: () => ({_fullLayout: {width: 1000, height: 600}}),
  createElement: () => ({style: {}, remove() {}, click() {}}),
  body: {appendChild() {}}
};
global.window = {Plotly: {
  newPlot: (holder, data, layout) => {captured = {data, layout}; return Promise.resolve();},
  toImage: () => Promise.resolve('data:image/png;base64,'), purge() {}
}};
""" + "window.autoqyDownloadPlot =" + script + """
(async () => {
  const figure = {data: [{name: 'Metastable ε', x: [400], y: [100], line: {color: '#d67b36'}}],
    layout: {title: {text: 'Reconstructed ε'}, xaxis: {title: {text: 'Wavelength'}},
      yaxis: {title: {text: 'ε'}}, showlegend: true}};
  const original = JSON.stringify(figure);
  const message = await window.autoqyDownloadPlot('nmr-plot', figure, 'png', [], [], ['on'], [], ['on'], []);
  assert.equal(message, 'PNG download started.');
  assert.equal(captured.layout.title, null);
  assert.equal(captured.layout.showlegend, false);
  assert.equal(captured.layout.xaxis.title.text, '');
  assert.equal(captured.layout.yaxis.title.text, '');
  assert.equal(captured.data.length, 1);
  assert.equal(JSON.stringify(figure), original);
  await window.autoqyDownloadPlot('nmr-plot', figure, 'svg', ['on'], ['on'], [], [], [], ['on']);
  assert.equal(captured.layout.title.text, 'Reconstructed ε');
  assert.equal(captured.layout.showlegend, true);
  assert.equal(captured.layout.yaxis.title.text, 'ε');
  let nativeSave;
  window.autoqySaveImage = async (filename, blob) => {nativeSave = {filename, blob}; return 'Saved in shared folder.';};
  window.fetch = async () => ({blob: async () => ({type: 'image/png'})});
  window.showSaveFilePicker = () => {throw new Error('Unexpected separate browser save dialog');};
  const nativeMessage = await window.autoqyDownloadPlot('nmr-plot', figure, 'png', [], [], [], [], [], []);
  assert.equal(nativeMessage, 'Saved in shared folder.');
  assert.equal(nativeSave.filename, 'reconstructed.png');
  assert.equal(nativeSave.blob.type, 'image/png');
})().catch(error => {console.error(error); process.exitCode = 1;});
"""
        result = subprocess.run([shutil.which("node"), "-"], input=program,
                                text=True, encoding="utf-8", capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("window.showSaveFilePicker", self.app.index_string)
        self.assertIn("exportLayout.showlegend = keepLegend", self.app.index_string)
        self.assertIn("exportLayout.title = null", self.app.index_string)
        self.assertIn("width = 1300", self.app.index_string)
        self.assertIn("height = 1000", self.app.index_string)
        self.assertIn("axis.ticks = 'outside'", self.app.index_string)
        self.assertIn("axis.mirror = true", self.app.index_string)
        self.assertIn("axis.linewidth = 3", self.app.index_string)
        self.assertIn("Math.max(4.2, lineWidth * 2)", self.app.index_string)
        self.assertIn("size: 40", self.app.index_string)
        self.assertIn("orientation: 'v'", self.app.index_string)
        self.assertIn("x: 0.98", self.app.index_string)
        self.assertIn("useOriginStyle ? 1.5 : 2", self.app.index_string)
        self.assertIn("showgrid: includeGrid", self.app.index_string)
        self.assertIn("axis.rangemode = 'nonnegative'", self.app.index_string)

    def test_main_plot_uses_colorblind_palette_by_default(self):
        self.assertEqual(_colors(), list(COLORBLIND_TRACE_PALETTE))
        self.assertEqual(_colors()[0], "#0072B2")
        self.assertEqual(_spectrum_colors(8)[7], _colors()[0])

    def test_legend_edit_toggle_is_in_legend_options_menu(self):
        menu = next(item for item in _components(self.app.layout)
                    if isinstance(item, html.Details) and isinstance(item.children[0], html.Summary)
                    and item.children[0].children == "Legend options")
        self.assertIn("edit-legend", [getattr(item, "id", None) for item in _components(menu)])

    def test_preview_exports_the_current_legend_names(self):
        data = _pack(SpectralDataset(np.array([400., 401.]), np.array([0.]),
                                    np.array([[0.2], [0.4]]), "csv"), ["original"], ["source.csv"])
        preview = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                       if "epsilon-plot.figure" in key)
        args = (data, 400, 401, [], None, None, "off", 11, 2, [], None,
                [0.001], [1], [["on"]], ["Edited name"], [], [],
                None, None, None, [{"type": "legend-spectrum", "index": 0}],
                [{"type": "legend-name", "index": 0}])
        result = preview(*args)
        self.assertEqual(result[4]["labels"], ["Edited name"])
        self.assertEqual(result[5]["legend_names"], ["Edited name"])
        for epsilon in (result[4], None):
            _, text = _export_csv_payload(epsilon, result[5])
            self.assertIn("Edited name", pd.read_csv(StringIO(text)).columns[1])

    def test_sheet_tracks_legend_clicks_without_losing_previous_visibility(self):
        callback = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                        if "displayed-data-sheet.columns" in key)
        figure = {"data": [{"x": [400], "y": [1], "name": "First"},
                           {"x": [400], "y": [2], "name": "Last"}]}
        token = context_value.set(AttributeDict(triggered_inputs=[{"prop_id": "epsilon-plot.restyleData"}]))
        try:
            result = callback("main", 0, 25, figure, None, None,
                              [{"visible": ["legendonly"]}, [0]], None, None, {})
            self.assertEqual(result[0][1]["name"], "Last")
            result = callback("main", 0, 25, figure, None, None,
                              [{"visible": ["legendonly"]}, [1]], None, None, result[-1])
            self.assertEqual(result[0], [])
        finally:
            context_value.reset(token)

    def test_sheet_callback_returns_all_rows_without_pagination(self):
        callback = next(value["callback"].__wrapped__ for key, value in self.app.callback_map.items()
                        if "displayed-data-sheet.columns" in key)
        figure = {"data": [{"x": list(range(305)), "y": list(range(305)), "name": "All rows"}]}
        token = context_value.set(AttributeDict(triggered_inputs=[{"prop_id": "epsilon-plot.figure"}]))
        try:
            result = callback("main", 3, 25, figure, None, None, None, None, None, {})
            self.assertEqual(len(result[1]), 305)
            self.assertEqual(result[1][-1]["x"], 304)
            self.assertEqual(result[4], 0)
        finally:
            context_value.reset(token)

    def test_sheet_has_whole_row_column_and_cell_copy_controls(self):
        by_id = {item.id: item for item in _components(self.app.layout)
                 if isinstance(getattr(item, "id", None), str)}
        sheet = by_id["displayed-data-sheet"]
        self.assertEqual(sheet.row_selectable, "multi")
        self.assertEqual(sheet.column_selectable, "multi")
        self.assertTrue(sheet.cell_selectable)
        self.assertTrue(sheet.include_headers_on_copy_paste)
        self.assertFalse(sheet.editable)
        self.assertIn("copy-sheet-selection", by_id)
        self.assertIn("copy-sheet-all", by_id)

    def test_calibrated_request_keeps_explicit_zero_nmr_error(self):
        callback = self.app.callback_map["isomer-request-store.data"]["callback"].__wrapped__
        token = context_value.set(AttributeDict(triggered_inputs=[{"prop_id": "reconstruct-isomer.n_clicks"}]))
        try:
            request = callback(1, {"labels": ["initial", "PSS"]}, 0, 1, 77,
                               None, 0, "sd", {"reference": {"filename": "reference.csv"}})
            self.assertEqual(request["error_percent"], 0)
            self.assertEqual(request["reference"], 0)
            self.assertEqual(request["mixture"], 1)
            self.assertEqual(request["mode"], "calibrated")
        finally:
            context_value.reset(token)

    def test_legend_can_show_only_selected_spectrum_traces(self):
        import plotly.graph_objects as go

        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 1.0, 2.0]),
            np.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]]),
            source_format="csv",
        )
        figure = _absorbance_figure(
            go, dataset, dataset.absorbance, dataset.absorbance,
            ["one", "two", "three"], "off",
            legend_visibility=[True, False, True],
        )
        self.assertEqual(
            [trace.showlegend for trace in figure.data], [True, False, True]
        )
        self.assertTrue(figure.layout.showlegend)

    def test_original_spectra_are_hidden_by_default_and_can_be_enabled(self):
        import plotly.graph_objects as go

        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 1.0]),
            np.array([[1.0, 2.0], [1.5, 2.5]]), source_format="csv",
        )
        processed = dataset.absorbance - 0.25
        labels = ["one", "two"]
        hidden = _absorbance_figure(
            go, dataset, dataset.absorbance, processed, labels, "off"
        )
        visible = _absorbance_figure(
            go, dataset, dataset.absorbance, processed, labels, "off",
            show_original=True,
        )
        self.assertEqual(len(hidden.data), 2)
        self.assertEqual(len(visible.data), 4)
        self.assertFalse(any(
            trace.name == "Uploaded absorbance" for trace in hidden.data
        ))
        self.assertEqual(
            sum(trace.name == "Uploaded absorbance" for trace in visible.data), 2
        )

    def test_slice_trace_remains_available_to_the_saved_image_legend(self):
        import plotly.graph_objects as go

        figure = _wavelength_slice_figure(
            go, np.array([0.0, 1.0]), np.array([1.0, 0.5]), 450.0
        )
        self.assertFalse(figure.layout.showlegend)
        self.assertTrue(figure.data[0].showlegend)

    def test_loaded_spectrum_rows_include_editable_legend_names(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 1.0]),
            np.array([[1.0, 2.0], [1.5, 2.5]]), source_format="csv",
        )
        rows = _loaded_spectrum_rows(
            html, dcc, _pack(dataset, ["one", "two"], ["one.csv", "two.csv"]),
            edit_legend=True,
        )
        components = list(_components(rows))
        legend_inputs = [
            component for component in components
            if isinstance(getattr(component, "id", None), dict)
            and component.id.get("type") == "legend-name"
        ]
        self.assertEqual([component.value for component in legend_inputs], ["one", "two"])

    def test_large_loaded_spectrum_manager_only_builds_first_and_last_controls(self):
        count = 571
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.arange(count, dtype=float),
            np.vstack((np.arange(count), np.arange(count) + 0.5)),
            source_format="csv",
        )
        labels = [f"spectrum {index}" for index in range(count)]
        rows = _loaded_spectrum_rows(
            html, dcc, _pack(dataset, labels, ["series.csv"]), edit_legend=True,
        )
        component_ids = [
            component.id for component in _components(rows)
            if isinstance(getattr(component, "id", None), dict)
        ]
        self.assertEqual(len(rows), 2)
        self.assertEqual(
            [item["index"] for item in component_ids if item["type"] == "legend-name"],
            [0, count - 1],
        )
        self.assertEqual(
            [item["index"] for item in component_ids
             if item["type"] == "legend-spectrum"],
            [0, count - 1],
        )
        self.assertFalse(any(
            item["type"] in {"move-spectrum-up", "move-spectrum-down", "remove-spectrum"}
            for item in component_ids
        ))

    def test_large_datasets_do_not_create_beer_lambert_inputs(self):
        labels = [f"spectrum {index}" for index in range(571)]
        cards = _parameter_cards(html, dcc, labels)
        serialized = str(cards[0].to_plotly_json())
        self.assertEqual(len(cards), 1)
        self.assertIn("controls are hidden", serialized)
        control_types = {
            component.id.get("type")
            for component in _components(cards)
            if isinstance(getattr(component, "id", None), dict)
        }
        self.assertNotIn("direct-concentration", control_types)
        self.assertNotIn("path-length", control_types)

    def test_sparse_legend_controls_merge_by_component_index(self):
        count = 571
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.arange(count, dtype=float),
            np.vstack((np.arange(count), np.arange(count) + 0.5)),
            source_format="csv",
        )
        labels = [f"spectrum {index}" for index in range(count)]
        packed = _pack(dataset, labels, ["series.csv"])
        packed["legend_visibility"] = [False] * count
        packed["legend_names"] = [f"stored {index}" for index in range(count)]
        updated = _with_spectrum_state(
            packed,
            legend_values=[["on"], []],
            legend_names=["Initial", "Final"],
            legend_value_ids=[
                {"type": "legend-spectrum", "index": 0},
                {"type": "legend-spectrum", "index": count - 1},
            ],
            legend_name_ids=[
                {"type": "legend-name", "index": 0},
                {"type": "legend-name", "index": count - 1},
            ],
        )
        self.assertTrue(updated["legend_visibility"][0])
        self.assertFalse(updated["legend_visibility"][1])
        self.assertFalse(updated["legend_visibility"][-1])
        self.assertEqual(updated["legend_names"][0], "Initial")
        self.assertEqual(updated["legend_names"][1], "stored 1")
        self.assertEqual(updated["legend_names"][-1], "Final")

    def test_large_absorbance_plot_has_sixty_spectra_and_two_legend_entries(self):
        import plotly.graph_objects as go

        count = 571
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.arange(count, dtype=float),
            np.vstack((np.arange(count), np.arange(count) + 0.5)),
            source_format="csv",
        )
        labels = [f"spectrum {index}" for index in range(count)]
        figure = _absorbance_figure(
            go, dataset, dataset.absorbance, dataset.absorbance,
            labels, "off", legend_visibility=[True] * count,
        )
        self.assertEqual(len(figure.data), MAX_INTERACTIVE_SPECTRA)
        self.assertEqual(
            [trace.name for trace in figure.data if trace.showlegend],
            [labels[0], labels[-1]],
        )
        self.assertTrue(figure.layout.showlegend)

    def test_large_epsilon_plot_has_sixty_spectra_and_two_spectrum_legends(self):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        count = 571
        absorbance = np.vstack((np.arange(count), np.arange(count) + 0.5))
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.arange(count, dtype=float),
            absorbance, source_format="csv",
        )
        result = EpsilonResult(
            wavelengths=dataset.wavelengths,
            absorbance=absorbance,
            concentrations_m=np.ones(count),
            path_lengths_cm=np.ones(count),
            individual=absorbance,
            mean=np.mean(absorbance, axis=1),
            standard_deviation=np.zeros(2),
            standard_error=np.zeros(2),
        )
        labels = [f"spectrum {index}" for index in range(count)]
        figure = _epsilon_figure(
            go, make_subplots, dataset, absorbance, result, labels, "off",
            legend_visibility=[True] * count,
        )
        spectrum_traces = [
            trace for trace in figure.data
            if str(trace.legendgroup).startswith("spectrum-")
        ]
        self.assertEqual(len(spectrum_traces), 2 * MAX_INTERACTIVE_SPECTRA)
        self.assertEqual(len({trace.legendgroup for trace in spectrum_traces}),
                         MAX_INTERACTIVE_SPECTRA)
        self.assertEqual(
            [trace.name for trace in spectrum_traces if trace.showlegend],
            [labels[0], labels[-1]],
        )

    @unittest.skipUnless(shutil.which("node"), "Node is needed to exercise figure selection")
    def test_two_figure_buttons_use_the_csv_content_selection(self):
        script = next(text for text in self.app._inline_scripts if "function(pngClicks, svgClicks, kind," in text)
        body = script.split("function(pngClicks, svgClicks, kind,", 1)[1].split("\n        }", 1)[0]
        function = "function(pngClicks, svgClicks, kind," + body + "\n        }"
        program = """
const assert = require('node:assert/strict');
let exported;
global.window = {
  dash_clientside: {callback_context: {triggered: [{prop_id: 'save-epsilon-png.n_clicks'}]}},
  autoqyDownloadPlot: (id, figure, format) => {exported = {id, figure, format}; return 'ok';}
};
const save = (""" + function + """ );
const main = {data: [{name: 'Measured'}, {name: 'Reconstructed pure metastable isomer'}],
  layout: {xaxis: {}, yaxis: {}, xaxis2: {title: {text: 'Wavelength'}}, yaxis2: {}}};
const slice = {data: [{name: 'Slice'}]}, epsilon = {data: [{name: 'Product ε'}]};
for (const [kind, request, expected] of [['processed', null, 'epsilon-plot'],
    ['slice', null, 'wavelength-slice-plot'], ['isomer', {mode: 'calibrated'}, 'nmr-plot'],
    ['nmr', null, 'nmr-plot']]) {
  assert.equal(save(1, 0, kind, main, slice, epsilon, request), 'ok');
  assert.equal(exported.id, expected);
}
save(1, 0, 'isomer', main, slice, epsilon, {});
assert.equal(exported.figure.data.length, 1);
assert.equal(exported.figure.layout.yaxis2, undefined);
assert.equal(main.data.length, 2);
window.dash_clientside.callback_context.triggered = [{prop_id: 'save-epsilon-svg.n_clicks'}];
save(0, 1, 'slice', main, slice, epsilon, null);
assert.equal(exported.format, 'svg');
"""
        result = subprocess.run([shutil.which("node"), "-"], input=program,
                                text=True, encoding="utf-8", capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_main_minimal_palette_includes_mean_and_isomer_without_purple(self):
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from autoqy_core.tools.smoother_gui import _overlay_isomer
        data = SpectralDataset(np.array([400., 401.]), np.array([0., 1.]),
                               np.array([[.2, .5], [.8, .4]]), "csv")
        result = calculate_epsilon_statistics(data.wavelengths, data.absorbance, [.001, .001], [1, 1])
        figure = _epsilon_figure(go, make_subplots, data, data.absorbance, result, ["stable", "PSS"], "off", minimal_colors=True)
        request = {"reference": 0, "mixture": 1, "percent": 80, "error_percent": 1}
        _overlay_isomer(go, figure, data, data.absorbance, request, ["on"], epsilon=True, minimal_colors=True)
        traces = {trace.name: trace for trace in figure.data}
        self.assertEqual(traces["Mean ε"].line.color, "#d67b36")
        self.assertEqual(traces["Reconstructed pure metastable isomer"].line.color, "#000000")
        self.assertEqual(traces["Reconstructed pure metastable isomer"].line.dash, "dot")
        self.assertFalse(any(trace.line.color == "#8a5a9b" or "138,90,155" in (trace.fillcolor or "")
                             for trace in figure.data))

    def test_minimal_palette_marks_initial_intermediate_and_final_spectra(self):
        self.assertEqual(_spectrum_colors(1, True), ["#2d6f8e"])
        self.assertEqual(
            _spectrum_colors(4, True),
            [
                "#2d6f8e",
                "rgba(108,114,128,0.38)",
                "rgba(108,114,128,0.38)",
                "#d67b36",
            ],
        )
        self.assertEqual(_spectrum_colors(7), [
            COLORBLIND_TRACE_PALETTE[index % len(COLORBLIND_TRACE_PALETTE)]
            for index in range(7)
        ])

    def test_single_time_series_keeps_coordinates_and_interpolates_a_slice(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]),
            np.array([0.0, 2.5, 5.0]),
            np.array([[0.0, 1.0, 2.0], [10.0, 11.0, 12.0]]),
            source_format="spectragryph",
        )
        combined, _, _ = _combine_loaded([(dataset, "series.dat")])
        np.testing.assert_array_equal(combined.coordinates, dataset.coordinates)
        selected, coordinates, values = _wavelength_slice(combined, 405.0)
        self.assertEqual(selected, 405.0)
        np.testing.assert_array_equal(coordinates, dataset.coordinates)
        np.testing.assert_allclose(values, [5.0, 6.0, 7.0])

    def test_multiple_files_are_naturally_sorted_by_filename(self):
        def spectrum(value):
            return SpectralDataset(
                np.array([400.0, 410.0]), np.array([0.0]),
                np.array([[value], [value + 0.5]]), source_format="csv",
            )

        combined, labels, _ = _combine_loaded([
            (spectrum(10.0), "sample10.csv"),
            (spectrum(2.0), "sample2.csv"),
            (spectrum(1.0), "Sample1.csv"),
            (spectrum(3.0), "sampleA.csv"),
            (spectrum(4.0), "sampleB.csv"),
        ])
        self.assertEqual(
            labels, ["Sample1", "sample2", "sample10", "sampleA", "sampleB"]
        )
        np.testing.assert_array_equal(
            combined.absorbance[0], [1.0, 2.0, 10.0, 3.0, 4.0]
        )

    def test_repeated_single_spectrum_drops_ignore_elapsed_time(self):
        first = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0]),
            np.array([[1.0], [2.0]]), source_format="csv",
        )
        second = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0]),
            np.array([[3.0], [4.0]]), source_format="csv",
        )
        combined = _append_packed(
            _pack(first, ["spectrum"], ["spectrum.csv"]),
            _pack(second, ["spectrum"], ["spectrum.csv"]),
            elapsed_seconds=2.5,
        )
        self.assertEqual(combined["labels"], ["spectrum", "spectrum (2)"])
        self.assertEqual(combined["coordinates"], [0.0, 1.0])
        np.testing.assert_array_equal(
            np.asarray(combined["absorbance"]),
            np.array([[1.0, 3.0], [2.0, 4.0]]),
        )

    def test_new_spectra_receive_the_same_baseline_and_smoothing(self):
        wavelengths = np.arange(400.0, 421.0)
        first_values = 0.4 + np.exp(-((wavelengths - 410.0) / 3.0) ** 2)
        second_values = 0.8 + 0.6 * np.exp(-((wavelengths - 412.0) / 4.0) ** 2)
        first = SpectralDataset(
            wavelengths, np.array([0.0]), first_values[:, None], source_format="csv"
        )
        second = SpectralDataset(
            wavelengths, np.array([0.0]), second_values[:, None], source_format="csv"
        )
        first_packed = _pack(first, ["first"], ["first.csv"])
        second_packed = _pack(second, ["second"], ["second.csv"])
        combined = _append_packed(first_packed, second_packed, elapsed_seconds=30.0)

        settings = (400.0, 420.0, ["on"], 400.0, 403.0, "savgol", 5.0, 2)
        _, _, combined_processed, _ = _prepare_processing(combined, *settings)
        _, _, first_processed, _ = _prepare_processing(first_packed, *settings)
        _, _, second_processed, _ = _prepare_processing(second_packed, *settings)

        np.testing.assert_allclose(combined_processed[:, 0], first_processed[:, 0])
        np.testing.assert_allclose(combined_processed[:, 1], second_processed[:, 0])
        self.assertFalse(np.allclose(combined_processed[:, 1], second_values))

    def test_dataset_updates_preserve_selected_wavelength_range(self):
        callback = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "wavelength-low.value" in key and "slice-wavelength.max" in key
        )
        dataset = SpectralDataset(
            np.arange(400.0, 411.0), np.array([0.0]),
            np.arange(11.0)[:, None], source_format="csv",
        )
        result = callback(_pack(dataset, ["one"], ["one.csv"]), 402.0, 408.0, 405.0)
        self.assertEqual(result[0:2], (402.0, 408.0))
        self.assertEqual(result[9], 405.0)

        narrower = SpectralDataset(
            np.arange(404.0, 407.0), np.array([0.0]),
            np.arange(3.0)[:, None], source_format="csv",
        )
        clipped = callback(
            _pack(narrower, ["one"], ["one.csv"]), 402.0, 408.0, 405.0
        )
        self.assertEqual(clipped[0:2], (404.0, 406.0))

    def test_selected_spectra_can_be_removed_without_retiming_the_rest(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 2.5, 8.0]),
            np.array([[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]]),
            source_format="csv",
        )
        packed = _pack(
            dataset, ["first", "middle", "last"], ["series.csv"],
            concentrations=[1.0, 2.0, 3.0], path_lengths=[1.0, 1.0, 2.0],
        )
        reduced = _remove_packed(
            packed, [1], concentrations=[10.0, 20.0, 30.0],
            path_lengths=[0.5, 1.0, 1.5],
            legend_values=[["on"], [], ["on"]],
            legend_names=["Start", "Middle", "Finish"],
        )
        self.assertEqual(reduced["labels"], ["first", "last"])
        self.assertEqual(reduced["coordinates"], [0.0, 8.0])
        self.assertEqual(reduced["concentrations"], [10.0, 30.0])
        self.assertEqual(reduced["path_lengths"], [0.5, 1.5])
        self.assertEqual(reduced["legend_visibility"], [True, True])
        self.assertEqual(reduced["legend_names"], ["Start", "Finish"])
        np.testing.assert_array_equal(
            np.asarray(reduced["absorbance"]),
            np.array([[1.0, 5.0], [2.0, 6.0]]),
        )
        self.assertIsNone(_remove_packed(packed, [0, 1, 2]))

    def test_loaded_spectra_can_be_reordered_with_their_settings(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0, 4.0, 9.0]),
            np.array([[1.0, 2.0, 3.0], [11.0, 12.0, 13.0]]),
            source_format="csv",
        )
        packed = _pack(
            dataset, ["one", "two", "three"],
            ["one.csv", "two.csv", "three.csv"],
        )
        reordered = _reorder_packed(
            packed, [1, 0, 2], concentrations=[1.0, 2.0, 3.0],
            path_lengths=[0.5, 1.0, 1.5],
            legend_values=[["on"], [], ["on"]],
            legend_names=["Initial", "Intermediate", "Final"],
        )
        self.assertEqual(reordered["labels"], ["two", "one", "three"])
        self.assertEqual(reordered["filenames"], ["two.csv", "one.csv", "three.csv"])
        self.assertEqual(reordered["coordinates"], [0.0, 4.0, 9.0])
        self.assertEqual(reordered["concentrations"], [2.0, 1.0, 3.0])
        self.assertEqual(reordered["legend_visibility"], [False, True, True])
        self.assertEqual(
            reordered["legend_names"], ["Intermediate", "Initial", "Final"]
        )
        np.testing.assert_array_equal(
            np.asarray(reordered["absorbance"])[0], [2.0, 1.0, 3.0]
        )

    def test_wavelength_slice_rejects_values_outside_the_processed_range(self):
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), np.array([0.0]),
            np.array([[0.0], [1.0]]), source_format="csv",
        )
        with self.assertRaisesRegex(ValueError, "between 400 and 410"):
            _wavelength_slice(dataset, 420.0)

    def test_slice_callback_builds_plot_export_data_and_custom_axes(self):
        callback = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "wavelength-slice-store.data" in key
        )
        dataset = SpectralDataset(
            np.array([400.0, 410.0]),
            np.array([0.0, 5.0]),
            np.array([[0.0, 2.0], [10.0, 12.0]]),
            source_format="spectragryph",
        )
        figure, data, message, fit_message, fit_class = callback(
            _pack(dataset, ["start", "end"], ["series.dat"]),
            405.0, "Elapsed time", "Optical density", 1, [],
        )
        self.assertEqual(data["values"], [5.0, 7.0])
        self.assertEqual(data["x_label"], "Elapsed time (s)")
        self.assertEqual(data["y_label"], "Optical density at 405 nm")
        self.assertEqual(figure.layout.xaxis.title.text, "Elapsed time (s)")
        self.assertEqual(figure.layout.yaxis.title.text, "Optical density")
        self.assertEqual(message, "")
        self.assertEqual(fit_message, "")
        self.assertEqual(fit_class, "slice-fit-result")

    def test_exponential_slice_fit_reports_lifetime_error_and_duration_flag(self):
        coordinates = np.linspace(0.0, 20.0, 21)
        values = 0.2 + 1.5 * np.exp(-coordinates / 4.0)
        values += np.array([0.002, -0.001, 0.001, 0.0, -0.002, 0.001, 0.0,
                            0.001, -0.001, 0.0, 0.001, 0.0, -0.001, 0.001,
                            0.0, -0.001, 0.0, 0.001, 0.0, -0.001, 0.0])
        fit = _fit_exponential_decay(coordinates, values)
        self.assertAlmostEqual(fit["lifetime"], 4.0, delta=0.1)
        self.assertGreaterEqual(fit["lifetime_error"], 0.0)
        self.assertLess(fit["lifetime_error"], 0.1)
        self.assertGreater(fit["duration"], fit["lifetime"])
        import plotly.graph_objects as go
        figure = _wavelength_slice_figure(go, coordinates, values, 450, fit=fit)
        self.assertEqual(figure.data[-1].line.dash, "dot")

    def test_slice_callback_scales_time_before_exponential_fit(self):
        callback = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "wavelength-slice-store.data" in key
        )
        coordinates = np.linspace(0.0, 20.0, 21)
        decay = 0.1 + 0.9 * np.exp(-coordinates / 5.0)
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), coordinates,
            np.vstack((decay, decay)), source_format="spectragryph",
        )
        figure, data, _, fit_message, fit_class = callback(
            _pack(dataset, [f"point {index}" for index in range(len(coordinates))],
                  ["series.dat"]),
            405.0, "Time", "Absorbance", 30, ["on"],
        )
        self.assertEqual(len(figure.data), 2)
        self.assertEqual(data["time_multiplier"], 30)
        self.assertEqual(data["coordinates"][-1], 600)
        self.assertEqual(data["time_unit"], "s")
        self.assertEqual(data["x_label"], "Time (s)")
        self.assertAlmostEqual(data["lifetime"], 150, delta=0.1)
        self.assertIn("fit_values", data)
        self.assertIn("Lifetime τ", fit_message)
        self.assertIn("extends beyond one fitted lifetime", fit_message)
        self.assertIn("status-ok", fit_class)

    def test_slice_fit_warns_when_recording_is_shorter_than_lifetime(self):
        callback = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "wavelength-slice-store.data" in key
        )
        coordinates = np.linspace(0.0, 20.0, 21)
        slow_decay = 0.1 + 0.9 * np.exp(-coordinates / 50.0)
        dataset = SpectralDataset(
            np.array([400.0, 410.0]), coordinates,
            np.vstack((slow_decay, slow_decay)), source_format="spectragryph",
        )
        _, _, _, fit_message, fit_class = callback(
            _pack(dataset, [f"point {index}" for index in range(len(coordinates))],
                  ["series.dat"]),
            405.0, "Time", "Absorbance", 1, ["on"],
        )
        self.assertIn("shorter than one fitted lifetime", fit_message)
        self.assertIn("status-warning", fit_class)

    def test_preview_and_save_work_without_concentrations(self):
        callbacks = self.app.callback_map
        preview = next(
            value["callback"].__wrapped__
            for key, value in callbacks.items()
            if "processed-store.data" in key
        )
        save = next(
            value["callback"].__wrapped__
            for key, value in callbacks.items()
            if "epsilon-save-message" in key
        )
        dataset = SpectralDataset(
            np.array([400.0, 401.0, 402.0]),
            np.array([0.0]),
            np.array([[0.1], [-0.2], [0.3]]),
            source_format="csv",
        )
        packed = _pack(dataset, ["irradiation"], ["irradiation.csv"])
        preview_result = preview(
            packed, 400, 402, [], None, None, "off", 5, 3,
            [], 1, [None], [1.0], [[]], ["Custom legend"], [],
            [], "Custom wavelength", "Custom OD", "Custom epsilon",
            [{"type": "legend-spectrum", "index": 0}],
            [{"type": "legend-name", "index": 0}],
        )
        self.assertIsNone(preview_result[4])
        self.assertIsNotNone(preview_result[5])
        self.assertFalse(preview_result[6])
        self.assertFalse(preview_result[0].layout.showlegend)
        self.assertEqual(preview_result[0].data[0].name, "Custom legend")
        self.assertEqual(preview_result[0].layout.xaxis.title.text, "Custom wavelength")
        self.assertEqual(preview_result[0].layout.yaxis.title.text, "Custom OD")
        self.assertEqual(tuple(preview_result[0].layout.xaxis.range), (400.0, 402.0))
        self.assertFalse(preview_result[0].layout.xaxis.autorange)
        self.assertEqual(
            preview_result[0].layout.uirevision, "wavelength-range:400:402"
        )
        self.assertEqual(preview_result[5]["labels"], ["irradiation"])

        with TemporaryDirectory() as temporary:
            token = context_value.set(AttributeDict(
                triggered_inputs=[{"prop_id": "export-epsilon.n_clicks"}]
            ))
            try:
                with patch("autoqy_core.tools.smoother_gui._choose_save_file",
                           return_value=str(Path(temporary) / "my-treated-irradiation.csv")):
                    save_result = save(1, None, preview_result[5], ["on"])
            finally:
                context_value.reset(token)
            destination = Path(temporary) / "my-treated-irradiation.csv"
            self.assertTrue(destination.is_file())
            self.assertIn(str(destination), save_result[0])
            frame = pd.read_csv(destination)
            self.assertEqual(frame["Custom legend"].min(), 0.0)

    def test_epsilon_plot_locks_both_x_axes_to_selected_wavelengths(self):
        preview = next(
            value["callback"].__wrapped__
            for key, value in self.app.callback_map.items()
            if "processed-store.data" in key
        )
        dataset = SpectralDataset(
            np.arange(400.0, 411.0), np.array([0.0]),
            np.linspace(0.1, 0.5, 11)[:, None], source_format="csv",
        )
        result = preview(
            _pack(dataset, ["sample"], ["sample.csv"]),
            402.0, 408.0, [], None, None, "off", 5, 2,
            [], 1, [1e-5], [1.0], [["on"]], ["sample"], [],
            [], "Wavelength (nm)", "Absorbance", "ε (M⁻¹ cm⁻¹)",
            [{"type": "legend-spectrum", "index": 0}],
            [{"type": "legend-name", "index": 0}],
        )
        figure = result[0]
        self.assertEqual(tuple(figure.layout.xaxis.range), (402.0, 408.0))
        self.assertEqual(tuple(figure.layout.xaxis2.range), (402.0, 408.0))
        self.assertFalse(figure.layout.xaxis.autorange)
        self.assertFalse(figure.layout.xaxis2.autorange)


if __name__ == "__main__":
    unittest.main()
