"""Regression coverage for non-finite instrument readings."""

import base64
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from autoqy_core.power import (
    REGION_NAMES, baseline_power_trace, load_generic_power_csv,
    load_thorlabs_opm, load_thorlabs_opm_text, process_power_trace,
)


def opm_text(values):
    return "\n" * 14 + "Samples,Power (W)\n" + "\n".join(
        f"{i},{value}" for i, value in enumerate(values)
    )


class InvalidPowerTests(unittest.TestCase):
    def setUp(self):
        self.regions = dict(zip(REGION_NAMES, (
            [0, 4], [4, 8], [8, 12], [12, 16], [16, 20], [20, 24],
        )))
        self.values = 0.1 * np.arange(24) + 2
        self.values[4:8] += 10
        self.values[16:20] += 8

    def test_instrument_file_and_upload_keep_positions_and_finite_values(self):
        text = opm_text(["0.001", "Infinity", "0.003", "-inf", "NaN"])
        with TemporaryDirectory() as directory:
            path = Path(directory) / "power.csv"
            path.write_text(text)
            for load in (lambda: load_thorlabs_opm(path),
                         lambda: load_thorlabs_opm_text(text)):
                with self.assertWarnsRegex(UserWarning, "excluded 3"):
                    values = load()
                np.testing.assert_array_equal(values, [1, np.nan, 3, np.nan, np.nan])

    def test_generic_csv_excludes_infinities(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "power.csv"
            path.write_text("power_mw\n1\ninf\n3\n-inf\n")
            with self.assertWarnsRegex(UserWarning, "excluded 2"):
                values = load_generic_power_csv(path)
            np.testing.assert_array_equal(values, [1, np.nan, 3, np.nan])

    def test_all_invalid_and_nonnumeric_uploads_still_fail(self):
        with self.assertRaisesRegex(ValueError, "No finite"):
            load_thorlabs_opm_text(opm_text(["Infinity", "-Infinity", "NaN"]))
        with self.assertRaises(ValueError):
            load_thorlabs_opm_text(opm_text(["0.001", "broken"]))

    def test_finite_input_is_unchanged(self):
        result = process_power_trace(self.values, self.regions, 1)
        self.assertAlmostEqual(result.open_beam_power_mw, 10)
        self.assertAlmostEqual(result.cuvette_power_mw, 8)
        self.assertAlmostEqual(result.power_mw, 9)

    def test_baseline_and_signal_omit_invalid_values_at_original_indices(self):
        values = self.values.copy()
        values[[1, 5, 10, 13, 18, 22]] = [np.inf, np.inf, -np.inf, np.nan, np.nan, np.inf]
        baseline = baseline_power_trace(values, self.regions, 1)
        np.testing.assert_allclose(baseline.open_baseline_mw, .1 * np.arange(24) + 2)
        np.testing.assert_allclose(baseline.cuvette_baseline_mw, .1 * np.arange(24) + 2)
        self.assertTrue(np.isnan(baseline.open_corrected_mw[5]))
        self.assertTrue(np.isinf(values[1]))  # Caller data must remain intact.
        result = process_power_trace(values, self.regions, 1)
        self.assertAlmostEqual(result.power_mw, 9)
        self.assertTrue(np.isfinite(result.standard_error_mw))

    def test_uncertainties_use_only_surviving_signal_readings(self):
        values = np.zeros(24)
        values[4:8] = [8, np.inf, 10, 12]
        values[16:20] = [6, 8, -np.inf, 10]
        result = process_power_trace(values, self.regions, 0)
        self.assertAlmostEqual(result.power_mw, 9)
        self.assertAlmostEqual(result.repeatability_sd_mw, 2)
        self.assertAlmostEqual(result.standard_error_mw, np.sqrt(2 / 3))

    def test_insufficient_valid_signal_and_baseline_points_fail_clearly(self):
        for bad_indices, degree, message in (
            ([4, 5, 6], 1, "open_on"),
            (list(range(16, 20)), 1, "cuvette_on"),
            (list(range(4)), 1, "Each power baseline region"),
            ([0, 1, 2, 8, 9, 10], 3, "too few finite points"),
        ):
            with self.subTest(bad_indices=bad_indices):
                values = self.values.copy()
                values[bad_indices] = np.inf
                with self.assertRaisesRegex(ValueError, message):
                    process_power_trace(values, self.regions, degree)

    def test_gui_upload_survives_json_storage_and_can_calculate(self):
        try:
            from autoqy_core import power_web
        except ImportError:
            self.skipTest("Dash GUI dependencies are not installed")
        app = power_web.create_app()
        callback = next(item["callback"].__wrapped__
                        for key, item in app.callback_map.items()
                        if "power-store.data" in key)
        text = opm_text(["Infinity" if i == 5 else str(value / 1000)
                         for i, value in enumerate(self.values)])
        content = "data:text/csv;base64," + base64.b64encode(text.encode()).decode()
        context = next(cell.cell_contents for name, cell in
                       zip(callback.__code__.co_freevars, callback.__closure__)
                       if name == "ctx")
        with patch.object(type(context), "triggered_id", property(lambda _: "upload-power")):
            with self.assertWarns(UserWarning):
                store, message = callback([content], None, 0, 0, ["trace.csv"],
                                          {"datasets": []}, None, 1)
        self.assertIn("trace.csv: 1", message)
        dataset = json.loads(json.dumps(store, allow_nan=False))["datasets"][0]
        self.assertIsNone(dataset["values_mw"][5])
        dataset["regions"] = self.regions
        figure = power_web._raw_figure(dataset)
        self.assertEqual(len(figure.data[0].x), 24)
        with patch.object(type(context), "triggered_id", property(lambda _: "calculate-power")):
            updated, _ = callback(None, None, 1, 0, None, {"datasets": [dataset]},
                                  dataset["id"], 1)
        self.assertAlmostEqual(updated["datasets"][0]["result"]["power_mw"], 9)
        json.dumps(updated, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
