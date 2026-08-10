from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import pandas as pd

from autoqy_core import (AnalysisConfig, format_value_uncertainty, load_config,
                         run_analysis, run_power_analysis)
from autoqy_core.cli import main
from autoqy_core.io import load_spectrum, load_timestamps
from autoqy_core.power import (REGION_NAMES, baseline_power_trace,
                               load_thorlabs_opm, process_power_trace)
from autoqy_core.power_web import (_initial_regions, _positions_from_relayout,
                                   _positions_to_regions)
from autoqy_core.spectra import fit_concentrations, fit_concentrations_regularized


ROOT = Path(__file__).parents[1]


class CoreTests(unittest.TestCase):
    def test_concentrations_are_recovered_from_synthetic_spectra(self):
        epsilon_r = np.array([1.0, 2.0, 3.0])
        epsilon_p = np.array([3.0, 1.0, 2.0])
        expected = np.array([[0.8, 0.2], [0.3, 0.7]]) * 0.001
        absorbance = (expected @ np.vstack((epsilon_r, epsilon_p))).T
        result = fit_concentrations(absorbance, np.arange(3), epsilon_r, epsilon_p)
        np.testing.assert_allclose(result.concentrations, expected)
        np.testing.assert_allclose(result.residuals, 0, atol=1e-15)

    def test_regularization_avoids_consecutive_boundary_clipping(self):
        wavelengths = np.arange(5.0)
        epsilon_r = np.array([2.0, 3.0, 4.0, 3.0, 2.0])
        epsilon_p = np.array([2.0, 2.5, 3.8, 4.0, 4.5])
        product_fraction = np.array([0, 0.015, 0.03, 0.05, 0.08, 0.12])
        mismatch = 0.04
        absorbance = np.column_stack([
            1e-5 * ((1 - fraction + mismatch) * epsilon_r
                    + (fraction - mismatch) * epsilon_p)
            for fraction in product_fraction
        ])
        legacy = fit_concentrations(absorbance, wavelengths, epsilon_r, epsilon_p)
        regularized = fit_concentrations_regularized(
            absorbance, wavelengths, epsilon_r, epsilon_p,
            np.arange(len(product_fraction)),
        )
        np.testing.assert_allclose(legacy.fractions[:3, 1], 0)
        self.assertGreater(regularized.fractions[1, 1], 0)
        self.assertTrue(np.all(np.diff(regularized.fractions[:, 1]) > 0))
        np.testing.assert_allclose(regularized.concentrations.sum(axis=1),
                                   regularized.concentrations[0].sum())

    def test_bundled_concentration_example_and_detailed_outputs(self):
        folder = ROOT / "ExampleData" / "Example-2_AB_455nm-100mA"
        config = load_config(folder / "analysis.json")
        with TemporaryDirectory() as temporary_directory:
            run = run_analysis(config, temporary_directory)
            result = run.result
            stem = Path(temporary_directory) / config.values["outputs"]["stem"]
            expected = folder / "Results_AB_455nm_first40_ODEMethod_Concentrations.txt"
            self.assertEqual(stem.with_suffix(".txt").read_text().splitlines(),
                             expected.read_text().splitlines())
            self.assertGreater(stem.with_suffix(".png").stat().st_size, 100_000)
            self.assertGreater(stem.with_suffix(".svg").stat().st_size, 100_000)
            summary = json.loads((Path(temporary_directory) /
                                  f"{stem.name}_results.json").read_text())
            self.assertEqual(summary["schema_version"], 2)
            self.assertEqual(summary["fit_method"], "concentrations")
            self.assertEqual(summary["experiment"]["volume_ml"], 3.02)
            self.assertNotIn("photostationary_state_percent", summary)
            last = summary["composition_at_last_timestamp_percent"]
            pss = summary["extrapolated_pss_percent"]
            self.assertAlmostEqual(last["time_s"], 1094.469)
            self.assertAlmostEqual(last["reactant"], 77.7141, places=3)
            self.assertAlmostEqual(pss["reactant"], 77.5634, places=3)
            self.assertEqual(len(run.files), 7)

            traces = pd.read_csv(Path(temporary_directory) / f"{stem.name}_traces.tsv", sep="\t")
            np.testing.assert_allclose(
                traces["product_concentration_residual_M"],
                traces["product_concentration_data_M"] - traces["product_concentration_fit_M"],
            )
            spectra = pd.read_csv(Path(temporary_directory) / f"{stem.name}_spectra.tsv", sep="\t")
            np.testing.assert_allclose(
                spectra["kinetic_fit_residual"],
                spectra["absorbance_measured"] - spectra["absorbance_kinetic_fit"],
            )
            np.testing.assert_allclose(
                result.absorbance,
                result.concentration_fit.reconstructed_absorbance.T
                + result.concentration_fit.residuals.T,
            )
        np.testing.assert_allclose(result.yield_fit.values * 100, [36.8, 49.9], atol=0.2)
        np.testing.assert_allclose(result.yield_errors * 100, [0.4, 0.6], atol=0.2)

    def test_emission_method_and_generic_spectral_input(self):
        folder = ROOT / "ExampleData" / "Example-1b_340nm_ManualPower_blcorrLED"
        config = _without_outputs(load_config(folder / "analysis.json"))
        result = run_analysis(config).result
        self.assertEqual(result.fit_method, "emission")
        np.testing.assert_allclose(result.yield_fit.values * 100, [17.2, 41.2], atol=0.3)
        wavelengths, epsilon = load_spectrum(
            folder / "cis_azobenzene_epsilons.csv",
            {"type": "generic_delimited", "delimiter": ",", "header": False,
             "wavelength_column": 0, "value_columns": [1]},
        )
        self.assertEqual(len(wavelengths), len(epsilon))
        self.assertAlmostEqual(wavelengths[0], 250.0)

    def test_new_methods_agree_on_azobenzene_concentration_trajectory(self):
        folder = ROOT / "ExampleData" / "Example-2_AB_455nm-100mA"
        config = _without_outputs(load_config(folder / "analysis.json"))
        results = {}
        for method in ("concentrations", "emission", "regularized_concentrations",
                       "ode_absorbance"):
            values = deepcopy(config.values)
            values["fit"]["method"] = method
            selected = AnalysisConfig(values, config.base_directory, config.source)
            results[method] = run_analysis(selected).result.yield_fit.concentrations

        reference = results["concentrations"]
        concentration_scale = reference[0].sum()
        for method in ("regularized_concentrations", "ode_absorbance"):
            normalized_rmse = (np.sqrt(np.mean((results[method] - reference) ** 2))
                               / concentration_scale)
            self.assertLess(normalized_rmse, 0.04)

        legacy_emission_rmse = (
            np.sqrt(np.mean((results["emission"] - reference) ** 2))
            / concentration_scale
        )
        self.assertGreater(legacy_emission_rmse, 0.05)

    def test_standalone_power_processing(self):
        config = (ROOT / "ExampleData" / "Example-2_AB_455nm-100mA" /
                  "Power" / "power_analysis.json")
        with TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "power.json"
            result = run_power_analysis(config, output)
            self.assertTrue(output.is_file())
            self.assertEqual(len(result.traces), 3)
            self.assertAlmostEqual(result.power_mw, 2.0666, places=3)
            self.assertAlmostEqual(result.repeatability_sd_mw, 0.0226, places=3)
            self.assertAlmostEqual(result.standard_error_mw, 0.0062, places=3)
            self.assertEqual(main(["power", str(config), "--output",
                                   str(Path(temporary_directory) / "cli-power.json")]), 0)

    def test_power_gui_helpers_use_the_core_baseline(self):
        folder = ROOT / "ExampleData" / "Example-2_AB_455nm-100mA" / "Power"
        values = load_thorlabs_opm(folder / "Power455nm_CH3CN_10pct_1.csv")
        regions = _initial_regions(len(values))
        self.assertEqual(tuple(regions), REGION_NAMES)
        positions = [boundary for name in REGION_NAMES for boundary in regions[name]]
        moved = _positions_from_relayout({"shapes[2].x0": positions[2] + 3}, regions)
        updated = _positions_to_regions(moved, len(values))
        self.assertEqual(updated["open_on"][0], regions["open_on"][0] + 3)

        configured = json.loads((folder / "power_analysis.json").read_text())
        selected = configured["measurements"][0]["regions"]
        detail = baseline_power_trace(values, selected, 3)
        trace = process_power_trace(values, selected, 3)
        open_values = detail.open_corrected_mw[slice(*detail.regions["open_on"])]
        cuvette_values = detail.cuvette_corrected_mw[slice(*detail.regions["cuvette_on"])]
        self.assertAlmostEqual(trace.open_beam_power_mw, np.mean(open_values))
        self.assertAlmostEqual(trace.cuvette_power_mw, np.mean(cuvette_values))

    def test_generic_timestamp_input(self):
        with TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "times.csv"
            path.write_text("measurement,irradiation_seconds\n1,0\n2,12.5\n", encoding="utf-8")
            timestamps = load_timestamps(
                path,
                {"type": "generic_delimited", "delimiter": ",", "header": True,
                 "time_column": "irradiation_seconds"},
            )
            np.testing.assert_allclose(timestamps, [0, 12.5])

    def test_uncertainty_controls_printed_precision(self):
        self.assertEqual(format_value_uncertainty(36.7737, 0.4408), ("36.8", "0.4"))
        self.assertEqual(format_value_uncertainty(1.234567, 0.000012),
                         ("1.234567", "0.000012"))
        self.assertEqual(format_value_uncertainty(1234, 120), ("1230", "120"))

    def test_cli_validates_example(self):
        config = ROOT / "ExampleData" / "Example-2_AB_455nm-100mA" / "analysis.json"
        self.assertEqual(main(["validate", str(config)]), 0)


def _without_outputs(config):
    values = deepcopy(config.values)
    values["outputs"].update(
        write_text=False, write_figures=False, write_json=False,
        write_config=False, write_detailed_data=False,
    )
    return AnalysisConfig(values, config.base_directory, config.source)


if __name__ == "__main__":
    unittest.main()
