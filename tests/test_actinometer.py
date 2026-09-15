import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from autoqy_core.config import AnalysisConfig, ConfigError, validate_config
from autoqy_core.kinetics import (photon_flux_mol_s_to_power_mw,
                                  power_mw_to_photon_flux_mol_s)
from autoqy_core.output import (format_scientific_value_uncertainty,
                                format_value_uncertainty, result_summary)
from autoqy_core.runner import run_analysis
from autoqy_core.spectra import monochromatic_emission


EXAMPLE_DIRECTORY = (
    Path(__file__).parents[1]
    / "ExampleData"
    / "Example-4_395nm-EpsilonError"
    / "generic_inputs"
)


def _example_actinometer_config(output_directory):
    values = json.loads((EXAMPLE_DIRECTORY / "analysis.json").read_text(encoding="utf-8"))
    experiment = values["experiment"]
    wavelength = experiment["irradiation_wavelength_nm"]
    photon_flux = power_mw_to_photon_flux_mol_s(experiment.pop("power_mw"), wavelength)
    photon_flux_error = power_mw_to_photon_flux_mol_s(
        experiment.pop("power_error_mw"), wavelength
    )
    experiment.update({
        "irradiation_source": "chemical_actinometer",
        "photon_flux_mol_s": photon_flux,
        "photon_flux_error_mol_s": photon_flux_error,
    })
    values["inputs"].pop("led_emission")
    values["inputs"]["formats"].pop("led_emission")
    values["outputs"].update({
        "directory": str(output_directory),
        "write_text": False,
        "write_figures": False,
        "write_json": False,
        "write_config": False,
        "write_detailed_data": False,
    })
    return values


class ActinometerTests(unittest.TestCase):
    def test_uncertainty_rollover_removes_spurious_value_digits(self):
        self.assertEqual(format_value_uncertainty(4.821, 0.099), ("4.8", "0.1"))
        self.assertEqual(
            format_scientific_value_uncertainty(4.821e-9, 0.099e-9),
            ("4.8", "0.1", -9),
        )

    def test_power_conversion_round_trip_at_nominal_wavelength(self):
        photon_flux = power_mw_to_photon_flux_mol_s(1.46, 395)
        self.assertAlmostEqual(photon_flux, 4.820999270427747e-9, places=20)
        self.assertAlmostEqual(photon_flux_mol_s_to_power_mw(photon_flux, 395), 1.46)

    def test_monochromatic_profile_has_one_exact_100_percent_wavelength(self):
        wavelengths, intensity = monochromatic_emission(
            np.array([394.0, 396.0]), 395.0
        )
        np.testing.assert_array_equal(wavelengths, [394.0, 395.0, 396.0])
        np.testing.assert_array_equal(intensity, [0.0, 100.0, 0.0])

    def test_actinometer_config_does_not_require_an_led_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            values = _example_actinometer_config(temporary)
            validate_config(AnalysisConfig(values, EXAMPLE_DIRECTORY))

    def test_actinometer_config_requires_the_nominal_wavelength(self):
        with tempfile.TemporaryDirectory() as temporary:
            values = _example_actinometer_config(temporary)
            values["experiment"].pop("irradiation_wavelength_nm")
            with self.assertRaisesRegex(ConfigError, "irradiation_wavelength_nm is required"):
                validate_config(AnalysisConfig(values, EXAMPLE_DIRECTORY))

    def test_example_4_power_converted_to_actinometer_flux_gives_similar_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            values = _example_actinometer_config(temporary)
            output = run_analysis(AnalysisConfig(values, EXAMPLE_DIRECTORY))

        np.testing.assert_allclose(
            output.result.yield_fit.values, [0.25, 0.21], atol=0.006, rtol=0
        )
        np.testing.assert_allclose(
            output.result.yield_errors, [0.012, 0.029], atol=0.002, rtol=0
        )
        reactant_pss = output.result.extrapolated_pss[0] / output.result.extrapolated_pss.sum()
        self.assertAlmostEqual(reactant_pss, 0.233, delta=0.002)
        summary = result_summary(output.result, output.data, 395)
        self.assertEqual(summary["experiment"]["irradiation_source"],
                         "chemical_actinometer")
        self.assertIn("optimizer_and_photon_flux",
                      summary["quantum_yield_error_components_percent"])
        self.assertEqual(
            summary["experiment"]["photon_flux_formatted_mol_s"],
            {"value": "4.8", "error": "0.1", "exponent": -9},
        )
        self.assertEqual(output.files, ())

    def test_nipe_pre_plateau_uncertainty_includes_epsilon_range(self):
        with tempfile.TemporaryDirectory() as temporary:
            values = json.loads(
                (EXAMPLE_DIRECTORY / "analysis.json").read_text(encoding="utf-8")
            )
            values["fit"]["method"] = "nipe"
            values["outputs"].update({
                "directory": temporary,
                "write_text": False,
                "write_figures": False,
                "write_json": False,
                "write_config": False,
                "write_detailed_data": False,
            })
            output = run_analysis(AnalysisConfig(values, EXAMPLE_DIRECTORY))

        uncertainty = output.result.epsilon_uncertainty
        self.assertGreater(uncertainty.nipe_window_bound_combination_count, 1)
        self.assertIsNotNone(uncertainty.nipe_window_combined_errors)
        summary = result_summary(output.result, output.data, 395)
        windows = summary["nipe"]["pre_plateau_window_analysis"]
        self.assertEqual(
            windows["reported_error_source"], "window_fit_and_epsilon_range"
        )
        self.assertIn("epsilon_range", windows)
        self.assertEqual(summary["ab_model_assessment"]["status"], "warning")
        counts = summary["ab_model_assessment"]["epsilon_bound_status_counts"]
        self.assertEqual(sum(counts.values()), uncertainty.bound_combination_count)


if __name__ == "__main__":
    unittest.main()
