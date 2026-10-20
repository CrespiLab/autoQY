"""Supported fit methods and full-spectrum ODE initialization."""

from dataclasses import replace
import json
from pathlib import Path
import unittest

import numpy as np

from autoqy_core.config import AnalysisConfig, ConfigError, validate_config
from autoqy_core.kinetics import _solve
from autoqy_core.pipeline import AnalysisInput, run_analysis_pipeline
from autoqy_core.spectra import fit_concentrations


EXAMPLE = (Path(__file__).parents[1] / "ExampleData"
           / "Example-4_395nm-EpsilonError/generic_inputs")


def synthetic_input():
    wavelengths = np.arange(360.0, 502.0, 2.0)
    epsilon_r = 5000 * np.exp(-((wavelengths - 390) / 22) ** 2)
    epsilon_p = 3500 * np.exp(-((wavelengths - 440) / 25) ** 2)
    timestamps = np.linspace(0, 500, 16)
    wavelength = 395.0
    photon_flux = 1e-9
    concentrations = _solve(
        (0.3, 0.18), np.array([1e-4, 0.0]), timestamps,
        np.array([wavelength * 1e-9]),
        np.array([np.interp(wavelength, wavelengths, epsilon_r)]),
        np.array([np.interp(wavelength, wavelengths, epsilon_p)]),
        1, np.array([photon_flux * 1000]), 1, 0,
    )
    absorbance = (concentrations @ np.vstack((epsilon_r, epsilon_p))).T
    return AnalysisInput(
        wavelengths=wavelengths, absorbance=absorbance, timestamps=timestamps,
        epsilon_r=(wavelengths, epsilon_r), epsilon_p=(wavelengths, epsilon_p),
        led=(wavelengths, np.zeros_like(wavelengths)),
        power_mw=None, power_error_mw=0, volume_ml=1,
        wavelength_limits=(360, 500), fit_method="ode_absorbance",
        absorbance_baseline_order=-1, photon_flux_mol_s=photon_flux,
        irradiation_wavelength_nm=wavelength,
    )


class FitMethodTests(unittest.TestCase):
    def test_configuration_accepts_only_supported_methods(self):
        values = json.loads((EXAMPLE / "analysis.json").read_text(encoding="utf-8"))
        for method in ("concentrations", "emission", "ode_absorbance"):
            with self.subTest(method=method):
                values["fit"]["method"] = method
                validate_config(AnalysisConfig(values, EXAMPLE))
        values["fit"]["method"] = "unsupported"
        with self.assertRaisesRegex(
                ConfigError, "fit.method must be concentrations, emission, or ode_absorbance"):
            validate_config(AnalysisConfig(values, EXAMPLE))

    def test_ode_uses_independent_spectral_recovery_and_recovers_yields(self):
        data = synthetic_input()
        result = run_analysis_pipeline(data)
        independent = fit_concentrations(
            result.absorbance, result.wavelengths, result.epsilon_r,
            result.epsilon_p, data.path_length_cm,
        )
        np.testing.assert_allclose(
            result.concentration_fit.concentrations, independent.concentrations
        )
        self.assertTrue(result.yield_fit.optimizer_success)
        np.testing.assert_allclose(result.yield_fit.values, [0.3, 0.18], atol=2e-4)
        np.testing.assert_allclose(
            result.yield_fit.concentrations.sum(axis=1), 1e-4, rtol=1e-5
        )

    def test_pipeline_rejects_unsupported_method(self):
        with self.assertRaisesRegex(ValueError, "Unsupported fit method: unsupported"):
            run_analysis_pipeline(replace(synthetic_input(), fit_method="unsupported"))


if __name__ == "__main__":
    unittest.main()
