import unittest

import numpy as np

from autoqy_core.kinetics import (
    _irradiation_inputs,
    _nipe_plateau_index,
    _solve,
    analyze_nipe_time_windows,
    fit_quantum_yields_nipe,
)
from autoqy_core.spectra import (
    assess_ab_model,
    fit_concentrations_variable_total,
)


class NIPETests(unittest.TestCase):
    def setUp(self):
        self.wavelengths = np.array([350.0, 365.0, 400.0])
        self.epsilon_r = np.array([12_000.0, 10_000.0, 8_000.0])
        self.epsilon_p = np.array([6_000.0, 8_000.0, 10_000.0])
        self.emission = np.array([0.0, 100.0, 0.0])
        self.timestamps = np.linspace(0, 40, 41)
        self.photon_flux = 1e-8
        irradiation = _irradiation_inputs(
            self.wavelengths, self.emission, self.epsilon_r, self.epsilon_p,
            None, self.photon_flux, 365.0,
        )
        wavelengths_m, irradiated_epsilon_r, irradiated_epsilon_p, photons = irradiation
        self.expected_yields = np.array([0.22, 0.08])
        self.concentrations = _solve(
            self.expected_yields, np.array([1e-4, 0.0]), self.timestamps,
            wavelengths_m, irradiated_epsilon_r, irradiated_epsilon_p,
            1.0, photons, 2.0, 0.0,
        )
        self.absorbance = (
            self.concentrations
            @ np.vstack((self.epsilon_r, self.epsilon_p))
        ).T

    def _fit(self, concentrations, absorbance):
        return fit_quantum_yields_nipe(
            self.wavelengths, self.emission, absorbance, concentrations,
            self.timestamps, self.epsilon_r, self.epsilon_p, None, 2.0, 0.0,
            photon_flux_mol_s=self.photon_flux,
            irradiation_wavelength_nm=365.0,
        )

    def test_nipe_recovers_reversible_yields_without_boundary_assumptions(self):
        fit = self._fit(self.concentrations, self.absorbance)
        np.testing.assert_allclose(fit.values, self.expected_yields, atol=3e-3)
        self.assertGreater(fit.nipe.interval_count, len(self.timestamps))
        self.assertLess(fit.nipe.normalized_residual_rmse, 2e-3)

    def test_automatic_windows_estimate_zero_exposure_before_plateau(self):
        full = self._fit(self.concentrations, self.absorbance)
        analysis = analyze_nipe_time_windows(
            self.wavelengths, self.emission, self.absorbance,
            self.concentrations, self.timestamps, self.epsilon_r,
            self.epsilon_p, None, 2.0, 0.0,
            photon_flux_mol_s=self.photon_flux,
            irradiation_wavelength_nm=365.0,
            full_trace_values=full.values,
        )
        self.assertIsNotNone(analysis)
        self.assertGreaterEqual(len(analysis.window_start_s), 3)
        self.assertGreaterEqual(analysis.window_point_count, 4)
        self.assertLessEqual(analysis.window_point_count, 6)
        np.testing.assert_allclose(
            analysis.extrapolated_values, self.expected_yields, atol=5e-3
        )

    def test_plateau_requires_sustained_small_composition_changes(self):
        product_fraction = np.array([
            0.00, 0.25, 0.40, 0.49, 0.54, 0.56, 0.568, 0.573, 0.577,
        ])
        index = _nipe_plateau_index(product_fraction)
        self.assertEqual(index, 6)
        self.assertLess(index, len(product_fraction))

    def test_clean_two_species_series_passes_balance_check(self):
        recovered = fit_concentrations_variable_total(
            self.absorbance, self.epsilon_r, self.epsilon_p
        )
        diagnostic = assess_ab_model(recovered, self.absorbance)
        self.assertEqual(diagnostic.level, "ok")
        self.assertTrue(diagnostic.model_supported)

    def test_degrading_series_is_red_but_retains_apparent_yield_estimate(self):
        remaining = np.linspace(1.0, 0.82, len(self.timestamps))
        concentrations = self.concentrations * remaining[:, None]
        absorbance = (
            concentrations @ np.vstack((self.epsilon_r, self.epsilon_p))
        ).T
        recovered = fit_concentrations_variable_total(
            absorbance, self.epsilon_r, self.epsilon_p
        )
        diagnostic = assess_ab_model(recovered, absorbance)
        fit = self._fit(recovered.concentrations, absorbance)

        self.assertEqual(diagnostic.level, "stop")
        self.assertFalse(diagnostic.model_supported)
        self.assertGreater(diagnostic.balance_span_fraction, 0.15)
        self.assertTrue(np.isfinite(fit.values).all())
        self.assertTrue(((fit.values >= 0) & (fit.values <= 1)).all())


if __name__ == "__main__":
    unittest.main()
