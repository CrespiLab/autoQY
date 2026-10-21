"""Independent numerical checks for the repeated-fit optimizations."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import odeint

from autoqy_core.kinetics import _baseline_correction, _rates, _solve
from autoqy_core.output import monte_carlo_messages
from autoqy_core.pipeline import run_analysis_pipeline
from test_fit_methods import synthetic_input


class KineticOptimizationTests(unittest.TestCase):
    def test_cached_rate_integral_matches_reference_for_broadband_and_single_wavelength(self):
        rng = np.random.default_rng(4)
        for size in (1, 71):
            wavelengths = np.cumsum(rng.uniform(.2, 2, size)) * 1e-9 + 400e-9
            er, ep = rng.uniform(0, 5000, (2, size))
            er[0] = ep[0] = 0  # No absorption at this wavelength.
            if size == 1:
                er[0], ep[0] = 4000, 500
            photons = rng.uniform(1e-6, 2e-6, size) / (1e-9 if size > 1 else 1)
            for yields in ((0, 0), (.3, .18), (1, 1)):
                for thermal in ((0, 0), (.002, .001)):
                    with self.subTest(size=size, yields=yields, thermal=thermal):
                        times = np.linspace(0, 400, 16)
                        initial = np.array([6e-5, 4e-5])
                        args = (wavelengths, *yields, er, ep, .8, photons, 2, *thermal)
                        expected = odeint(_rates, initial, times, args=args)
                        actual = _solve(yields, initial, times, wavelengths, er, ep,
                                        .8, photons, 2, *thermal)
                        np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=1e-12)
                        np.testing.assert_allclose(actual.sum(axis=1), initial.sum(), atol=1e-12)

    def test_cached_baseline_projection_matches_least_squares(self):
        rng = np.random.default_rng(7)
        for x in (np.linspace(250, 800, 250), np.ones(10)):
            delta = rng.normal(size=(16, len(x)))
            coordinate = (x - x.mean()) / np.ptp(x) if np.ptp(x) else np.zeros_like(x)
            for order in (-1, 0, 1):
                if order == -1:
                    expected = np.zeros_like(delta)
                else:
                    design = np.column_stack([coordinate ** degree for degree in range(order + 1)])
                    expected = (design @ np.linalg.lstsq(design, delta.T, rcond=None)[0]).T
                np.testing.assert_allclose(_baseline_correction(delta, x, order), expected, atol=1e-14)

    def test_full_fit_agrees_with_reference_ode_integration(self):
        data = synthetic_input()
        data = replace(data, thermal_rate=.001, thermal_forward_rate=.0002,
                       absorbance_baseline_order=1)
        actual = run_analysis_pipeline(data, include_power_variations=False)

        def reference(yields, initial, timestamps, wavelengths, er, ep, path, photons,
                      volume, back, forward=0):
            return odeint(_rates, initial, timestamps,
                          args=(wavelengths, *yields, er, ep, path, photons, volume, back, forward))

        with patch('autoqy_core.kinetics._solve', side_effect=reference):
            expected = run_analysis_pipeline(data, include_power_variations=False)
        np.testing.assert_allclose(actual.yield_fit.values, expected.yield_fit.values, atol=2e-4)
        np.testing.assert_allclose(actual.yield_fit.concentrations,
                                   expected.yield_fit.concentrations, atol=1e-8)


class MonteCarloMessageTests(unittest.TestCase):
    def summary(self, **changes):
        values = dict(failed_samples=0, bound_combination_count=256,
                      attempted_samples=256, requested_samples=256, converged=False,
                      boundary_fractions=np.zeros(2))
        values.update(changes)
        return SimpleNamespace(**values)

    def test_limit_reports_remaining_precision_without_automatic_doubling(self):
        text = ' '.join(monte_carlo_messages(self.summary()))
        self.assertIn('All 256 simulations finished successfully.', text)
        self.assertIn('simulation limit was reached', text)
        self.assertIn('remaining simulation error', text)
        self.assertNotIn('Set Maximum simulations', text)
        self.assertNotIn('provisional', text)
        self.assertIn('1 SD', text)
        self.assertNotIn('endpoint', text)

    def test_stable_run_does_not_ask_for_more_simulations(self):
        text = ' '.join(monte_carlo_messages(self.summary(converged=True)))
        self.assertIn('passed the precision check', text)
        self.assertNotIn('run again', text)

    def test_failed_fits_need_review_before_more_simulations(self):
        text = ' '.join(monte_carlo_messages(self.summary(failed_samples=2,
                                                        bound_combination_count=254)))
        self.assertIn('2 fits failed', text)
        self.assertIn('Technical details', text)
        self.assertNotIn('Set Maximum simulations', text)

    def test_no_repeated_doubling_at_large_limits(self):
        for limit in (4096, 8192, 100000):
            text = ' '.join(monte_carlo_messages(self.summary(requested_samples=limit)))
            self.assertIn('Technical details', text)
            self.assertNotIn('run again', text)
            self.assertNotIn('Set Maximum simulations', text)
