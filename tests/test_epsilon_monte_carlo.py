"""Joint sampling, uncertainty reporting, and legacy ε-range compatibility."""

from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import warnings
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import truncnorm

from autoqy_core.config import AnalysisConfig, ConfigError, validate_config
from autoqy_core.epsilon import (calculate_epsilon_statistics, export_epsilon_csv,
                                export_nmr_subtraction_csv, reconstruct_product_from_nmr)
from autoqy_core.epsilon_monte_carlo import (
    _endpoint_error, _summary_precision, _positive_normal, _reference_sampler, _reconstruction_sampler,
    paired_epsilon_sampler, run_monte_carlo,
)
from autoqy_core.epsilon_uncertainty import load_epsilon_envelope
from autoqy_core.output import result_summary, write_detailed_data, write_results
from autoqy_core.pipeline import run_analysis_pipeline
from autoqy_core.plotting import write_figure
from autoqy_core.tools.analysis_gui import (
    _configuration, _epsilon_reliability_check, _form_values, _interactive_figures,
    _yield_card, create_app,
)
from test_fit_methods import EXAMPLE, synthetic_input


def reference_frame(mean, error=.02):
    mean = np.asarray(mean, float)
    return pd.DataFrame({
        'Wavelength_nm': np.arange(len(mean), dtype=float) + 400,
        'Epsilon_mean_M-1_cm-1': mean,
        'Epsilon_SD_M-1_cm-1': np.abs(mean) * error,
        'Epsilon_SEM_M-1_cm-1': np.abs(mean) * error / np.sqrt(3),
        'Epsilon_M-1_cm-1__1': mean * (1 - error),
        'Epsilon_M-1_cm-1__2': mean,
        'Epsilon_M-1_cm-1__3': mean * (1 + error),
    })


class MonteCarloSamplingTests(unittest.TestCase):
    def test_mean_and_sd_precision_is_separate_from_interval_endpoint_precision(self):
        values = np.random.default_rng(4).normal([.5, .75], [.0075, .05], size=(2048, 2))
        stable, errors, tolerances = _summary_precision(values)
        self.assertTrue(stable)
        self.assertTrue(np.all(2 * errors <= tolerances))
        self.assertFalse(_summary_precision(values[:256])[0])
        self.assertTrue(np.any(2 * _endpoint_error(values, .95) > 1e-6))

    def test_sd_precision_accounts_for_non_gaussian_fourth_moment(self):
        values = np.random.default_rng(2).standard_t(3, size=(2048, 2)) * .01 + .5
        _, errors, _ = _summary_precision(values)
        sd = values.std(axis=0, ddof=1)
        expected = ((values - values.mean(axis=0)) ** 2).std(axis=0, ddof=1) / (2 * sd * np.sqrt(len(values)))
        self.assertTrue(np.all(errors[1] >= expected * (1 - 1e-14)))
        self.assertTrue(np.any(errors[1] > sd / np.sqrt(2 * (len(values) - 1)) * 2))

    def test_identical_draws_have_no_simulation_error(self):
        stable, errors, tolerances = _summary_precision(np.tile([.3, .18], (128, 1)))
        self.assertTrue(stable)
        np.testing.assert_array_equal(errors, np.zeros((2, 2)))

    def test_positive_normal_is_truncated_not_clipped_and_matches_quantiles(self):
        probabilities = (np.arange(10000) + .5) / 10000
        samples = _positive_normal(-.5, 1., probabilities)
        self.assertTrue(np.all(samples > 0))
        expected_mean, expected_var = truncnorm.stats(.5, np.inf, loc=-.5, scale=1, moments='mv')
        self.assertAlmostEqual(float(samples.mean()), float(expected_mean), places=3)
        self.assertAlmostEqual(float(samples.var()), float(expected_var), places=3)
        self.assertEqual(float(_positive_normal(0, 0, .3)), 0)
        with self.assertRaisesRegex(ValueError, 'zero uncertainty'):
            _positive_normal(-1, 0, .5)

    def test_endpoint_stability_does_not_depend_only_on_small_batches(self):
        # Identical batch quantiles can arise by chance; unresolved tails still
        # need an order-statistic check before claiming endpoint stability.
        one_batch = np.column_stack([np.linspace(.2,.8,32),np.linspace(.1,.9,32)])
        error = _endpoint_error(np.tile(one_batch,(4,1)),.95)
        self.assertTrue(np.all(error > .005))
        fine = np.column_stack([np.linspace(.2,.8,4096),np.linspace(.1,.9,4096)])
        fine = fine[np.random.default_rng(1).permutation(len(fine))]
        self.assertTrue(np.all(_endpoint_error(fine,.95) < error))

    def test_reference_draws_preserve_whole_curve_rank_dependence(self):
        frame = reference_frame([100, 200, 300])
        sample = _reference_sampler(frame, 'sd')
        rng = np.random.default_rng(3)
        draws = np.array([sample(rng) for _ in range(100)])
        np.testing.assert_allclose(draws[:, 1], draws[:, 0] * 2)
        np.testing.assert_allclose(draws[:, 2], draws[:, 0] * 3)
        # SD versus SEM retains the specified scale rather than silently switching.
        sd = np.std(draws[:, 0], ddof=1)
        rng = np.random.default_rng(3)
        sem_draws = np.array([_reference_sampler(frame, 'sem')(rng) for _ in range(100)])
        self.assertAlmostEqual(np.std(sem_draws[:, 0], ddof=1) / sd, 1/np.sqrt(3), places=5)

    def test_missing_replicates_is_actionable_instead_of_assuming_independence(self):
        frame = reference_frame([100, 200]).drop(columns=['Epsilon_M-1_cm-1__1',
                                                       'Epsilon_M-1_cm-1__2'])
        with self.assertRaisesRegex(ValueError, 'ε ranges'):
            _reference_sampler(frame, 'sd')

    def reconstruction(self):
        x = np.arange(400., 405.)
        stable = np.array([1, 2, 3, 2, 1.])
        mixture = np.array([.7, 1.4, 2.5, 2.2, 1.4])
        epsilon = calculate_epsilon_statistics(x, np.column_stack([stable*.98, stable, stable*1.02]),
                                             np.full(3, .001), np.ones(3))
        result = reconstruct_product_from_nmr(x, stable, mixture, epsilon, 40, 1)
        from io import StringIO
        reference = pd.read_csv(StringIO(export_epsilon_csv(epsilon, ['a','b','c'])))
        derived = pd.read_csv(StringIO(export_nmr_subtraction_csv(result)))
        return x, reference, derived

    def test_reconstruction_metadata_and_legacy_inference_reuse_reference_scale(self):
        x, reference, derived = self.reconstruction()
        before = derived.copy(deep=True)
        # A deterministic NMR fraction isolates the shared-reference relationship.
        derived['Composition_uncertainty_percentage_points'] = 0
        a, f = _reconstruction_sampler(derived, reference, x)
        old = derived.drop(columns='Stable_fraction_in_PSS')
        b, inferred = _reconstruction_sampler(old, reference, x)
        self.assertAlmostEqual(f, .4)
        self.assertAlmostEqual(inferred, .4)
        nominal_reference = reference['Epsilon_mean_M-1_cm-1'].to_numpy()
        first = a(np.random.default_rng(0), nominal_reference)
        second = b(np.random.default_rng(0), nominal_reference * 1.1)
        np.testing.assert_allclose(second, first * 1.1)
        pd.testing.assert_series_equal(derived['PSS_normalized_to_reactant'],
                                      before['PSS_normalized_to_reactant'])

    def test_pairing_works_in_either_kinetic_direction_and_rejects_wrong_reference(self):
        x, reference, derived = self.reconstruction()
        with TemporaryDirectory() as folder:
            rp, dp = Path(folder)/'reference.csv', Path(folder)/'derived.csv'
            reference.to_csv(rp,index=False)
            derived.to_csv(dp,index=False)
            r, d = load_epsilon_envelope(rp), load_epsilon_envelope(dp)
            forward, _ = paired_epsilon_sampler(r,d)
            reverse, _ = paired_epsilon_sampler(d,r)
            a,b = forward(np.random.default_rng(10)), reverse(np.random.default_rng(10))
            np.testing.assert_array_equal(a[0],b[1])
            np.testing.assert_array_equal(a[1],b[0])
            reference['Epsilon_mean_M-1_cm-1'] *= 2
            with self.assertRaisesRegex(ValueError, 'shared reference'):
                _reconstruction_sampler(derived, reference, x)

    def test_saved_reference_curve_rejects_a_different_shape_with_the_same_peak(self):
        x, reference, derived = self.reconstruction()
        derived['Calibration_reference_epsilon_M-1_cm-1'] = reference['Epsilon_mean_M-1_cm-1']
        reference.loc[0, 'Epsilon_mean_M-1_cm-1'] *= .5
        with self.assertRaisesRegex(ValueError, 'differs from the shared reference'):
            _reconstruction_sampler(derived,reference,x)

    def test_missing_composition_error_does_not_become_zero(self):
        x, reference, derived = self.reconstruction()
        messages = []
        with self.assertWarnsRegex(RuntimeWarning, 'assumed ±1 percentage point \\(1 SD\\)'):
            fallback, _ = _reconstruction_sampler(
                derived.drop(columns='Composition_uncertainty_percentage_points'), reference, x,
                assumed_input_warnings=messages)
        derived['Composition_uncertainty_percentage_points'] = 1
        explicit, _ = _reconstruction_sampler(derived, reference, x)
        values = reference['Epsilon_mean_M-1_cm-1'].to_numpy()
        rng_fallback, rng_explicit = np.random.default_rng(8), np.random.default_rng(8)
        a = np.array([fallback(rng_fallback, values) for _ in range(64)])
        b = np.array([explicit(rng_explicit, values) for _ in range(64)])
        np.testing.assert_array_equal(a, b)
        self.assertGreater(float(a[:, 0].std()), 0)
        self.assertEqual(len(messages), 1)
        self.assertIn('NMR composition uncertainty is missing', messages[0])

    def test_saved_composition_error_including_zero_is_kept_without_assumption(self):
        x, reference, derived = self.reconstruction()
        for value in (0., .25, 2.):
            derived['Composition_uncertainty_percentage_points'] = value
            messages = []
            with warnings.catch_warnings(record=True) as caught:
                _reconstruction_sampler(derived, reference, x, assumed_input_warnings=messages)
            self.assertEqual(messages, [])
            self.assertEqual(caught, [])
        for invalid in (-1., float('nan'), float('inf')):
            derived['Composition_uncertainty_percentage_points'] = invalid
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                _reconstruction_sampler(derived, reference, x)

    def test_incomplete_reconstruction_names_the_file_and_each_missing_detail(self):
        x, reference, derived = self.reconstruction()
        derived = derived.drop(columns=['Reactant_normalized', 'Composition_uncertainty_percentage_points'])
        with TemporaryDirectory() as folder:
            rp, dp = Path(folder)/'stable.csv', Path(folder)/'student-reconstruction.csv'
            reference.to_csv(rp,index=False)
            derived.to_csv(dp,index=False)
            with self.assertRaises(ValueError) as caught:
                paired_epsilon_sampler(load_epsilon_envelope(dp), load_epsilon_envelope(rp))
        message = str(caught.exception)
        self.assertIn("reactant ε file 'student-reconstruction.csv'", message)
        self.assertIn('missing the saved calibration', message)
        self.assertIn('Missing CSV columns: Reactant_normalized.', message)

    def test_bundled_older_reconstruction_uses_fallback_and_warns(self):
        from autoqy_core.epsilon_uncertainty import _read_epsilon_frame
        folder = EXAMPLE
        reference = load_epsilon_envelope(folder/'reactant_absorptivity.csv')
        frame = _read_epsilon_frame(folder/'product_absorptivity.csv')
        with self.assertWarnsRegex(RuntimeWarning, 'assumed ±1 percentage point'):
            sample, fraction = _reconstruction_sampler(
                frame, _read_epsilon_frame(reference.source_path), reference.wavelengths)
        self.assertTrue(0 <= fraction < 1)
        self.assertTrue(np.isfinite(sample(np.random.default_rng(8), reference.nominal)).all())


class MonteCarloPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = TemporaryDirectory()
        cls.data = synthetic_input()
        cls.envelopes = []
        for name, (_, curve) in zip(('reactant','product'), (cls.data.epsilon_r, cls.data.epsilon_p)):
            frame = reference_frame(curve)
            frame['Wavelength_nm'] = cls.data.wavelengths
            path = Path(cls.temporary.name) / f'{name}.csv'
            frame.to_csv(path,index=False)
            cls.envelopes.append(load_epsilon_envelope(path))
        cls.output = run_monte_carlo(cls.data, *cls.envelopes, samples=32, seed=11)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_real_fits_recover_known_yields_and_keep_nominal_data(self):
        result = self.output
        u = result.epsilon_uncertainty
        self.assertEqual(u.method, 'monte_carlo')
        self.assertEqual(u.attempted_samples,32)
        self.assertEqual(u.failed_samples,0)
        self.assertFalse(u.converged)
        np.testing.assert_allclose(result.yield_fit.values,[.3,.18],atol=2e-4)
        self.assertTrue(np.all(u.interval_lower <= [.3,.18]))
        self.assertTrue(np.all(u.interval_upper >= [.3,.18]))
        np.testing.assert_allclose(result.yield_errors,np.std(u.bound_yields,axis=0,ddof=1))
        self.assertIn('PSS reconstruction absorbance is fixed',u.scope)

    def test_seed_is_reproducible_and_local_optimizer_errors_are_not_added(self):
        again = run_monte_carlo(self.data,*self.envelopes,samples=32,seed=11)
        np.testing.assert_array_equal(self.output.epsilon_uncertainty.bound_yields,
                                      again.epsilon_uncertainty.bound_yields)
        np.testing.assert_array_equal(self.data.absorbance,synthetic_input().absorbance)

    def mocked_run(self, samples=256, fail_at=None, at_boundary=False):
        nominal = run_analysis_pipeline(self.data)
        trial = replace(nominal,yield_fit=replace(nominal.yield_fit,
                        values=np.array([1.,.18]) if at_boundary else nominal.yield_fit.values,
                        active_bounds=(at_boundary,False)))
        calls = []
        def pipeline(data,**kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                return nominal
            self.assertEqual(kwargs,{'include_power_variations':False})
            self.assertEqual(data.photon_flux_error_mol_s,0)
            self.assertEqual(data.power_error_mw,0)
            if fail_at and len(calls) == fail_at:
                raise ValueError('deliberate failure')
            return trial
        with patch('autoqy_core.pipeline.run_analysis_pipeline',side_effect=pipeline):
            result = run_monte_carlo(self.data,*self.envelopes,samples=samples,seed=2)
        return result,calls

    def test_early_stability_is_not_rejected_against_maximum_sample_count(self):
        result,calls = self.mocked_run(samples=512)
        u = result.epsilon_uncertainty
        self.assertEqual(u.attempted_samples,128)
        self.assertEqual(len(calls),129)
        self.assertTrue(u.converged)

    def test_failed_fits_and_boundary_counts_are_reported(self):
        result,_ = self.mocked_run(samples=32,fail_at=5,at_boundary=True)
        u = result.epsilon_uncertainty
        self.assertEqual(u.failed_samples,1)
        self.assertEqual(u.bound_combination_count,31)
        self.assertEqual(u.boundary_fractions[0],1)
        self.assertFalse(u.converged)
        check = _epsilon_reliability_check(result)
        self.assertEqual(check['level'],'stop')
        self.assertIn('incomplete',check['body'])
        self.assertIn('limits affect the uncertainty range',check['body'])

    def test_exports_show_coverage_scope_and_simulations(self):
        summary = result_summary(self.output,self.data,395)
        self.assertNotIn('quantum_yield_error_components_percent',summary)
        self.assertEqual(summary['quantum_yield_interval_percent']['coverage_probability'],.95)
        self.assertEqual(summary['monte_carlo']['successful_simulations'],32)
        from dash import html
        cards = _yield_card(html,'R → P',summary,'R_to_P')
        self.assertEqual(cards[-1].children,'Mean ± SD')
        self.assertIn('±',cards[1].children)
        np.testing.assert_array_equal(
            [summary['quantum_yield_percent'][k] for k in ('R_to_P','P_to_R')],
            self.output.epsilon_uncertainty.sample_mean * 100)
        np.testing.assert_array_equal(
            [summary['monte_carlo']['nominal_quantum_yield_percent'][k] for k in ('R_to_P','P_to_R')],
            self.output.yield_fit.values * 100)
        with TemporaryDirectory() as folder:
            stem = Path(folder)/'results'
            write_results(stem.with_suffix('.txt'),self.output,self.data,395)
            text = stem.with_suffix('.txt').read_text(encoding='utf-8')
            self.assertIn('PSS reconstruction absorbance is held fixed',text)
            self.assertIn('simulation limit was reached',text)
            self.assertNotIn('Set Maximum simulations',text)
            self.assertNotIn('provisional',text)
            paths = write_detailed_data(stem,self.output,self.data)
            self.assertEqual(len(paths),3)
            self.assertEqual(len(pd.read_csv(paths[-1])),32)
            self.assertIn('reactant_concentration_fit_mc_interval_lower_M',pd.read_csv(paths[0]).columns)
            captured = []
            with patch('matplotlib.figure.Figure.savefig',side_effect=lambda *a,**k:captured.append(plt.gcf())):
                write_figure(stem,self.output,self.data)
                texts = ' '.join(t.get_text() for ax in captured[0].axes for t in ax.texts)
                self.assertIn('±',texts)
                self.assertNotIn('(',texts)
                self.assertNotIn('interval',texts)
                self.assertEqual(len(captured[0].axes[1].texts),1)
                self.assertIsNone(captured[0]._suptitle)
                residual = captured[0].axes[2]
                self.assertIn('95% interval',residual.get_legend_handles_labels()[1])
                np.testing.assert_allclose(
                    residual.collections[0].get_paths()[0].vertices[:,1].min(),
                    self.output.epsilon_uncertainty.fraction_residual_minimum.min())
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        figures = _interactive_figures(go,make_subplots,self.output,self.data,99,395)
        self.assertEqual(figures[0].layout.title.text,'Concentration fit')
        self.assertTrue(any(t.name == '95% interval' and t.fill == 'tonexty' for t in figures[1].data))

    def test_approximate_interval_does_not_mark_stable_mean_and_sd_as_provisional(self):
        u = replace(self.output.epsilon_uncertainty, converged=True, interval_converged=False,
                    bound_jacobian_conditions=np.ones(32), boundary_fractions=np.zeros(2))
        result = replace(self.output, epsilon_uncertainty=u)
        check = _epsilon_reliability_check(result)
        self.assertEqual(check['level'], 'ok')
        self.assertIn('mean and error passed', check['body'])
        self.assertNotIn('More simulations', check['body'])
        self.assertIn('interval limits are approximate', check['details'])
        summary = result_summary(result, self.data, 395)
        self.assertTrue(summary['monte_carlo']['mean_and_error_stable'])
        self.assertFalse(summary['monte_carlo']['endpoints_stable'])
        self.assertFalse(summary['quantum_yield_interval_percent']['endpoints_stable'])

    def test_missing_nmr_error_runs_and_records_a_visible_assumption(self):
        # The fixture's calibrated curves cover the same axis as its reference.
        # Match that axis to a real synthetic kinetic input for an end-to-end run.
        x = self.data.wavelengths
        stable = self.data.epsilon_r[1]
        mixture = .4 * stable + .6 * self.data.epsilon_p[1]
        epsilon = calculate_epsilon_statistics(x, np.column_stack([stable*.98, stable, stable*1.02]),
                                             np.ones(3), np.ones(3))
        reconstructed = reconstruct_product_from_nmr(x, stable, mixture, epsilon, 40, 1)
        from io import StringIO
        derived = pd.read_csv(StringIO(export_nmr_subtraction_csv(reconstructed)))
        derived = derived.drop(columns='Composition_uncertainty_percentage_points')
        with TemporaryDirectory() as folder:
            rp, dp = Path(folder)/'reference.csv', Path(folder)/'old-reconstruction.csv'
            rp.write_text(export_epsilon_csv(epsilon, ['a','b','c']), encoding='utf-8')
            derived.to_csv(dp, index=False)
            r, p = load_epsilon_envelope(rp), load_epsilon_envelope(dp)
            nominal = run_analysis_pipeline(self.data)
            with self.assertWarnsRegex(RuntimeWarning, "product ε file 'old-reconstruction.csv'"), \
                 patch('autoqy_core.pipeline.run_analysis_pipeline', return_value=nominal):
                result = run_monte_carlo(self.data, r, p, samples=32)
            u = result.epsilon_uncertainty
            self.assertEqual(u.attempted_samples, 32)
            self.assertEqual(u.failed_samples, 0)
            self.assertEqual(len(u.assumed_input_warnings), 1)
            visible = replace(u, converged=True, boundary_fractions=np.zeros(2),
                              bound_jacobian_conditions=np.ones(32))
            check = _epsilon_reliability_check(replace(result, epsilon_uncertainty=visible))
            self.assertEqual(check['level'], 'warning')
            self.assertIn('assumed ±1 percentage point', check['body'])
            summary = result_summary(result, self.data, 395)
            self.assertEqual(summary['monte_carlo']['assumed_input_warnings'], list(u.assumed_input_warnings))
            self.assertIn(u.assumed_input_warnings[0], summary['monte_carlo']['sampling_assumptions'])
            path = Path(folder)/'result.txt'
            write_results(path, result, self.data, 395)
            self.assertIn('assumed ±1 percentage point', path.read_text(encoding='utf-8'))

    def test_sampling_stops_on_summary_precision_even_with_unresolved_endpoints(self):
        nominal = run_analysis_pipeline(self.data)
        rng = np.random.default_rng(4)
        draws = iter(rng.normal([.5, .75], [.0075, .05], size=(2048, 2)))
        def pipeline(data, **kwargs):
            return (replace(nominal, yield_fit=replace(nominal.yield_fit, values=next(draws)))
                    if kwargs else nominal)
        with patch('autoqy_core.pipeline.run_analysis_pipeline', side_effect=pipeline):
            result = run_monte_carlo(self.data, *self.envelopes, samples=2048, tolerance=1e-6)
        u = result.epsilon_uncertainty
        self.assertTrue(u.converged)
        self.assertFalse(u.interval_converged)
        self.assertLess(u.attempted_samples, 2048)


    def test_runner_wires_mc_settings_and_protects_additional_output_file(self):
        from autoqy_core.runner import run_analysis, _write_outputs
        values = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        values['uncertainty']['epsilon'].update(method='monte_carlo',samples=32,seed=14)
        for key in ('write_text','write_figures','write_json','write_config','write_detailed_data'):
            values['outputs'][key] = False
        config = AnalysisConfig(values,EXAMPLE)
        with TemporaryDirectory() as folder, \
             patch('autoqy_core.epsilon_monte_carlo.run_monte_carlo',return_value=self.output) as mc:
            run_analysis(config,folder)
            self.assertEqual(mc.call_args.kwargs['samples'],32)
            self.assertEqual(mc.call_args.kwargs['seed'],14)
            values['outputs']['overwrite'] = False
            values['outputs']['write_detailed_data'] = True
            path = Path(folder)/f"{values['outputs']['stem']}_uncertainty.csv"
            path.write_text('preserve')
            with self.assertRaises(FileExistsError):
                _write_outputs(config,self.data,self.output,folder)
            self.assertEqual(path.read_text(),'preserve')


class MonteCarloConfigurationTests(unittest.TestCase):
    def test_gui_preserves_mc_settings_and_keeps_range_option(self):
        values = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        settings = {'method':'monte_carlo','error_metric':'sem','samples':512,'seed':13,
                    'coverage_probability':.9,'endpoint_tolerance':.002}
        values['uncertainty']['epsilon'] = settings
        form = _form_values(values)
        self.assertEqual(_configuration(form)['uncertainty']['epsilon'],settings | {'workers': 1})
        for method in ('none','deterministic_extremes','monte_carlo'):
            values['uncertainty']['epsilon']['method'] = method
            validate_config(AnalysisConfig(values,EXAMPLE))
        def walk(root):
            if isinstance(root,(list,tuple)):
                for child in root: yield from walk(child)
            else:
                yield root
                if getattr(root,'children',None) is not None:
                    yield from walk(root.children)
        dropdown = next(c for c in walk(create_app().layout) if getattr(c,'id',None) ==
                        {'type':'analysis-field','name':'epsilon_method'})
        self.assertEqual({o['value'] for o in dropdown.options},
                         {'none','deterministic_extremes','monte_carlo'})

    def test_browser_number_controls_accept_defaults_and_saved_positive_decimals(self):
        from decimal import Decimal
        def walk(root):
            if isinstance(root,(list,tuple)):
                for child in root: yield from walk(child)
            else:
                yield root
                if getattr(root,'children',None) is not None:
                    yield from walk(root.children)
        controls = {component.id['name']:component for component in walk(create_app().layout)
                    if isinstance(getattr(component,'id',None),dict)
                    and component.id.get('type') == 'analysis-field'}
        def browser_accepts(control,value):
            # HTML number inputs use min as the step base when min is present.
            numeric = Decimal(str(value))
            minimum, maximum = getattr(control,'min',None),getattr(control,'max',None)
            if minimum is not None and numeric < Decimal(str(minimum)): return False
            if maximum is not None and numeric > Decimal(str(maximum)): return False
            step = getattr(control,'step',1)
            if step == 'any': return True
            base = minimum if minimum is not None else getattr(control,'value',0)
            steps = (numeric-Decimal(str(base)))/Decimal(str(step))
            return steps == steps.to_integral_value()
        document = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        form = _form_values(document)  # Older ε-range JSON has no Monte Carlo fields.
        self.assertEqual(form['mc_tolerance'],.005)
        self.assertEqual(form['mc_samples'],2048)
        self.assertEqual(controls['mc_samples'].value,2048)
        for tolerance in (.005,.002,.000001,.00035):
            for coverage in (.95,.955,.9999):
                with self.subTest(tolerance=tolerance,coverage=coverage):
                    self.assertTrue(browser_accepts(controls['mc_tolerance'],tolerance))
                    self.assertTrue(browser_accepts(controls['mc_coverage'],coverage))
                    form.update(epsilon_method='monte_carlo',mc_tolerance=tolerance,mc_coverage=coverage)
                    saved = _configuration(form)
                    validate_config(AnalysisConfig(saved,EXAMPLE))
                    self.assertEqual(_form_values(saved)['mc_tolerance'],tolerance)
                    self.assertEqual(_form_values(saved)['mc_coverage'],coverage)

    def test_empty_tolerance_has_a_useful_message_without_hiding_zero_or_negative_values(self):
        document = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        form = _form_values(document)
        form['epsilon_method'] = 'monte_carlo'
        for invalid in (None,0,-.001):
            form['mc_tolerance'] = invalid
            with self.assertRaisesRegex(ConfigError,'Endpoint tolerance greater than zero'):
                validate_config(AnalysisConfig(_configuration(form),EXAMPLE))

    def test_invalid_mc_settings_are_rejected(self):
        baseline = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        for setting in ({'samples':True},{'samples':0},{'samples':32.5},{'seed':-1},
                        {'coverage_probability':1},{'coverage_probability':float('nan')},
                        {'endpoint_tolerance':0},{'endpoint_tolerance':float('inf')}):
            values = deepcopy(baseline)
            values['uncertainty']['epsilon'].update(method='monte_carlo',**setting)
            with self.subTest(setting=setting),self.assertRaises(ConfigError):
                validate_config(AnalysisConfig(values,EXAMPLE))


if __name__ == '__main__':
    unittest.main()
