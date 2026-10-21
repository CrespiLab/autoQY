"""Spawned workers must preserve the draw sequence and resumable checkpoints."""

from dataclasses import replace
from contextlib import contextmanager
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import time
import unittest
from unittest.mock import patch

import numpy as np

from autoqy_core.config import AnalysisConfig, ConfigError, validate_config
from autoqy_core.epsilon_monte_carlo import _summary_precision, run_monte_carlo
from autoqy_core.epsilon_uncertainty import load_epsilon_envelope
from autoqy_core.monte_carlo_execution import cpu_worker_limits, process_pool, resolve_workers, trial_results
from autoqy_core.tools.analysis_gui import _configuration, _form_values, create_app
from test_epsilon_monte_carlo import reference_frame
from test_fit_methods import EXAMPLE, synthetic_input


def worker_details():
    from threadpoolctl import threadpool_info
    time.sleep(.2)
    return os.getpid(), [item['num_threads'] for item in threadpool_info()]


def components(root):
    if isinstance(root, (list, tuple)):
        for child in root:
            yield from components(child)
    else:
        yield root
        if getattr(root, 'children', None) is not None:
            yield from components(root.children)


class CpuPolicyTests(unittest.TestCase):
    def test_defaults_leave_capacity_and_cap_the_initial_suggestion(self):
        for cpus, suggested, maximum in ((1,1,1), (2,1,1), (3,2,2), (4,2,3), (8,4,7), (16,4,15), (128,4,61)):
            with self.subTest(cpus=cpus), patch('autoqy_core.monte_carlo_execution.available_cpu_count', return_value=cpus):
                self.assertEqual(cpu_worker_limits(), (cpus, suggested, maximum))
                self.assertEqual(resolve_workers(61), (maximum, cpus))
        for invalid in (True, 0, -1, 1.5, 62):
            with self.assertRaises(ValueError):
                resolve_workers(invalid)

    def test_gui_is_opt_in_and_suggests_four_of_eight_cpus(self):
        with patch('autoqy_core.monte_carlo_execution.available_cpu_count', return_value=8):
            app = create_app()
            fields = {c.id['name']: c for c in components(app.layout)
                      if isinstance(getattr(c, 'id', None), dict) and c.id.get('type') == 'analysis-field'}
            self.assertEqual(fields['mc_parallel'].value, [])
            self.assertEqual(fields['mc_workers'].value, 4)
            self.assertEqual(fields['mc_workers'].max, 7)
            self.assertTrue(fields['mc_workers'].disabled)
            document = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
            form = _form_values(document)
            form.update(epsilon_method='monte_carlo', mc_parallel=['on'], mc_workers=4)
            saved = _configuration(form)
            self.assertEqual(saved['uncertainty']['epsilon']['workers'], 4)
            restored = _form_values(saved)
            self.assertEqual(restored['mc_parallel'], ['on'])
            self.assertEqual(restored['mc_workers'], 4)
            form['mc_parallel'] = []
            self.assertEqual(_configuration(form)['uncertainty']['epsilon']['workers'], 1)

    def test_bad_json_worker_counts_have_a_clear_validation_error(self):
        for workers in (True, 0, -1, 3.5, 62):
            document = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
            document['uncertainty']['epsilon'].update(method='monte_carlo', workers=workers)
            with self.subTest(workers=workers), self.assertRaisesRegex(ConfigError, 'workers must be an integer'):
                validate_config(AnalysisConfig(document, EXAMPLE))

    def test_non_power_of_two_simulation_limit_round_trips_unchanged(self):
        document = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        document['uncertainty']['epsilon'].update(method='monte_carlo', samples=1500)
        validate_config(AnalysisConfig(document, EXAMPLE))
        saved = _configuration(_form_values(document))
        self.assertEqual(saved['uncertainty']['epsilon']['samples'], 1500)

    def test_low_cpu_machine_keeps_parallel_option_disabled(self):
        with patch('autoqy_core.monte_carlo_execution.available_cpu_count', return_value=2):
            app = create_app()
            parallel = next(c for c in components(app.layout)
                            if getattr(c, 'id', None) == {'type':'analysis-field','name':'mc_parallel'})
            self.assertTrue(parallel.options[0]['disabled'])


@unittest.skipUnless(cpu_worker_limits()[2] >= 2, 'Parallel fit tests need room for two workers and the parent')
class ParallelMonteCarloTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = TemporaryDirectory()
        cls.data = synthetic_input()
        cls.envelopes = []
        for name, (_, curve) in zip(('r','p'), (cls.data.epsilon_r, cls.data.epsilon_p)):
            frame = reference_frame(curve)
            frame['Wavelength_nm'] = cls.data.wavelengths
            path = Path(cls.folder.name)/f'{name}.csv'
            frame.to_csv(path, index=False)
            cls.envelopes.append(load_epsilon_envelope(path))
        cls.serial = run_monte_carlo(cls.data, *cls.envelopes, samples=32, seed=19)

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def assert_same_results(self, result):
        a, b = result.epsilon_uncertainty, self.serial.epsilon_uncertainty
        for name in ('bound_yields','sample_mean','sample_standard_deviation','interval_lower','interval_upper',
                     'fraction_residual_minimum','fraction_residual_maximum'):
            np.testing.assert_allclose(getattr(a,name), getattr(b,name), rtol=1e-7, atol=1e-10, err_msg=name)

    def test_real_spawned_workers_match_serial_fits_and_keep_ordered_progress(self):
        progress = []
        result = run_monte_carlo(self.data, *self.envelopes, samples=32, seed=19,
                                workers=2, progress_callback=progress.append)
        self.assert_same_results(result)
        self.assertEqual(result.epsilon_uncertainty.workers, 2)
        self.assertEqual([p['attempted'] for p in progress], sorted(p['attempted'] for p in progress))
        self.assertEqual(progress[-1]['attempted'], 32)
        self.assertTrue(any(p['phase'] == 'starting_workers' for p in progress))

    def test_worker_threads_are_limited_and_distinct_processes_are_used(self):
        with process_pool(2) as executor:
            results = [future.result(timeout=30) for future in
                       [executor.submit(worker_details) for _ in range(8)]]
        self.assertEqual(len({pid for pid, _ in results}), 2)
        for pid, counts in results:
            self.assertNotEqual(pid, os.getpid())
            self.assertTrue(counts)
            self.assertTrue(all(count == 1 for count in counts), counts)

    def test_passed_precision_check_stops_parallel_submissions_and_cached_rerun(self):
        futures, checks, progress = [], [], []

        @contextmanager
        def counted_pool(workers):
            with process_pool(workers) as executor:
                submit = executor.submit
                def counted_submit(*args, **kwargs):
                    future = submit(*args, **kwargs)
                    futures.append(future)
                    return future
                executor.submit = counted_submit
                yield executor

        def precision_after_192(values):
            # Control when the check passes to test scheduling independently
            # of the numerical criterion, which has its own sampling tests.
            checks.append(len(values))
            _, errors, tolerances = _summary_precision(values)
            return len(values) >= 192, errors, tolerances

        with TemporaryDirectory() as folder:
            with patch('autoqy_core.monte_carlo_execution.process_pool', counted_pool), \
                 patch('autoqy_core.epsilon_monte_carlo._summary_precision', side_effect=precision_after_192):
                result = run_monte_carlo(self.data, *self.envelopes, samples=1000, seed=19,
                                        workers=2, checkpoint_directory=folder,
                                        progress_callback=progress.append)
            self.assertEqual(checks, [128, 192, 192])
            self.assertEqual(len(futures), 192)
            self.assertTrue(all(future.done() for future in futures))
            self.assertEqual(progress[-1]['attempted'], 192)
            self.assertEqual(len(result.epsilon_uncertainty.bound_yields), 192)
            self.assertTrue(result.epsilon_uncertainty.converged)

            with patch('autoqy_core.monte_carlo_execution.process_pool', side_effect=AssertionError('unneeded pool')), \
                 patch('autoqy_core.epsilon_monte_carlo._summary_precision', side_effect=precision_after_192):
                resumed = run_monte_carlo(self.data, *self.envelopes, samples=1500, seed=19,
                                         workers=2, checkpoint_directory=folder)
            self.assertEqual(resumed.epsilon_uncertainty.reused_samples, 192)
            self.assertEqual(len(resumed.epsilon_uncertainty.bound_yields), 192)

    def test_switching_modes_after_interruption_keeps_completed_draws(self):
        with TemporaryDirectory() as folder:
            run_monte_carlo(self.data, *self.envelopes, samples=16, seed=19, checkpoint_directory=folder)
            def interrupt(progress):
                if progress['attempted'] == 20:
                    raise RuntimeError('simulated interruption')
            with self.assertRaisesRegex(RuntimeError, 'simulated interruption'):
                run_monte_carlo(self.data, *self.envelopes, samples=32, seed=19, workers=2,
                                checkpoint_directory=folder, progress_callback=interrupt)
            result = run_monte_carlo(self.data, *self.envelopes, samples=32, seed=19,
                                    checkpoint_directory=folder)
        self.assertEqual(result.epsilon_uncertainty.reused_samples, 20)
        self.assert_same_results(result)

    def test_fully_cached_run_does_not_start_worker_processes(self):
        with TemporaryDirectory() as folder:
            run_monte_carlo(self.data, *self.envelopes, samples=32, seed=19, checkpoint_directory=folder)
            with patch('autoqy_core.monte_carlo_execution.process_pool', side_effect=AssertionError('unneeded pool')):
                result = run_monte_carlo(self.data, *self.envelopes, samples=32, seed=19,
                                        workers=4, checkpoint_directory=folder)
        self.assertEqual(result.epsilon_uncertainty.reused_samples, 32)
        self.assert_same_results(result)

    def test_worker_fit_failure_preserves_draw_numbers_and_rng_states(self):
        outputs = []
        for workers in (1,2):
            rng = np.random.default_rng(9)
            index = 0
            def prepare():
                nonlocal index
                index += 1
                rng.random()
                return replace(self.data, fit_method='invalid' if index == 2 else self.data.fit_method)
            with process_pool(workers) as executor:
                outputs.append(list(trial_results(0, 4, workers, prepare, rng, executor)))
        for serial, parallel in zip(*outputs):
            self.assertEqual(serial[0], parallel[0])
            self.assertEqual(serial[2:], parallel[2:])
        self.assertIn('Simulation 2:', outputs[1][1][2])


if __name__ == '__main__':
    unittest.main()
