"""Saved Monte Carlo runs must continue the original random sequence exactly."""

from dataclasses import replace
from contextlib import closing
import json
import hashlib
import shutil
import sqlite3
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from autoqy_core.epsilon_monte_carlo import run_monte_carlo
from autoqy_core.epsilon_uncertainty import load_epsilon_envelope
from autoqy_core.monte_carlo_checkpoint import (MonteCarloCheckpoint, checkpoint_key,
                                               _calculation_signature, legacy_checkpoint_keys)
from autoqy_core.pipeline import run_analysis_pipeline
from autoqy_core.tools.analysis_gui import create_app
from test_epsilon_monte_carlo import reference_frame
from test_fit_methods import synthetic_input


class MonteCarloCheckpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = TemporaryDirectory()
        cls.data = synthetic_input()
        cls.envelopes = []
        for name, (_, curve) in zip(('r', 'p'), (cls.data.epsilon_r, cls.data.epsilon_p)):
            frame = reference_frame(curve)
            frame['Wavelength_nm'] = cls.data.wavelengths
            path = Path(cls.folder.name)/f'{name}.csv'
            frame.to_csv(path,index=False)
            cls.envelopes.append(load_epsilon_envelope(path))
        cls.fresh = run_monte_carlo(cls.data,*cls.envelopes,samples=32,seed=3)

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def run_saved(self, folder, samples=32, **kwargs):
        return run_monte_carlo(self.data,*self.envelopes,samples=samples,seed=3,
                               checkpoint_directory=folder,**kwargs)

    def assert_same_draws(self, actual):
        a, b = actual.epsilon_uncertainty, self.fresh.epsilon_uncertainty
        for name in ('bound_yields','bound_optimizer_errors','bound_active_bounds',
                     'bound_fraction_rmse','fraction_residual_minimum','fraction_residual_maximum',
                     'sample_mean','sample_standard_deviation','interval_lower','interval_upper'):
            np.testing.assert_array_equal(getattr(a,name),getattr(b,name),err_msg=name)

    def test_increasing_limit_keeps_previous_draws_and_exact_rng_sequence(self):
        with TemporaryDirectory() as folder:
            first = self.run_saved(folder,samples=16)
            progress = []
            second = self.run_saved(folder,progress_callback=progress.append)
            self.assert_same_draws(second)
            self.assertEqual(second.epsilon_uncertainty.reused_samples,16)
            np.testing.assert_array_equal(first.epsilon_uncertainty.bound_yields,
                                          second.epsilon_uncertainty.bound_yields[:16])
            self.assertEqual(progress[0]['attempted'],16)
            self.assertEqual(progress[-1]['attempted'],32)
            self.assertEqual(progress[-1]['phase'],'complete')
            self.assertTrue(all(x['attempted'] >=16 for x in progress))
            # A smaller limit reads only the corresponding prefix, without losing later draws.
            small = self.run_saved(folder,samples=16)
            self.assertEqual(small.epsilon_uncertainty.reused_samples,16)
            np.testing.assert_array_equal(small.epsilon_uncertainty.bound_yields,
                                          first.epsilon_uncertainty.bound_yields)
            with patch('autoqy_core.pipeline.run_analysis_pipeline', wraps=run_analysis_pipeline) as fit:
                again = self.run_saved(folder)
            self.assertEqual(fit.call_count,1)  # Only the nominal fit; all simulations are reused.
            self.assertEqual(again.epsilon_uncertainty.reused_samples,32)
            self.assert_same_draws(again)

    def test_interruption_keeps_every_completed_simulation(self):
        with TemporaryDirectory() as folder:
            def interrupt(progress):
                if progress['attempted'] == 10:
                    raise RuntimeError('simulated window closure')
            with self.assertRaisesRegex(RuntimeError,'window closure'):
                self.run_saved(folder,progress_callback=interrupt)
            resumed = self.run_saved(folder)
            self.assertEqual(resumed.epsilon_uncertainty.reused_samples,10)
            self.assert_same_draws(resumed)

    def test_reopening_in_a_new_python_process_resumes_saved_draws(self):
        with TemporaryDirectory() as folder:
            self.run_saved(folder,samples=16)
            root = str(Path(__file__).parents[1])
            output = Path(folder)/'continued.npz'
            code = f'''
import sys
sys.path[:0]=[{root!r}, {str(Path(__file__).parent)!r}]
import numpy as np
from test_fit_methods import synthetic_input
from autoqy_core.epsilon_uncertainty import load_epsilon_envelope
from autoqy_core.epsilon_monte_carlo import run_monte_carlo
r=run_monte_carlo(synthetic_input(),
    load_epsilon_envelope({self.envelopes[0].source_path!r}),
    load_epsilon_envelope({self.envelopes[1].source_path!r}),
    samples=32,seed=3,checkpoint_directory={folder!r})
assert r.epsilon_uncertainty.reused_samples == 16
np.savez({str(output)!r},yields=r.epsilon_uncertainty.bound_yields)
'''
            completed = subprocess.run([sys.executable,'-c',code],capture_output=True,text=True,timeout=90)
            self.assertEqual(completed.returncode,0,completed.stderr)
            with np.load(output,allow_pickle=False) as arrays:
                np.testing.assert_array_equal(arrays['yields'],self.fresh.epsilon_uncertainty.bound_yields)

    def test_fingerprint_changes_with_inputs_fit_settings_and_seed_but_not_paths(self):
        key = checkpoint_key(self.data,*self.envelopes,'sd',3)
        for data in (replace(self.data,thermal_rate=.001),
                     replace(self.data,absorbance=self.data.absorbance+1e-5),
                     replace(self.data,initial_yields=(.1,.2)),
                     replace(self.data,power_error_mw=.02),
                     replace(self.data,fit_method='concentrations')):
            self.assertNotEqual(checkpoint_key(data,*self.envelopes,'sd',3),key)
        self.assertNotEqual(checkpoint_key(self.data,*self.envelopes,'sem',3),key)
        self.assertNotEqual(checkpoint_key(self.data,*self.envelopes,'sd',4),key)
        with TemporaryDirectory() as folder:
            path = Path(folder)/'renamed-reference.csv'
            path.write_bytes(Path(self.envelopes[0].source_path).read_bytes())
            moved = replace(self.envelopes[0],source_path=str(path))
            self.assertEqual(checkpoint_key(self.data,moved,self.envelopes[1],'sd',3),key)
            path.write_bytes(path.read_bytes()+b'\n')
            self.assertNotEqual(checkpoint_key(self.data,moved,self.envelopes[1],'sd',3),key)

    def test_simultaneous_identical_runs_are_blocked_and_lock_is_released(self):
        with TemporaryDirectory() as folder:
            key = checkpoint_key(self.data,*self.envelopes,'sd',3)
            with MonteCarloCheckpoint(folder,key):
                with self.assertRaisesRegex(ValueError,'already running'):
                    with MonteCarloCheckpoint(folder,key):
                        self.fail('Both runs acquired the same checkpoint')
            with MonteCarloCheckpoint(folder,key) as checkpoint:
                self.assertEqual(checkpoint.load(32)[:3],([],[],0))

    def test_failed_draws_are_saved_and_do_not_shift_the_random_sequence(self):
        nominal = run_analysis_pipeline(self.data)
        def pipeline(data,**kwargs):
            if kwargs.get('include_power_variations') is False:
                if data.epsilon_r[1].sum() > self.data.epsilon_r[1].sum() * 1.02:
                    raise ValueError('simulated fit failure')
            return nominal
        with TemporaryDirectory() as folder, patch('autoqy_core.pipeline.run_analysis_pipeline',side_effect=pipeline):
            fresh = run_monte_carlo(self.data,*self.envelopes,samples=64,seed=8)
            run_monte_carlo(self.data,*self.envelopes,samples=32,seed=8,checkpoint_directory=folder)
            resumed = run_monte_carlo(self.data,*self.envelopes,samples=64,seed=8,checkpoint_directory=folder)
        self.assertEqual(fresh.epsilon_uncertainty.failure_messages,resumed.epsilon_uncertainty.failure_messages)
        self.assertGreater(resumed.epsilon_uncertainty.failed_samples,0)
        np.testing.assert_array_equal(fresh.epsilon_uncertainty.bound_yields,resumed.epsilon_uncertainty.bound_yields)

    def test_new_coverage_and_tolerance_reuse_saved_draws(self):
        with TemporaryDirectory() as folder:
            self.run_saved(folder)
            changed = self.run_saved(folder,coverage=.9,tolerance=.001)
        self.assertEqual(changed.epsilon_uncertainty.reused_samples,32)
        np.testing.assert_array_equal(changed.epsilon_uncertainty.bound_yields,
                                      self.fresh.epsilon_uncertainty.bound_yields)
        self.assertEqual(changed.epsilon_uncertainty.coverage_probability,.9)

    def test_verified_legacy_draws_migrate_exactly_without_refitting(self):
        with TemporaryDirectory() as original, TemporaryDirectory() as migrated:
            output = self.run_saved(original)
            source = Path(output.epsilon_uncertainty.checkpoint_path)
            with closing(sqlite3.connect(source)) as connection:
                rows = connection.execute('SELECT * FROM draws ORDER BY iteration').fetchall()
            legacy = legacy_checkpoint_keys(self.data, *self.envelopes, 'sd', 3)
            self.assertGreaterEqual(len(legacy), 2)  # Previous v2 plus the original v1.
            with MonteCarloCheckpoint(migrated, legacy[0]) as old:
                with old.connection:
                    old.connection.executemany('INSERT INTO draws VALUES (?, ?, ?, ?)', rows)
                old_path = old.path
            original_hash = hashlib.sha256(old_path.read_bytes()).hexdigest()
            with patch('autoqy_core.pipeline.run_analysis_pipeline', wraps=run_analysis_pipeline) as fit:
                output = self.run_saved(migrated)
            self.assertEqual(fit.call_count, 1)
            self.assertEqual(output.epsilon_uncertainty.reused_samples, 32)
            self.assert_same_draws(output)
            self.assertEqual(hashlib.sha256(old_path.read_bytes()).hexdigest(), original_hash)
            with closing(sqlite3.connect(output.epsilon_uncertainty.checkpoint_path)) as connection:
                self.assertEqual(connection.execute('SELECT * FROM draws ORDER BY iteration').fetchall(), rows)

    def test_fallback_metadata_cannot_reuse_pre_fallback_draws(self):
        import pandas as pd
        with TemporaryDirectory() as folder:
            path = Path(folder)/'derived.csv'
            # Header eligibility is enough here; normal sampling validates the
            # complete reconstruction before it opens any checkpoint.
            pd.DataFrame({'Product_epsilon_M-1_cm-1': [1]}).to_csv(path, index=False)
            derived = replace(self.envelopes[0], source_path=str(path))
            compatible = legacy_checkpoint_keys(self.data, derived, self.envelopes[1], 'sd', 3)
            self.assertEqual(len(compatible), 1)  # The previous fallback-capable sampler remains compatible.
            pd.DataFrame({'Product_epsilon_M-1_cm-1': [1],
                          'Composition_uncertainty_percentage_points': [1]}).to_csv(path, index=False)
            self.assertGreaterEqual(len(legacy_checkpoint_keys(self.data, derived, self.envelopes[1], 'sd', 3)), 3)

    def test_calculation_signature_ignores_stopping_but_detects_sampling_changes(self):
        package = Path(__file__).parents[1]/'autoqy_core'
        with TemporaryDirectory() as folder:
            for name in ('epsilon_monte_carlo.py', 'epsilon_uncertainty.py', 'kinetics.py', 'pipeline.py', 'spectra.py', 'monte_carlo_draw.py'):
                shutil.copy2(package/name, Path(folder)/name)
            path = Path(folder)/'epsilon_monte_carlo.py'
            text = path.read_text(encoding='utf-8')
            original = _calculation_signature(folder)
            path.write_text(text.replace('len(records) >= 128', 'len(records) >= 256'), encoding='utf-8')
            self.assertEqual(_calculation_signature(folder), original)
            draw_path = Path(folder)/'monte_carlo_draw.py'
            draw_path.write_text(draw_path.read_text(encoding='utf-8').replace('rng.random()', 'rng.random() / 2'), encoding='utf-8')
            self.assertNotEqual(_calculation_signature(folder), original)
        with patch('autoqy_core.monte_carlo_checkpoint._calculation_signature', return_value='changed'):
            self.assertEqual(legacy_checkpoint_keys(self.data, *self.envelopes, 'sd', 3), ())


class MonteCarloProgressGuiTests(unittest.TestCase):
    def test_counter_reads_progress_while_worker_is_active(self):
        app = create_app()
        state, lock = app.server.config['AUTOQY_MC_PROGRESS']
        key = 'analysis-mc-progress.children'
        callback = app.callback_map[key]['callback'].__wrapped__
        self.assertEqual(callback(0),'')
        with lock:
            state.update(kind='monte_carlo',phase='running',active=True,attempted=384,
                         successful=383,failed=1,maximum=4096,reused=256,
                         fit_method='ode_absorbance',checkpoint_path='saved.sqlite3')
        # Exercise the polling HTTP endpoint independently of the long analysis request.
        response = app.server.test_client().post('/_dash-update-component',json={
            'output':key,'outputs':{'id':'analysis-mc-progress','property':'children'},
            'inputs':[{'id':'analysis-progress-poll','property':'n_intervals','value':1}],
            'state':[],'changedPropIds':['analysis-progress-poll.n_intervals']})
        self.assertEqual(response.status_code,200,response.data)
        text = json.dumps(response.get_json())
        self.assertIn('384 of 4096',text)
        self.assertIn('256 reused',text)
        self.assertIn('1 failed',text)
        self.assertIn('Progress',text)
        self.assertIn('Monte Carlo running',text)
        with lock:
            state.update(phase='complete',active=False)
        self.assertIn('Monte Carlo finished',json.dumps(callback(2),default=lambda c:c.to_plotly_json()))
