"""Recover action controls from server state after success, failure or a lost response."""

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event
import unittest
from unittest.mock import patch

from autoqy_core.tools.analysis_gui import create_app, _form_values
from test_fit_methods import EXAMPLE


class AnalysisControlTests(unittest.TestCase):
    def setUp(self):
        self.app = create_app()
        self.client = self.app.server.test_client()
        self.action_key = next(key for key in self.app.callback_map
                               if key.startswith('..analysis-action-message.children'))
        self.control_key = next(key for key in self.app.callback_map
                                if key.startswith('..run-analysis.disabled'))
        document = json.loads((EXAMPLE/'analysis.json').read_text(encoding='utf-8'))
        self.values = _form_values(document)
        self.values['config_folder'] = str(EXAMPLE)

    def outputs(self, key):
        return [{'id': output.component_id, 'property': output.component_property}
                for output in self.app.callback_map[key]['output']]

    def controls(self):
        response = self.client.post('/_dash-update-component', json={
            'output': self.control_key, 'outputs': self.outputs(self.control_key),
            'inputs': [{'id': 'analysis-progress-poll', 'property': 'n_intervals', 'value': 1}],
            'state': [], 'changedPropIds': ['analysis-progress-poll.n_intervals'],
        })
        self.assertEqual(response.status_code, 200, response.data)
        data = response.get_json()['response']
        return tuple(data[key]['disabled'] for key in
                     ('run-analysis', 'compare-fit-methods', 'save-analysis-json'))

    def payload(self, action='run-analysis'):
        return {
            'output': self.action_key, 'outputs': self.outputs(self.action_key),
            'inputs': [{'id': name, 'property': 'n_clicks', 'value': int(name == action)}
                       for name in ('save-analysis-json', 'run-analysis', 'compare-fit-methods')],
            'state': [
                {'id': '{"name":["ALL"],"type":"analysis-field"}', 'property': 'value',
                 'value': list(self.values.values())},
                {'id': '{"name":["ALL"],"type":"analysis-field"}', 'property': 'id',
                 'value': [{'type': 'analysis-field', 'name': name} for name in self.values]},
                {'id': 'analysis-embedded-inputs', 'property': 'data', 'value': {}},
            ],
            'changedPropIds': [action + '.n_clicks'],
        }

    def perform(self, action='run-analysis'):
        # A separate client emulates the independent long action request.
        return self.app.server.test_client().post('/_dash-update-component', json=self.payload(action))

    def test_idle_poll_restores_controls_without_the_original_action_response(self):
        self.assertEqual(self.controls(), (False, False, False))
        lock = self.app.server.config['AUTOQY_ANALYSIS_LOCK']
        with lock:
            self.assertEqual(self.controls(), (True, True, True))
        # No action response is delivered to the browser. Its next independent
        # poll must still reset a previously disabled button.
        self.assertEqual(self.controls(), (False, False, False))

    def test_failed_request_releases_controls_and_duplicate_click_cannot_start_another_fit(self):
        started, release = Event(), Event()
        def fit(*args, **kwargs):
            started.set()
            if not release.wait(10):
                raise TimeoutError('test did not release fit')
            raise RuntimeError('deliberate fit failure')
        with patch('autoqy_core.tools.analysis_gui._input_checks', return_value=[]), \
             patch('autoqy_core.runner.run_analysis', side_effect=fit) as run, \
             ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(self.perform)
            try:
                self.assertTrue(started.wait(5), 'analysis did not reach fit')
                self.assertEqual(self.controls(), (True, True, True))
                duplicate = self.perform()
                self.assertEqual(duplicate.status_code, 200, duplicate.data)
                self.assertEqual(duplicate.get_json()['response'], {})
                self.assertEqual(run.call_count, 1)
                self.assertEqual(self.controls(), (True, True, True))
            finally:
                release.set()
            response = future.result(timeout=10)
        self.assertEqual(response.status_code, 200, response.data)
        self.assertIn('deliberate fit failure', response.get_data(as_text=True))
        self.assertEqual(self.controls(), (False, False, False))

    def test_save_dialog_cancel_releases_all_action_controls(self):
        started, release = Event(), Event()
        def choose(*args):
            started.set()
            if not release.wait(10):
                raise TimeoutError('test did not close save dialog')
            return None
        with patch('autoqy_core.tools.analysis_gui._input_checks', return_value=[]), \
             patch('autoqy_core.tools.analysis_gui._choose_save_json', side_effect=choose), \
             ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(self.perform, 'save-analysis-json')
            try:
                self.assertTrue(started.wait(5), 'save did not reach dialog')
                self.assertEqual(self.controls(), (True, True, True))
            finally:
                release.set()
            response = future.result(timeout=10)
        self.assertEqual(response.status_code, 200, response.data)
        self.assertIn('Save cancelled', response.get_data(as_text=True))
        self.assertEqual(self.controls(), (False, False, False))

    def test_configuration_error_also_releases_controls(self):
        self.values['mc_tolerance'] = 0
        self.values['epsilon_method'] = 'monte_carlo'
        response = self.perform()
        self.assertEqual(response.status_code, 200, response.data)
        self.assertIn('Endpoint tolerance', response.get_data(as_text=True))
        self.assertEqual(self.controls(), (False, False, False))

    def test_successful_save_releases_controls_without_starting_analysis(self):
        with TemporaryDirectory() as folder, \
             patch('autoqy_core.tools.analysis_gui._input_checks', return_value=[]), \
             patch('autoqy_core.tools.analysis_gui._choose_save_json',
                   return_value=str(Path(folder)/'analysis.json')), \
             patch('autoqy_core.runner.run_analysis') as run:
            response = self.perform('save-analysis-json')
            self.assertEqual(response.status_code, 200, response.data)
            self.assertTrue((Path(folder)/'analysis.json').is_file())
            self.assertIn('Saved and validated', response.get_data(as_text=True))
            run.assert_not_called()
        self.assertEqual(self.controls(), (False, False, False))

    def test_unhandled_callback_error_cannot_leave_controls_disabled(self):
        # Fail before the action's usual error renderer is entered.
        with patch('autoqy_core.tools.analysis_gui._empty_figure', side_effect=RuntimeError('lost action response')), \
             self.assertLogs(self.app.server.logger, level='ERROR'):
            response = self.perform()
        self.assertEqual(response.status_code, 500)
        self.assertEqual(self.controls(), (False, False, False))

    def test_clicks_disable_controls_immediately_and_independent_poll_can_restore_them(self):
        dependencies = self.client.get('/_dash-dependencies').get_json()
        action = next(item for item in dependencies if item['output'] == self.action_key)
        self.assertTrue(all(action['running']['running'][name + '.disabled'] for name in
                            ('run-analysis', 'compare-fit-methods', 'save-analysis-json')))
        self.assertEqual(self.controls(), (False, False, False))


if __name__ == '__main__':
    unittest.main()
