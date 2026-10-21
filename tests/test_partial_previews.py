import base64
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from autoqy_core.config import AnalysisConfig
from autoqy_core.io import load_timestamps
from autoqy_core.tools.analysis_gui import _configuration, _preprocessing_figure


class TimestampErrorTests(unittest.TestCase):
    def test_wrong_epsilon_file_identifies_file_header_and_correction(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'epsilon-spectra-reactant.csv'
            path.write_text('Wavelength_nm,Epsilon_mean_M-1_cm-1\n455,123\n')
            with self.assertRaises(ValueError) as caught:
                load_timestamps(path)
            message = str(caught.exception)
            for phrase in (path.name, 'line 1', 'Wavelength_nm', 'Event',
                           'Inspect this file', 'ε spectrum', 'Generic CSV'):
                self.assertIn(phrase, message)

    def test_embedded_source_filename_is_preserved(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'timestamps.csv'
            path.write_text('Wavelength_nm,value\n455,12\n')
            with self.assertRaisesRegex(ValueError, 'original-reference.csv'):
                load_timestamps(path, display_name='original-reference.csv')

    def test_valid_ahk_and_generic_timestamps_keep_their_meaning(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'log.csv'
            path.write_text('ElapsedTime (s),Event\n0,Measure\n1,LEDon\n11,LEDoff\n12,Measure\n')
            np.testing.assert_allclose(load_timestamps(path), [0, 10])
            path.write_text('irradiation_seconds\n0\n10\n')
            np.testing.assert_allclose(load_timestamps(path, 'generic_delimited'), [0, 10])
            with self.assertRaisesRegex(ValueError, 'Unexpected header'):
                load_timestamps(path, 'ahk_csv')


class PartialPreviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.files = {
            'measurement_spectra': 'Wavelength,0,1\n400,.8,.7\n450,.4,.5\n500,.1,.2\n',
            'reactant_absorptivity': 'Wavelength,value\n400,100\n450,200\n500,50\n',
            'product_absorptivity': 'Wavelength,value\n400,20\n450,100\n500,80\n',
            'led_emission': 'Wavelength,value\n' + ''.join(f'{w},{np.exp(-((w-450)/15)**2)}\n' for w in range(400,501,5)),
        }
        for name, content in self.files.items():
            (self.folder / (name + '.csv')).write_text(content)
        self.values = dict(wavelength_low=400, wavelength_high=500,
                           led_window=5, led_order=2, baseline_multiplier=10,
                           irradiation_wavelength_nm=455, volume_value='unfinished',
                           expected_pss_product_percent='unfinished')
        for name in (*self.files, 'timestamps'):
            self.values[name] = ''
            self.values['format_' + name] = 'generic_delimited'
            self.values['header_' + name] = ['on']

    def preview(self, values, embedded=None):
        document = _configuration(values, embedded, preprocessing_only=True)
        return _preprocessing_figure(go, make_subplots, AnalysisConfig(document, self.folder))

    def test_each_input_can_be_previewed_on_its_own(self):
        expected = {'measurement_spectra': 'Initial spectrum',
                    'reactant_absorptivity': 'Reactant ε',
                    'product_absorptivity': 'Product ε',
                    'led_emission': 'Processed LED (normalized)'}
        for name, trace in expected.items():
            with self.subTest(name=name):
                figure = self.preview(self.values | {name: name + '.csv'})
                self.assertIn(trace, [t.name for t in figure.data])

    def test_broken_input_does_not_hide_other_loaded_inputs(self):
        values = self.values | {name: name + '.csv' for name in self.files}
        values['product_absorptivity'] = 'missing-product.csv'
        figure = self.preview(values)
        names = [t.name for t in figure.data]
        for expected in ('Reactant ε', 'Initial spectrum', 'Processed LED (normalized)'):
            self.assertIn(expected, names)
        self.assertNotIn('Product ε', names)
        self.assertTrue(any('Product ε' in a.text and 'missing-product.csv' in a.text
                            for a in figure.layout.annotations))

    def test_incomplete_experiment_and_invalid_timestamps_do_not_block_preview(self):
        values = self.values | {'measurement_spectra': 'measurement_spectra.csv',
                                'timestamps': 'not-a-log.csv'}
        figure = self.preview(values)
        self.assertIn('Final spectrum', [t.name for t in figure.data])

    def test_bad_led_treatment_shows_original_led_and_the_problem(self):
        values = self.values | {'led_emission': 'led_emission.csv', 'led_window': 1000}
        figure = self.preview(values)
        self.assertIn('Original LED (normalized; processing unavailable)', [t.name for t in figure.data])
        self.assertTrue(any('LED processing' in a.text for a in figure.layout.annotations))

    def test_embedded_reference_can_be_previewed_without_other_inputs(self):
        payload = self.files['reactant_absorptivity'].encode()
        item = dict(filename='reference.csv', encoding='base64',
                    data=base64.b64encode(payload).decode(), sha256=hashlib.sha256(payload).hexdigest())
        embedded = {'reactant_absorptivity': dict(label='Embedded: reference.csv', input=item,
                                                format={'type': 'generic_delimited'})}
        figure = self.preview(self.values | {'reactant_absorptivity': 'Embedded: reference.csv'}, embedded)
        self.assertEqual([t.name for t in figure.data], ['Reactant ε'])

    def test_embedded_data_changes_refresh_the_callback(self):
        from autoqy_core.tools.analysis_gui import create_app
        app = create_app()
        entry = next(value for key, value in app.callback_map.items()
                     if key.startswith("spectra-figure.figure@"))
        self.assertIn({"id": "analysis-embedded-inputs", "property": "data"}, entry["inputs"])
        preview = entry["callback"].__wrapped__
        values = self.values | {"reactant_absorptivity": "Embedded: ref.csv"}
        ids = [{"type": "analysis-field", "name": key} for key in values]
        traces = []
        for peak in (100, 200):
            payload = f"Wavelength,value\n400,10\n450,{peak}\n500,20\n".encode()
            item = dict(filename="ref.csv", encoding="base64", data=base64.b64encode(payload).decode(),
                        sha256=hashlib.sha256(payload).hexdigest())
            embedded = {"reactant_absorptivity": dict(label="Embedded: ref.csv", input=item,
                                                    format={"type": "generic_delimited"})}
            figure = preview(list(values.values()), embedded, ids)
            traces.append(list(figure.data[0].y))
        self.assertEqual(traces, [[10, 100, 20], [10, 200, 20]])
