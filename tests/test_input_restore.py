import base64
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import struct
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from dash import no_update
import numpy as np
import pandas as pd

from autoqy_core.config import AnalysisConfig, self_contained_document
from autoqy_core.epsilon_uncertainty import load_epsilon_envelope
from autoqy_core.input_restore import restore_embedded_inputs
from autoqy_core.io import load_spectra, load_timestamps
from autoqy_core.runner import run_analysis
from autoqy_core.tools.analysis_gui import _embedded_form_inputs, create_app


EXAMPLE = (Path(__file__).parents[1] / 'ExampleData'
           / 'Example-4_395nm-EpsilonError/generic_inputs')


def embedded(payload, filename='source.dat'):
    return {'filename': filename, 'encoding': 'base64',
            'data': base64.b64encode(payload).decode('ascii'),
            'sha256': hashlib.sha256(payload).hexdigest()}


class InputRestoreTests(unittest.TestCase):
    def config(self, name, payload, spec):
        return AnalysisConfig({'inputs': {name: embedded(payload),
                               'formats': {name: spec}}}, Path.cwd())

    def example(self):
        values = json.loads((EXAMPLE / 'analysis.json').read_text(encoding='utf-8'))
        values['uncertainty']['epsilon']['method'] = 'none'
        for name in ('write_text', 'write_figures', 'write_json', 'write_config', 'write_detailed_data'):
            values['outputs'][name] = False
        return AnalysisConfig(self_contained_document(AnalysisConfig(values, EXAMPLE)), Path.cwd())

    def test_restored_example_retains_all_columns_errors_and_analysis(self):
        config = self.example()
        before = deepcopy(config.values)
        with TemporaryDirectory() as folder:
            paths = restore_embedded_inputs(config, folder)
            self.assertEqual(len(paths), 5)
            for name, path in paths.items():
                pd.testing.assert_frame_equal(
                    pd.read_csv(config.input_path(name), float_precision='round_trip'),
                    pd.read_csv(path, float_precision='round_trip'))
            for name in ('reactant_absorptivity', 'product_absorptivity'):
                original = load_epsilon_envelope(config.input_path(name))
                restored = load_epsilon_envelope(paths[name])
                for attr in ('wavelengths', 'nominal', 'lower', 'upper'):
                    np.testing.assert_array_equal(getattr(original, attr), getattr(restored, attr))
            linked = deepcopy(config.values)
            linked['inputs'].update({name: str(path) for name, path in paths.items()})
            original = run_analysis(config, Path(folder) / 'original')
            restored = run_analysis(AnalysisConfig(linked, Path(folder)), Path(folder) / 'restored')
            np.testing.assert_array_equal(original.data.absorbance, restored.data.absorbance)
            np.testing.assert_allclose(original.result.yield_fit.values, restored.result.yield_fit.values)
        self.assertEqual(config.values, before)

    def test_tab_tables_keep_extra_columns_and_ahk_events(self):
        config = self.config('measurement_spectra',
            b'Wavelength\tSignal\tWavenumbers [1/cm]\n400\t0.12345678901234567\t25000\n',
            'spectragryph_tsv')
        config.values['inputs']['timestamps'] = embedded(
            b'Event,ElapsedTime (s),Note\nMeasure,0,NA\nLEDon,2,light\nLEDoff,5,off\nMeasure,6,last\n')
        config.values['inputs']['formats']['timestamps'] = 'ahk_csv'
        with TemporaryDirectory() as folder:
            paths = restore_embedded_inputs(config, folder)
            frame = pd.read_csv(paths['measurement_spectra'], float_precision='round_trip')
            self.assertIn('Wavenumbers [1/cm]', frame)
            self.assertEqual(frame['Signal'][0], float('0.12345678901234567'))
            frame = pd.read_csv(paths['timestamps'], keep_default_na=False)
            self.assertEqual(frame['Note'][0], 'NA')
            self.assertEqual(list(frame['Event']), ['Measure', 'LEDon', 'LEDoff', 'Measure'])
            np.testing.assert_array_equal(load_timestamps(paths['timestamps']), [0, 3])

    def test_headerless_decimal_comma_table_is_converted(self):
        config = self.config('led_emission', b'ignored preamble\n400;0,2\n401;0,3\n',
            {'type': 'generic_delimited', 'delimiter': ';', 'decimal': ',',
             'header': False, 'skip_rows': 1})
        with TemporaryDirectory() as folder:
            paths = restore_embedded_inputs(config, folder)
            frame = pd.read_csv(paths['led_emission'])
            self.assertEqual(list(frame), ['Column_1', 'Column_2'])
            np.testing.assert_array_equal(frame.to_numpy(), [[400, .2], [401, .3]])

    def test_real_binary_readers_convert_specord_and_cary_to_csv(self):
        wavelengths = np.array([400, 401, 402, 403], dtype='<f4')
        signals = np.array([[.1, .2, .3, .4], [.5, .6, .7, .8]], dtype='<f4')
        specord = (b'[GENERAL]\r\nXUNITS=nm\r\nNPOINTS=4\r\nNCYCL=2\r\n'
                   b'[MESS]\r\nORIGIN=SPECORD 200 PLUS\r\n[DATA]\r\nXDATA=\r\n'
                   + wavelengths.tobytes() + b'\r\nYDATA=\r\n' + signals.tobytes())
        cary = (bytes([17]) + b'Varian UV-VIS-NIR' + bytes(110)
                + b''.join(struct.pack('<ff', w, s) for w, s in
                           zip(np.arange(430., 399., -1.), np.linspace(.1, .4, 31)))
                + bytes(32))
        for kind, payload in [('specord', specord), ('agilent_cary', cary)]:
            with self.subTest(kind=kind), TemporaryDirectory() as folder:
                config = self.config('measurement_spectra', payload, kind)
                original = load_spectra(config.input_path('measurement_spectra'), kind)
                paths = restore_embedded_inputs(config, folder)
                restored = load_spectra(paths['measurement_spectra'], 'generic_delimited')
                for a, b in zip(original, restored):
                    np.testing.assert_array_equal(a, b)

    def test_partial_bundle_does_not_need_linked_files_or_experiment_settings(self):
        config = self.config('timestamps', b'time\n0\n3\n', 'generic_delimited')
        config.values['inputs']['measurement_spectra'] = 'does-not-exist.csv'
        with TemporaryDirectory() as folder:
            paths = restore_embedded_inputs(config, folder)
            self.assertEqual(list(paths), ['timestamps'])
            np.testing.assert_array_equal(load_timestamps(paths['timestamps'], 'generic_delimited'), [0, 3])

    def test_existing_files_and_untrusted_original_names_are_safe(self):
        config = self.config('timestamps', b'time\n0\n3\n', 'generic_delimited')
        config.values['inputs']['timestamps']['filename'] = '../../outside.csv'
        with TemporaryDirectory() as folder:
            destination = Path(folder)
            existing = destination / 'timestamps.csv'
            existing.write_text('keep me')
            first = restore_embedded_inputs(config, destination)
            second = restore_embedded_inputs(config, destination)
            self.assertEqual(existing.read_text(), 'keep me')
            self.assertEqual(first['timestamps'].name, 'timestamps_2.csv')
            self.assertEqual(second['timestamps'].name, 'timestamps_3.csv')
            self.assertEqual(first['timestamps'].parent, destination.resolve())

    def test_corrupt_input_is_reported_before_any_csv_is_written(self):
        config = self.example()
        config.values['inputs']['timestamps']['sha256'] = 'bad'
        with TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, 'timestamps.*integrity check'):
                restore_embedded_inputs(config, folder)
            self.assertEqual(list(Path(folder).iterdir()), [])

    def test_no_embedded_data_explains_how_to_enable_restoration(self):
        config = AnalysisConfig({'inputs': {'timestamps': 'missing.csv'}}, Path.cwd())
        with TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, 'Include all input data in JSON'):
                restore_embedded_inputs(config, folder)

    def test_gui_restore_uses_folder_choice_without_running_analysis(self):
        config = self.example()
        store = _embedded_form_inputs(config.values)
        app = create_app()
        callback = app.callback_map['restore-analysis-message.children']['callback'].__wrapped__
        enabled = app.callback_map['restore-analysis-inputs.disabled']['callback'].__wrapped__
        self.assertTrue(enabled({}))
        self.assertFalse(enabled(store))
        with TemporaryDirectory() as folder, \
             patch('autoqy_core.tools.analysis_gui._choose_folder', return_value=folder) as chooser:
            message = callback(1, store, 'base-folder')
            chooser.assert_called_once_with('base-folder')
            self.assertIn('Restored 5 CSV files', message)
            self.assertEqual(len(list(Path(folder).glob('*.csv'))), 5)
        self.assertIn('Load a JSON', callback(1, {}, 'base-folder'))
        with patch('autoqy_core.tools.analysis_gui._choose_folder', return_value=None):
            self.assertIs(callback(1, store, 'base-folder'), no_update)
        store['timestamps']['input']['sha256'] = 'bad'
        with TemporaryDirectory() as folder, \
             patch('autoqy_core.tools.analysis_gui._choose_folder', return_value=folder):
            message = callback(1, store, 'base-folder')
            self.assertIn('Could not restore input files', message)
            self.assertEqual(list(Path(folder).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
