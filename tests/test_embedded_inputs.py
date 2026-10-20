import base64
from copy import deepcopy
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np
from dash import ctx

from autoqy_core.config import (AnalysisConfig, ConfigError, load_config,
                               self_contained_document, validate_config)
from autoqy_core.runner import run_analysis
from autoqy_core.tools.analysis_gui import (
    _configuration, _embedded_form_inputs, _form_values, _portable_document, create_app,
)


EXAMPLE = (Path(__file__).parents[1] / 'ExampleData'
           / 'Example-4_395nm-EpsilonError/generic_inputs')


class EmbeddedInputsTests(unittest.TestCase):
    def example(self):
        values = json.loads((EXAMPLE / 'analysis.json').read_text(encoding='utf-8'))
        values['uncertainty']['epsilon']['method'] = 'none'
        for name in ('write_text', 'write_figures', 'write_json', 'write_config', 'write_detailed_data'):
            values['outputs'][name] = False
        return AnalysisConfig(values, EXAMPLE)

    def test_shared_json_runs_without_any_original_files_and_retains_epsilon_errors(self):
        with TemporaryDirectory() as temporary:
            base = Path(temporary)
            raw = base / 'raw'
            raw.mkdir()
            # Copy contents without OneDrive directory/read-only attributes.
            for name in ('measurement_spectra.csv', 'reactant_absorptivity.csv',
                         'product_absorptivity.csv', 'led_emission.csv', 'timestamps.csv'):
                shutil.copyfile(EXAMPLE / name, raw / name)
            config = self.example()
            config = AnalysisConfig(config.values, raw)
            original = run_analysis(config, base / 'original')
            document = self_contained_document(config)
            for name in ('measurement_spectra', 'reactant_absorptivity', 'product_absorptivity',
                         'led_emission', 'timestamps'):
                self.assertEqual(base64.b64decode(document['inputs'][name]['data']),
                                 config.input_path(name).read_bytes())
            shared = base / 'recipient'
            shared.mkdir()
            saved = shared / 'experiment.json'
            saved.write_text(json.dumps(document), encoding='utf-8')
            shutil.rmtree(raw)
            loaded = load_config(saved)
            updated = run_analysis(loaded, shared / 'results')
            np.testing.assert_array_equal(original.data.absorbance, updated.data.absorbance)
            np.testing.assert_allclose(original.result.yield_fit.values, updated.result.yield_fit.values)
            # Exact source bytes retain every uncertainty column, even with propagation off.
            from autoqy_core.epsilon_uncertainty import load_epsilon_envelope
            envelope = load_epsilon_envelope(loaded.input_path('product_absorptivity'), 'sd')
            self.assertGreater(len(envelope.lower), 0)

    def test_inline_inputs_survive_gui_reload_resave_and_path_rebasing(self):
        document = self_contained_document(self.example())
        values = _form_values(document)
        store = _embedded_form_inputs(document)
        rebuilt = _configuration(values, store)
        self.assertEqual(rebuilt['schema_version'], 2)
        self.assertEqual(rebuilt['inputs'], document['inputs'])
        with TemporaryDirectory() as folder:
            rebased = _portable_document(rebuilt, Path(folder), EXAMPLE)
            self.assertEqual(rebased['inputs'], document['inputs'])
        values['measurement_spectra'] = 'replacement.csv'
        rebuilt = _configuration(values, store)
        self.assertEqual(rebuilt['inputs']['measurement_spectra'], 'replacement.csv')

    def test_result_and_snapshot_json_are_both_rerunnable(self):
        config = self.example()
        config.values['outputs'].update(embed_inputs=True, write_json=True, write_config=True)
        with TemporaryDirectory() as folder:
            result = run_analysis(config, folder)
            for path in result.files:
                loaded = load_config(path)
                self.assertEqual(loaded.values['schema_version'], 2)
                self.assertIsInstance(loaded.values['inputs']['timestamps'], dict)

    def test_invalid_or_modified_data_is_rejected(self):
        document = self_contained_document(self.example())
        for mutation in ({'data': '!'}, {'sha256': 'bad'}, {'encoding': 'pickle'}, {'filename': ''}):
            bad = deepcopy(document)
            bad['inputs']['timestamps'].update(mutation)
            with self.subTest(mutation=mutation), self.assertRaises(ConfigError):
                validate_config(AnalysisConfig(bad, Path.cwd()))
        document['schema_version'] = 1
        with self.assertRaisesRegex(ConfigError, 'schema_version 2'):
            validate_config(AnalysisConfig(document, Path.cwd()))

    def test_arbitrary_binary_bytes_roundtrip_and_cannot_escape_private_directory(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'data.dat'
            payload = bytes(range(256))
            path.write_bytes(payload)
            config = self.example()
            config.values['inputs']['measurement_spectra'] = str(path)
            document = self_contained_document(config)
            document['inputs']['measurement_spectra']['filename'] = '../../outside.dat'
            loaded = AnalysisConfig(document, Path(folder))
            extracted = loaded.input_path('measurement_spectra')
            self.assertEqual(extracted.read_bytes(), payload)
            self.assertEqual(extracted.parent, Path(loaded._embedded_files.name))
            self.assertEqual(extracted, loaded.input_path('measurement_spectra'))

    def test_actinometer_bundle_does_not_require_an_led_file(self):
        config = self.example()
        config.values['experiment'].update(irradiation_source='chemical_actinometer',
                                           photon_flux_mol_s=1e-9)
        config.values['inputs']['led_emission'] = 'missing-led.csv'
        document = self_contained_document(config)
        self.assertNotIn('led_emission', document['inputs'])
        validate_config(AnalysisConfig(document, Path.cwd()))

    def test_legacy_files_load_and_linked_saving_remains_available(self):
        config = load_config(EXAMPLE / 'analysis.json')
        self.assertEqual(config.values['schema_version'], 1)
        self.assertIsInstance(config.values['inputs']['measurement_spectra'], str)
        values = _form_values(config.values)
        values['embed_inputs'] = []
        document = _configuration(values)
        self.assertFalse(document['outputs']['embed_inputs'])
        self.assertEqual(document['schema_version'], 1)

    def test_gui_load_callback_preserves_embedded_data(self):
        document = self_contained_document(self.example())
        app = create_app()
        callback = next(item['callback'].__wrapped__ for key, item in app.callback_map.items()
                        if 'load-analysis-message.children' in key)
        ids = [{'type': 'analysis-field', 'name': name} for name in _form_values(document)]
        values = list(_form_values(document).values())
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'shared.json'
            path.write_text(json.dumps(document))
            with patch.object(type(ctx), 'triggered_id', property(lambda _: 'load-analysis-json')), \
                 patch('autoqy_core.tools.analysis_gui._choose_file', return_value=str(path)):
                fields, message, embedded = callback(1, 0, [], values, ids, {})
        self.assertIn('Loaded', message)
        rebuilt = _configuration(dict(zip((item['name'] for item in ids), fields)), embedded)
        self.assertEqual(rebuilt['inputs'], document['inputs'])

    def test_gui_save_button_writes_a_self_contained_json(self):
        config = self.example()
        app = create_app()
        callback = next(item['callback'].__wrapped__ for key, item in app.callback_map.items()
                        if 'analysis-action-message.children' in key)
        values = _form_values(config.values)
        values['config_folder'] = str(EXAMPLE)
        ids = [{'type': 'analysis-field', 'name': name} for name in values]
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'shared.json'
            with patch.object(type(ctx), 'triggered_id', property(lambda _: 'save-analysis-json')), \
                 patch('autoqy_core.tools.analysis_gui._choose_save_json', return_value=str(path)), \
                 patch('autoqy_core.tools.analysis_gui._input_checks', return_value=[]):
                result = callback(1, 0, 0, list(values.values()), ids, {})
            self.assertIn('Saved and validated', result[0])
            loaded = load_config(path)
            self.assertEqual(loaded.values['schema_version'], 2)
            self.assertTrue(all(isinstance(loaded.values['inputs'][name], dict)
                                for name in ('measurement_spectra', 'timestamps', 'led_emission')))


if __name__ == '__main__':
    unittest.main()
