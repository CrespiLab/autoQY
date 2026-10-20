"""GUI composition uses product percent; schema v1 retains reactant percent."""

import json
from pathlib import Path
import unittest
from tempfile import TemporaryDirectory

import numpy as np

from autoqy_core.tools.analysis_gui import (
    _complement_percent, _configuration, _form_values,
)


class CompositionTests(unittest.TestCase):
    def document(self):
        path = (Path(__file__).parents[1] / 'ExampleData'
                / 'Example-4_395nm-EpsilonError/generic_inputs/analysis.json')
        return json.loads(path.read_text(encoding='utf-8'))

    def test_product_percent_round_trips_through_legacy_schema(self):
        for product in (0, 23, 77, 100, 76.25, None):
            with self.subTest(product=product):
                values = _form_values(self.document())
                values['expected_pss_product_percent'] = product
                saved = _configuration(values)
                expected = None if product is None else 100 - product
                self.assertEqual(saved['analysis']['expected_pss_reactant_percent'], expected)
                self.assertNotIn('expected_pss_product_percent', saved['analysis'])
                loaded = _form_values(saved)
                self.assertEqual(loaded['expected_pss_product_percent'], product)
                self.assertNotIn('expected_pss_reactant_percent', loaded)

    def test_existing_configs_display_complement_without_changing_the_document(self):
        document = self.document()
        document['analysis']['expected_pss_reactant_percent'] = 23
        self.assertEqual(_form_values(document)['expected_pss_product_percent'], 77)
        self.assertEqual(document['analysis']['expected_pss_reactant_percent'], 23)
        document['analysis'].pop('expected_pss_reactant_percent')
        self.assertIsNone(_form_values(document)['expected_pss_product_percent'])

    def test_invalid_compositions_are_rejected_server_side(self):
        for value in (-1, 101, float('nan'), float('inf'), True, 'invalid'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _complement_percent(value)

    def test_form_preserves_all_other_configuration_values(self):
        document = self.document()
        before = _configuration(_form_values(document))
        values = _form_values(document)
        values['expected_pss_product_percent'] = 77
        after = _configuration(values)
        before['analysis']['expected_pss_reactant_percent'] = 23
        self.assertEqual(before, after)

    def test_example_fit_and_pss_diagnostic_use_the_same_composition(self):
        from autoqy_core.config import AnalysisConfig
        from autoqy_core.runner import run_analysis
        from autoqy_core.tools.analysis_gui import _fit_diagnostic_checks

        document = self.document()
        base = (Path(__file__).parents[1] / 'ExampleData'
                / 'Example-4_395nm-EpsilonError/generic_inputs')
        document['uncertainty']['epsilon']['method'] = 'none'
        for key in ('write_text', 'write_figures', 'write_json', 'write_config', 'write_detailed_data'):
            document['outputs'][key] = False
        with TemporaryDirectory() as temporary:
            original = run_analysis(AnalysisConfig(document, base), temporary)
            values = _form_values(document)
            values['expected_pss_product_percent'] = 77
            converted = _configuration(values)
            updated = run_analysis(AnalysisConfig(converted, base), temporary)
        np.testing.assert_allclose(original.result.yield_fit.values, updated.result.yield_fit.values)
        np.testing.assert_allclose(original.result.extrapolated_pss, updated.result.extrapolated_pss)
        checks = _fit_diagnostic_checks(updated.result, updated.data, converted)
        pss = next(c for c in checks if c['title'] == 'Expected versus fitted PSS')
        self.assertIn('photoproduct expected 77.00%', pss['body'])
        # The expected composition remains a diagnostic, never an optimizer constraint.
        self.assertEqual(converted['analysis']['expected_pss_reactant_percent'], 23)


if __name__ == '__main__':
    unittest.main()
