from types import SimpleNamespace
import unittest

from autoqy_core.tools.analysis_gui import (
    _render_checks, _yield_reliability_check, create_app,
)


def components(root):
    if isinstance(root, (list, tuple)):
        for item in root:
            yield from components(item)
    else:
        yield root
        children = getattr(root, 'children', None)
        if children is not None:
            yield from components(children)


class WarningTests(unittest.TestCase):
    def test_fit_checks_preserve_thresholds_and_explain_the_problem(self):
        cases = [
            (1e6, True, (False, False), 'ok', 'passed'),
            (1e6 + 1, True, (False, False), 'warning', 'measurement errors'),
            (1e8, True, (False, False), 'warning', 'measurement errors'),
            (1e8 + 1, True, (False, False), 'stop', 'cannot reliably separate'),
            (float('inf'), True, (False, False), 'stop', 'cannot reliably separate'),
            (float('nan'), True, (False, False), 'stop', 'cannot reliably separate'),
            (1, False, (False, False), 'stop', 'did not finish'),
            (1, True, (True, False), 'warning', 'minimum or maximum'),
            (1e9, True, (True, False), 'stop', 'minimum or maximum'),
        ]
        for condition, success, bounds, level, phrase in cases:
            with self.subTest(condition=condition, success=success, bounds=bounds):
                check = _yield_reliability_check(SimpleNamespace(
                    jacobian_condition=condition, optimizer_success=success, active_bounds=bounds))
                self.assertEqual(check['level'], level)
                self.assertIn(phrase, check['body'])
                self.assertNotIn('Jacobian', check['body'])
                self.assertIn('Jacobian', check['details'])

    def test_explanation_is_visible_and_technical_details_are_collapsed(self):
        from dash import html
        check = _yield_reliability_check(SimpleNamespace(
            jacobian_condition=1e9, optimizer_success=True, active_bounds=(False, False)))
        rendered = list(components(_render_checks(html, [check])))
        message = next(c for c in rendered if getattr(c, 'className', '') == 'diagnostic-message')
        self.assertEqual(message.children, check['body'])
        details = next(c for c in rendered if isinstance(c, html.Details))
        self.assertFalse(getattr(details, 'open', False))

    def test_analysis_labels_and_input_use_metastable_product_percent(self):
        app = create_app()
        layout = list(components(app.layout))
        text = ' '.join(c.children for c in layout if isinstance(getattr(c, 'children', None), str))
        self.assertIn('Expected metastable photoproduct at PSS (%)', text)
        self.assertNotIn('Expected reactant at PSS', text)
        field = next(c for c in layout if getattr(c, 'id', None) == {
            'type': 'analysis-field', 'name': 'expected_pss_product_percent'})
        self.assertEqual((field.min, field.max, field.value), (0, 100, None))
        self.assertNotIn('name": "expected_pss_reactant_percent', text)


if __name__ == '__main__':
    unittest.main()
