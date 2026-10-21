"""Consistent values and precision across human-readable yield reports."""
from types import SimpleNamespace
import unittest
import numpy as np
from autoqy_core.output import format_value_uncertainty, formatted_quantum_yields, quantum_yield_values
from autoqy_core.tools.analysis_gui import _render_comparison


class ResultFormattingTests(unittest.TestCase):
    def test_error_significant_digits_determine_mean_decimal_place(self):
        for value,error,expected in (
            (50.7891677,.7699928,("50.8","0.8")),
            (75.922899,4.8728445,("76","5")),
            (50.86078,1.69507,("50.9","1.7")),
            (79.0718,21.9018,("79","22")),
            (.123456,.012345,("0.123","0.012")),
            (1234,123,("1230","120")),
            (4.821,.099,("4.8","0.1")),
        ):
            with self.subTest(value=value,error=error):
                self.assertEqual(format_value_uncertainty(value,error),expected)

    def test_mc_reporting_uses_mean_and_keeps_the_original_fit(self):
        result=SimpleNamespace(
            yield_fit=SimpleNamespace(values=np.array([.5086078,.790718])),
            yield_errors=np.array([.007699928,.048728445]),
            epsilon_uncertainty=SimpleNamespace(method='monte_carlo',sample_mean=np.array([.507891677,.75922899])))
        self.assertEqual(formatted_quantum_yields(result),[("50.8","0.8"),("76","5")])
        np.testing.assert_array_equal(result.yield_fit.values,[.5086078,.790718])
        result.epsilon_uncertainty=None
        np.testing.assert_array_equal(quantum_yield_values(result),result.yield_fit.values)

    def test_terminal_uses_the_same_mc_mean_and_rounding(self):
        from contextlib import redirect_stdout
        from io import StringIO
        from unittest.mock import patch
        from autoqy_core.cli import main
        result=SimpleNamespace(
            yield_fit=SimpleNamespace(values=np.array([.5086078,.790718])),
            yield_errors=np.array([.007699928,.048728445]),
            epsilon_uncertainty=SimpleNamespace(method='monte_carlo',sample_mean=np.array([.507891677,.75922899])))
        stream=StringIO()
        with patch('autoqy_core.cli.load_config'), \
             patch('autoqy_core.cli.run_analysis',return_value=SimpleNamespace(result=result,files=())), \
             redirect_stdout(stream):
            self.assertEqual(main(['run','analysis.json']),0)
        self.assertEqual(stream.getvalue().splitlines(),
                         ['Phi_R->P: 50.8 +/- 0.8%','Phi_P->R: 76 +/- 5%'])

    def test_comparison_reports_mean_and_error_instead_of_verbose_intervals(self):
        from dash import html
        row={'method':'ode_absorbance','values':np.array([50.7891677,75.922899]),
             'errors':np.array([.7699928,4.8728445]),'interval':([49,65],[53,85],.95),
             'fraction_rmse':.001,'absorbance_rmse':.001,'diagnostics':[]}
        def strings(root):
            if isinstance(root,str): yield root
            elif isinstance(root,(list,tuple)):
                for child in root: yield from strings(child)
            elif getattr(root,'children',None) is not None:
                yield from strings(root.children)
        text=' '.join(strings(_render_comparison(html,[row],'ode_absorbance')))
        self.assertIn('50.8 ± 0.8%',text)
        self.assertIn('76 ± 5%',text)
        self.assertNotIn('Nominal',text)
        self.assertNotIn('95%',text)


if __name__=='__main__': unittest.main()
