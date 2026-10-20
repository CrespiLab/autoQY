"""Regression checks for comparison uncertainty and exported result presentation."""

from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np

from autoqy_core.config import AnalysisConfig
from autoqy_core.plotting import write_figure
from autoqy_core.tools.analysis_gui import (
    _compare_fit_methods,
    _comparison_diagnostic_summary,
    _interactive_figures,
    _method_label,
)

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    go = make_subplots = None


def _result_and_data():
    concentrations = np.array([[1.0, 0.0], [0.7, 0.3], [0.5, 0.5]]) * 1e-4
    epsilon = np.array([[1000.0, 400.0, 100.0], [100.0, 400.0, 1000.0]])
    uncertainty = SimpleNamespace(
        concentration_data_minimum=concentrations * 0.9,
        concentration_data_maximum=concentrations * 1.1,
        concentration_fit_minimum=concentrations * 0.9,
        concentration_fit_maximum=concentrations * 1.1,
        fraction_residual_minimum=np.full(3, -0.01),
        fraction_residual_maximum=np.full(3, 0.01),
        absorbance_residual_rmse_minimum=0.001,
        absorbance_residual_rmse_maximum=0.002,
        bound_labels=("lower references",),
        bound_optimizer_errors=np.array([[0.02, 0.03]]),
    )
    result = SimpleNamespace(
        concentration_fit=SimpleNamespace(
            concentrations=concentrations,
            fractions=concentrations / concentrations.sum(axis=1, keepdims=True),
        ),
        yield_fit=SimpleNamespace(
            concentrations=concentrations.copy(),
            values=np.array([0.12, 0.08]),
            absorbance_correction=None,
            optimizer_success=True,
            jacobian_condition=20.0,
            active_bounds=(False, False),
        ),
        epsilon_r=epsilon[0], epsilon_p=epsilon[1],
        wavelengths=np.array([400.0, 450.0, 500.0]),
        absorbance=(concentrations @ epsilon).T,
        yield_errors=np.array([0.04, 0.05]),
        epsilon_uncertainty=uncertainty,
        extrapolated_pss=np.array([0.4, 0.6]),
    )
    data = SimpleNamespace(
        timestamps=np.array([0.0, 10.0, 20.0]), path_length_cm=1.0,
        photon_flux_mol_s=1e-9,
        led=(result.wavelengths, np.array([0.0, 1.0, 0.0])),
        wavelength_limits=(400.0, 500.0),
    )
    return result, data


class MethodComparisonTests(unittest.TestCase):
    def test_comparison_preserves_selected_epsilon_and_reports_combined_errors(self):
        epsilon_settings = {"method": "deterministic_extremes", "error_metric": "sem"}
        document = {
            "analysis": {"expected_pss_reactant_percent": None},
            "fit": {"method": "concentrations"},
            "uncertainty": {"epsilon": epsilon_settings},
            "outputs": {name: True for name in (
                "write_text", "write_figures", "write_json", "write_config",
                "write_detailed_data",
            )},
        }
        original = deepcopy(document)
        config = AnalysisConfig(document, Path.cwd())
        calls = []

        def analyze(candidate, output_directory):
            calls.append(candidate.values["fit"]["method"])
            self.assertEqual(candidate.values["uncertainty"]["epsilon"], epsilon_settings)
            self.assertFalse(any(candidate.values["outputs"].values()))
            self.assertTrue(Path(output_directory).is_dir())
            self.assertEqual(list(Path(output_directory).iterdir()), [])
            result, data = _result_and_data()
            return SimpleNamespace(result=result, data=data)

        with patch("autoqy_core.runner.run_analysis", side_effect=analyze):
            rows = _compare_fit_methods(config)

        self.assertEqual(set(calls), {"concentrations", "ode_absorbance", "emission"})
        self.assertEqual(len(calls), 3)
        for row in rows:
            np.testing.assert_allclose(row["errors"], [4.0, 5.0])
            self.assertIn(
                "Uncertainty when reference spectra are varied",
                [item["title"] for item in row["diagnostics"]],
            )
        self.assertEqual(document, original)

    def test_comparison_flags_report_highest_diagnostic_level(self):
        checks = [
            {"level": "warning", "title": "Initial product"},
            {"level": "stop", "title": "Quantum-yield uncertainty"},
        ]
        self.assertEqual(
            _comparison_diagnostic_summary(checks),
            ("stop", "Red flag: Quantum-yield uncertainty"),
        )
        self.assertEqual(_method_label("ode_absorbance"), "Full spectrum ODE")


class ResultPresentationTests(unittest.TestCase):
    def test_png_and_svg_export_simplify_labels_and_keep_epsilon_ranges(self):
        result, data = _result_and_data()
        with TemporaryDirectory() as temporary:
            with patch("autoqy_core.plotting.plt.close", wraps=plt.close) as closed:
                png, svg = write_figure(Path(temporary) / "result", result, data)
            self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertIn("<svg", svg.read_text(encoding="utf-8"))
            figure = closed.call_args.args[0]
            concentration, spectra, residual, heatmap = figure.axes[:4]
            self.assertEqual(concentration.get_title(), "Concentrations with ε ranges")
            labels = concentration.get_legend_handles_labels()[1]
            self.assertEqual(labels, ["Reactant data", "Product data", "Reactant fit", "Product fit"])
            for line in concentration.lines:
                if line.get_label().endswith(" fit"):
                    self.assertGreaterEqual(line.get_linewidth(), 2)
                    self.assertGreater(line.get_zorder(), 3)
            self.assertIn("Residual", residual.get_legend_handles_labels()[1])
            self.assertEqual(heatmap.get_title(), "Absorbance residuals")
            self.assertEqual(len(spectra.texts), 1)
            self.assertIn("±", spectra.texts[0].get_text())

    @unittest.skipUnless(go is not None, "Plotly is not installed")
    def test_gui_draws_fitted_lines_above_all_measured_points(self):
        result, data = _result_and_data()
        with patch("autoqy_core.tools.analysis_gui._spectra_led_figure", return_value=go.Figure()):
            concentration, fraction, _, _, heatmap = _interactive_figures(
                go, make_subplots, result, data, 99, 450,
            )
        names = [trace.name for trace in concentration.data]
        last_data = max(names.index("Reactant data"), names.index("Product data"))
        for name in ("Reactant fit", "Product fit"):
            trace_index = names.index(name)
            self.assertGreater(trace_index, last_data)
            self.assertEqual(concentration.data[trace_index].line.width, 3)
        self.assertEqual(fraction.data[-1].name, "Residual")
        self.assertEqual(heatmap.layout.title.text, "Absorbance residuals")


if __name__ == "__main__":
    unittest.main()
