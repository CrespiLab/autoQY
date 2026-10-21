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
    _epsilon_reliability_check,
    _fit_diagnostic_checks,
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
        bound_yields=np.array([[0.12, 0.08]]),
        bound_jacobian_conditions=np.array([20.0]),
        bound_active_bounds=np.array([[False, False]]),
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


class EpsilonReliabilityTests(unittest.TestCase):
    def test_precise_refits_with_large_shifts_warn_in_analysis_and_comparison(self):
        result, data = _result_and_data()
        result.yield_fit.values = np.array([0.348, 0.422])
        result.yield_errors = np.array([0.655, 0.033])
        u = result.epsilon_uncertainty
        u.bound_labels = ("reactant lower / product nominal", "nominal", "reactant upper / product nominal")
        u.bound_yields = np.array([[1.0, 0.395], [0.348, 0.422], [0.156, 0.424]])
        u.bound_optimizer_errors = np.full((3, 2), 0.003)
        u.bound_active_bounds = np.array([[True, False], [False, False], [False, False]])
        u.bound_jacobian_conditions = np.array([86.0, 34.0, 25.0])
        checks = _fit_diagnostic_checks(result, data, {"analysis": {}})
        sweep = next(c for c in checks if c["title"] == "Uncertainty when reference spectra are varied")
        self.assertEqual(sweep["level"], "warning")
        self.assertIn("65.5", sweep["body"])
        self.assertIn("15.6–100%", sweep["body"])
        self.assertIn("truncated", sweep["body"])
        self.assertIn("reactant lower / product nominal", sweep["details"])
        overall = next(c for c in checks if c["title"] == "Can I trust the quantum yields?")
        self.assertEqual(overall["level"], "warning")
        self.assertEqual(_comparison_diagnostic_summary(checks)[0], "warning")

    def test_combined_error_warns_even_without_active_limits(self):
        result, _ = _result_and_data()
        result.yield_errors = np.array([0.3, 0.01])
        result.epsilon_uncertainty.bound_optimizer_errors[:] = 0.001
        check = _epsilon_reliability_check(result)
        self.assertEqual(check["level"], "warning")
        self.assertIn("30", check["body"])

    def test_small_absolute_but_large_relative_error_warns(self):
        result, _ = _result_and_data()
        result.yield_fit.values = np.array([0.01, 0.08])
        result.yield_errors = np.array([0.02, 0.005])
        self.assertEqual(_epsilon_reliability_check(result)["level"], "warning")

    def test_each_sweep_failure_mode_prevents_a_green_check(self):
        cases = (
            ("bound_active_bounds", np.array([[False, True]]), "warning", "minimum or maximum"),
            ("bound_jacobian_conditions", np.array([1e7]), "warning", "sensitive"),
            ("bound_jacobian_conditions", np.array([1e9]), "stop", "cannot reliably separate"),
            ("bound_jacobian_conditions", np.array([np.nan]), "stop", "cannot reliably separate"),
            ("bound_optimizer_success", (False,), "stop", "did not finish"),
            ("bound_optimizer_errors", np.array([[np.nan, 0.001]]), "stop", "not finite"),
            ("bound_yields", np.array([[np.nan, 0.08]]), "stop", "not finite"),
        )
        for name, value, level, phrase in cases:
            with self.subTest(name=name, value=value):
                result, _ = _result_and_data()
                setattr(result.epsilon_uncertainty, name, value)
                check = _epsilon_reliability_check(result)
                self.assertEqual(check["level"], level)
                self.assertIn(phrase, check["body"])

    def test_large_and_nonfinite_combined_errors_are_red(self):
        for error in (1.0, np.inf, np.nan):
            with self.subTest(error=error):
                result, _ = _result_and_data()
                result.yield_errors[0] = error
                self.assertEqual(_epsilon_reliability_check(result)["level"], "stop")

    def test_small_stable_ranges_stay_green_and_disabled_epsilon_is_unchanged(self):
        result, _ = _result_and_data()
        self.assertEqual(_epsilon_reliability_check(result)["level"], "ok")
        result.epsilon_uncertainty = None
        self.assertIsNone(_epsilon_reliability_check(result))

    def test_epsilon_propagation_retains_unsuccessful_refit_status(self):
        from autoqy_core.epsilon_uncertainty import EpsilonEnvelope, run_with_epsilon_uncertainty
        from autoqy_core.pipeline import AnalysisInput, AnalysisResult

        result, data = _result_and_data()
        inputs = AnalysisInput(
            wavelengths=result.wavelengths, absorbance=result.absorbance,
            timestamps=data.timestamps,
            epsilon_r=(result.wavelengths, result.epsilon_r),
            epsilon_p=(result.wavelengths, result.epsilon_p),
            led=data.led, power_mw=1, power_error_mw=0, volume_ml=1,
        )
        def envelope(values):
            return EpsilonEnvelope(result.wavelengths, values, values * .9, values * 1.1,
                                   "epsilon_spectra", "sd", "test.csv")
        calls = []
        def fit(_):
            outcome = deepcopy(result.yield_fit)
            outcome.optimizer_success = len(calls) != 1
            calls.append(outcome)
            return AnalysisResult(result.concentration_fit, outcome, result.yield_errors,
                                  result.wavelengths, result.absorbance,
                                  result.epsilon_r, result.epsilon_p, "ode_absorbance",
                                  result.extrapolated_pss)
        with patch("autoqy_core.pipeline.run_analysis_pipeline", side_effect=fit):
            output = run_with_epsilon_uncertainty(inputs, envelope(result.epsilon_r), envelope(result.epsilon_p))
        self.assertEqual(len(calls), 9)
        u = output.epsilon_uncertainty
        self.assertEqual(len(u.bound_optimizer_success), 9)
        failed = u.bound_optimizer_success.index(False)
        self.assertEqual(u.bound_labels[failed], "reactant lower / product lower")
        self.assertEqual(_epsilon_reliability_check(output)["level"], "stop")


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
            annotation = spectra.texts[0].get_text()
            self.assertIn("±", annotation)
            self.assertNotIn("(", annotation)
            self.assertIsNone(figure._suptitle)

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
