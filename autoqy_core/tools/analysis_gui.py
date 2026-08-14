"""Local browser GUI for building and running AutoQY analyses."""

from argparse import ArgumentParser
from copy import deepcopy
import json
from pathlib import Path
import socket
import subprocess
import sys
from threading import Lock, Thread, Timer
from tempfile import TemporaryDirectory
import time
import warnings
import webbrowser

import numpy as np

from ..config import AnalysisConfig, input_format, load_config, validate_config
from ..epsilon_uncertainty import load_epsilon_nominal
from ..io import load_spectra, load_spectrum, load_timestamps
from ..output import result_summary
from ..runner import run_analysis
from ..spectra import process_led


SPECTRAL_INPUTS = (
    ("measurement_spectra", "Measurement spectra"),
    ("reactant_absorptivity", "Reactant molar absorptivity"),
    ("product_absorptivity", "Product molar absorptivity"),
    ("led_emission", "LED emission"),
)
TIMESTAMP_INPUT = ("timestamps", "Irradiation timestamps")


def create_app():
    try:
        from dash import ALL, Dash, Input, Output, State, ctx, dcc, html
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as error:
        raise RuntimeError(
            "The analysis GUI requires the 'power-gui' optional dependencies"
        ) from error

    assets = Path(__file__).parents[1] / "assets"
    app = Dash(__name__, assets_folder=str(assets), suppress_callback_exceptions=True)
    app.title = "AutoQY Analysis"
    window_state = {"close_requested_at": None}
    window_lock = Lock()
    app.server.config["AUTOQY_WINDOW_STATE"] = (window_state, window_lock)

    @app.server.post("/_autoqy_analysis_heartbeat")
    def autoqy_heartbeat():
        with window_lock:
            window_state["close_requested_at"] = None
        return "", 204

    @app.server.post("/_autoqy_analysis_window_closed")
    def autoqy_window_closed():
        with window_lock:
            window_state["close_requested_at"] = time.monotonic()
        return "", 204

    app.index_string = """<!DOCTYPE html>
<html><head>{%metas%}<title>{%title%}</title>{%favicon%}{%css%}</head>
<body>{%app_entry%}<footer>{%config%}{%scripts%}{%renderer%}</footer>
<script>
(() => {
  const heartbeat = () => fetch('/_autoqy_analysis_heartbeat', {method: 'POST', keepalive: true});
  heartbeat();
  window.addEventListener('pageshow', heartbeat);
  window.addEventListener('pagehide', () => navigator.sendBeacon('/_autoqy_analysis_window_closed'));
})();
</script></body></html>"""

    def field(name, **kwargs):
        return dcc.Input(id={"type": "analysis-field", "name": name}, **kwargs)

    def number(name, value, **kwargs):
        return field(name, type="number", value=value, **kwargs)

    def toggle(name, label, enabled=True):
        return dcc.Checklist(
            id={"type": "analysis-field", "name": name},
            options=[{"label": label, "value": "on"}],
            value=["on"] if enabled else [], className="toggle-control",
        )

    def dropdown(name, options, value, **kwargs):
        return dcc.Dropdown(
            id={"type": "analysis-field", "name": name},
            options=[{"label": label, "value": option} for option, label in options],
            value=value, clearable=False, **kwargs,
        )

    def pair(label_a, component_a, label_b, component_b):
        return html.Div(className="input-row", children=[
            html.Div([html.Label(label_a), component_a]),
            html.Div([html.Label(label_b), component_b]),
        ])

    def file_card(name, label, timestamp=False):
        formats = (("ahk_csv", "AHK CSV"), ("simple_csv", "Simple CSV"),
                   ("generic_delimited", "Generic delimited")) if timestamp else (
                       ("spectragryph_tsv", "SpectraGryph / AutoQY TSV"),
                       ("generic_delimited", "Generic delimited"),
                   )
        default = "ahk_csv" if timestamp else "spectragryph_tsv"
        return html.Div(className="analysis-file-card", children=[
            html.Label(label),
            html.Div(className="path-row", children=[
                field(name, type="text", value="", placeholder="Choose a local file"),
                html.Button("Browse", id={"type": "browse-analysis-file", "name": name},
                            n_clicks=0, className="button button-secondary compact-button"),
            ]),
            dropdown(f"format_{name}", formats, default),
            html.Details(className="format-options", children=[
                html.Summary("Generic-format options"),
                pair("Delimiter", field(f"delimiter_{name}", type="text", value=","),
                     "Header row", toggle(f"header_{name}", "Present", True)),
            ]),
        ])

    default_folder = str(Path.cwd())
    empty = _empty_figure(go, "Run a validated analysis to display this plot")
    app.layout = html.Div(className="app-shell analysis-app", children=[
        html.Header(className="app-header", children=[
            html.Div([
                html.P("AUTOQY CORE", className="eyebrow"),
                html.H1("Quantum-yield analysis"),
                html.P("Build, validate, run, and inspect a reproducible analysis.json.",
                       className="subtitle"),
            ]),
            html.Div("Local session", className="local-badge"),
        ]),
        html.Main(className="workspace analysis-workspace", children=[
            html.Aside(className="control-column analysis-controls", children=[
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("1 · Project", className="step-label"),
                                  html.Span("JSON and tools")]),
                    html.Button("Load existing analysis.json", id="load-analysis-json",
                                n_clicks=0, className="button button-secondary"),
                    html.Label("JSON base folder"),
                    html.Div(className="path-row", children=[
                        field("config_folder", type="text", value=default_folder),
                        html.Button("Choose", id="choose-analysis-folder", n_clicks=0,
                                    className="button button-secondary compact-button"),
                    ]),
                    html.P("Relative input and output paths are resolved from this folder.",
                           className="helper-text"),
                    html.Button("Open Spectral Treatment", id="launch-spectral-treatment",
                                n_clicks=0, className="button button-accent"),
                    html.Div(id="tool-message", className="message"),
                    html.Div(id="load-analysis-message", className="message"),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("2 · Identity and data", className="step-label"),
                                  html.Span("Experiment files")]),
                    pair("Analysis ID", field("analysis_id", type="text", value="new_analysis"),
                         "Output stem", field("output_stem", type="text", value="AutoQY_results")),
                    pair("Reactant name", field("reactant_name", type="text", value="reactant"),
                         "Product name", field("product_name", type="text", value="product")),
                    *[file_card(name, label) for name, label in SPECTRAL_INPUTS],
                    file_card(*TIMESTAMP_INPUT, timestamp=True),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("3 · Experiment", className="step-label"),
                                  html.Span("Physical parameters")]),
                    pair("Sample volume", html.Div(className="volume-input-row", children=[
                             number("volume_value", 3000, min=0, step="any"),
                             dropdown("volume_unit", (("ul", "µL"), ("ml", "mL")), "ul"),
                         ]),
                         "Path length (cm)", number("path_length_cm", 1.0, min=0, step="any")),
                    pair("Power (mW)", number("power_mw", 1.0, min=0, step="any"),
                         "Power error (mW)", number("power_error_mw", 0.0, min=0, step="any")),
                    pair("Irradiation wavelength (nm)",
                         number("irradiation_wavelength_nm", 455, min=0, step="any"),
                         "Thermal back-reaction (s⁻¹)",
                         number("thermal_rate", 0.0, min=0, step="any")),
                    html.Div(id="physical-parameter-preview", className="message"),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("4 · Processing", className="step-label"),
                                  html.Span("Wavelength and LED")]),
                    pair("Wavelength start (nm)", number("wavelength_low", 250, step="any"),
                         "Wavelength end (nm)", number("wavelength_high", 800, step="any")),
                    pair("LED Savitzky–Golay window (points)", number("led_window", 12, min=1, step=1),
                         "LED polynomial order", number("led_order", 3, min=0, step=1)),
                    toggle("led_baseline", "Baseline-correct LED emission", True),
                    html.Label("Baseline exclusion (FWHM multiplier)"),
                    number("baseline_multiplier", 10, min=0, step="any"),
                    html.P("This LED-only smoothing is separate from Spectral Treatment.",
                           className="helper-text"),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("5 · Fit", className="step-label"),
                                  html.Span("Kinetic model")]),
                    html.Label("Method"),
                    dropdown("fit_method", (
                        ("concentrations", "Concentrations (independent NNLS)"),
                        ("regularized_concentrations", "Regularized concentrations"),
                        ("ode_absorbance", "Full-spectrum ODE absorbance"),
                        ("emission", "Emission (legacy)"),
                    ), "concentrations"),
                    pair("Initial Φ R→P", number("initial_rp", 0.5, min=0, step="any"),
                         "Initial Φ P→R", number("initial_pr", 0.5, min=0, step="any")),
                    pair("Lower Φ bound", number("yield_min", 0.0, min=0, step="any"),
                         "Upper Φ bound", number("yield_max", 1.0, min=0, step="any")),
                    html.Details(className="nested-tool", children=[
                        html.Summary("Method-specific controls"),
                        html.Label("Emission threshold fraction"),
                        number("emission_threshold", 0.01, min=0, max=1, step="any"),
                        html.Label("Concentration regularization strength"),
                        number("regularization_strength", 1.0, min=0, step="any"),
                        pair("Absorbance baseline order",
                             dropdown("absorbance_baseline_order",
                                      ((-1, "None"), (0, "Constant"), (1, "Linear")), 1),
                             "Robust-loss scale", number("robust_loss_scale", 0.02, min=0, step="any")),
                        html.Label("Expected reactant at PSS (%) · optional"),
                        number("expected_pss_reactant_percent", None, min=0, max=100,
                               step="any", placeholder="e.g. 23"),
                    ]),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("6 · Uncertainty", className="step-label"),
                                  html.Span("Molar absorptivity")]),
                    html.Label("ε propagation"),
                    dropdown("epsilon_method", (
                        ("none", "Off (default)"),
                        ("deterministic_extremes", "Deterministic wavelength bounds"),
                    ), "none"),
                    html.Label("Repeat-spectrum error metric"),
                    dropdown("epsilon_metric", (("sd", "Standard deviation"),
                                                 ("sem", "Standard error")), "sd"),
                    html.P("Deterministic bounds require AutoQY Spectral Treatment TSV files for both ε spectra.",
                           className="helper-text"),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("7 · Output", className="step-label"),
                                  html.Span("Files and display")]),
                    html.Label("Results directory"),
                    field("output_directory", type="text", value="results"),
                    html.Label("Residual color percentile"),
                    number("residual_percentile", 99, min=0, max=100, step="any"),
                    toggle("write_text", "Write TXT summary", True),
                    toggle("write_figures", "Write PNG and SVG", True),
                    toggle("write_json", "Write result JSON", True),
                    toggle("write_config", "Write configuration snapshot", True),
                    toggle("write_detailed", "Write detailed TSV data", True),
                    toggle("overwrite", "Allow replacing existing results", False),
                    toggle("relative_paths", "Save absolute inputs as paths relative to JSON", True),
                ]),
                html.Section(className="panel action-panel", children=[
                    html.P("8 · Analyze", className="step-label"),
                    html.Div(className="action-grid", children=[
                        html.Button("Save JSON", id="save-analysis-json", n_clicks=0,
                                    className="button button-secondary"),
                        html.Button("Validate", id="validate-analysis", n_clicks=0,
                                    className="button button-secondary"),
                    ]),
                    html.Button("Run analysis", id="run-analysis", n_clicks=0,
                                className="button button-primary run-button"),
                    html.Button("Compare fit methods", id="compare-fit-methods", n_clicks=0,
                                className="button button-secondary"),
                    dcc.Loading(html.Div(id="analysis-action-message", className="message"),
                                type="circle"),
                ]),
            ]),
            html.Div(className="plot-column analysis-results", children=[
                html.Section(className="result-strip analysis-result-strip", children=[
                    html.Div([html.Span("R → P"), html.Strong("—"), html.Small("Quantum yield")],
                             id="result-rp", className="result-card"),
                    html.Div([html.Span("P → R"), html.Strong("—"), html.Small("Quantum yield")],
                             id="result-pr", className="result-card result-card-accent"),
                    html.Div([html.Span("Fit"), html.Strong("—"), html.Small("Not run")],
                             id="result-fit", className="result-card result-card-neutral"),
                ]),
                html.Section(className="panel diagnostic-panel", children=[
                    html.H2("Preflight and fit health"),
                    html.Div("Validate or run to inspect the inputs.", id="preflight-checks",
                             className="diagnostic-list"),
                ]),
                html.Section(className="panel comparison-panel", children=[
                    html.H2("Fit-method comparison"),
                    html.Div("Use Compare fit methods to run nominal fits without ε propagation.",
                             id="method-comparison", className="helper-text"),
                ]),
                html.Section(className="plot-panel analysis-plot-panel", children=[
                    dcc.Tabs(id="analysis-plot-tabs", value="concentrations", children=[
                        dcc.Tab(label="Concentrations", value="concentrations",
                                children=[dcc.Graph(id="concentration-figure", figure=empty)]),
                        dcc.Tab(label="Fraction residual", value="fraction-residual",
                                children=[dcc.Graph(id="fraction-figure", figure=empty)]),
                        dcc.Tab(label="Spectra", value="spectra",
                                children=[dcc.Graph(id="spectra-figure", figure=empty)]),
                        dcc.Tab(label="Input compatibility", value="input-compatibility",
                                children=[dcc.Graph(id="input-diagnostic-figure", figure=empty)]),
                        dcc.Tab(label="Absorbance residuals", value="absorbance-residuals",
                                children=[dcc.Graph(id="residual-heatmap", figure=empty)]),
                    ]),
                ]),
                html.Section(className="panel", children=[
                    html.H2("Generated files"),
                    html.Div("Run an analysis to list its outputs.", id="analysis-output-files",
                             className="helper-text output-file-list"),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary("Live analysis.json preview"),
                    dcc.Textarea(id="analysis-json-preview", readOnly=True,
                                 className="json-preview"),
                ]),
                html.Details(open=False, className="panel tool-details error-panel", children=[
                    html.Summary("Python errors and warnings"),
                    html.Pre(id="analysis-python-error"),
                    html.Small("An empty panel means no Python error is active."),
                ]),
            ]),
        ]),
    ])

    @app.callback(
        Output({"type": "analysis-field", "name": ALL}, "value"),
        Output("load-analysis-message", "children"),
        Input("load-analysis-json", "n_clicks"),
        Input("choose-analysis-folder", "n_clicks"),
        Input({"type": "browse-analysis-file", "name": ALL}, "n_clicks"),
        State({"type": "analysis-field", "name": ALL}, "value"),
        State({"type": "analysis-field", "name": ALL}, "id"),
        prevent_initial_call=True,
    )
    def select_paths(_, __, ___, current_values, field_ids):
        values = dict(zip((item["name"] for item in field_ids), current_values))
        triggered = ctx.triggered_id
        message = ""
        try:
            if triggered == "load-analysis-json":
                selected = _choose_file(values.get("config_folder"), "JSON files|*.json")
                if selected:
                    config = load_config(selected)
                    values.update(_form_values(config.values))
                    values["config_folder"] = str(config.base_directory)
                    message = f"Loaded {selected}"
            elif triggered == "choose-analysis-folder":
                selected = _choose_folder(values.get("config_folder"))
                if selected:
                    values["config_folder"] = selected
                    message = f"JSON base folder: {selected}"
            elif isinstance(triggered, dict):
                name = triggered["name"]
                selected = _choose_file(values.get("config_folder"), _file_filter(name))
                if selected:
                    values[name] = selected
                    message = f"Selected {Path(selected).name}"
        except Exception as error:
            message = f"Could not select or load the file: {error}"
        return [values.get(item["name"]) for item in field_ids], message

    @app.callback(
        Output("analysis-json-preview", "value"),
        Input({"type": "analysis-field", "name": ALL}, "value"),
        State({"type": "analysis-field", "name": ALL}, "id"),
    )
    def preview_json(current_values, field_ids):
        try:
            values = dict(zip((item["name"] for item in field_ids), current_values))
            return json.dumps(_configuration(values), indent=2) + "\n"
        except Exception as error:
            return f"Cannot build JSON preview: {error}"

    @app.callback(
        Output("physical-parameter-preview", "children"),
        Output("physical-parameter-preview", "className"),
        Input({"type": "analysis-field", "name": ALL}, "value"),
        State({"type": "analysis-field", "name": ALL}, "id"),
    )
    def preview_physical_parameters(current_values, field_ids):
        values = dict(zip((item["name"] for item in field_ids), current_values))
        return _physical_parameter_preview(values)

    @app.callback(
        Output("tool-message", "children"),
        Input("launch-spectral-treatment", "n_clicks"),
        prevent_initial_call=True,
    )
    def launch_spectral_treatment(_):
        try:
            started = _start_spectral_treatment()
            address = "http://127.0.0.1:8051/"
            webbrowser.open(address)
            state = "started and is ready" if started else "was already running"
            return html.Span([
                f"Spectral Treatment {state}. ",
                html.A("Open it", href=address, target="_blank"),
                ".",
            ])
        except Exception as error:
            return f"Could not open Spectral Treatment: {error}"

    @app.callback(
        Output("analysis-action-message", "children"),
        Output("analysis-action-message", "className"),
        Output("result-rp", "children"),
        Output("result-pr", "children"),
        Output("result-fit", "children"),
        Output("concentration-figure", "figure"),
        Output("fraction-figure", "figure"),
        Output("spectra-figure", "figure"),
        Output("input-diagnostic-figure", "figure"),
        Output("residual-heatmap", "figure"),
        Output("preflight-checks", "children"),
        Output("method-comparison", "children"),
        Output("analysis-output-files", "children"),
        Output("analysis-python-error", "children"),
        Input("save-analysis-json", "n_clicks"),
        Input("validate-analysis", "n_clicks"),
        Input("run-analysis", "n_clicks"),
        Input("compare-fit-methods", "n_clicks"),
        State({"type": "analysis-field", "name": ALL}, "value"),
        State({"type": "analysis-field", "name": ALL}, "id"),
        prevent_initial_call=True,
    )
    def analyze(_, __, ___, ____, current_values, field_ids):
        values = dict(zip((item["name"] for item in field_ids), current_values))
        blank_card = [html.Span("R → P"), html.Strong("—"), html.Small("Quantum yield")]
        blank_back = [html.Span("P → R"), html.Strong("—"), html.Small("Quantum yield")]
        blank_fit = [html.Span("Fit"), html.Strong("—"), html.Small("Not run")]
        blank_figures = [_empty_figure(go, "Run a validated analysis to display this plot")] * 5
        blank_comparison = "Use Compare fit methods to run nominal fits without ε propagation."
        try:
            document = _configuration(values)
            base = Path(values.get("config_folder") or Path.cwd()).expanduser().resolve()
            config = AnalysisConfig(document, base)
            validate_config(config)
            with warnings.catch_warnings(record=True) as input_warnings:
                warnings.simplefilter("always")
                input_checks = _input_checks(config)
            input_checks.extend({
                "level": "warning",
                "title": "Input processing warning",
                "body": f"{item.category.__name__}: {item.message}",
            } for item in input_warnings)
            preflight = _render_checks(html, input_checks)
            action = ctx.triggered_id
            if action == "save-analysis-json":
                selected = _choose_save_json(base, "analysis.json")
                if not selected:
                    return ("Save cancelled.", "message", blank_card, blank_back, blank_fit,
                            *blank_figures, preflight, blank_comparison,
                            "No analysis was run.", "")
                target = Path(selected)
                saved = _portable_document(document, target.parent) if _on(values, "relative_paths") else document
                validate_config(AnalysisConfig(saved, target.parent))
                target.write_text(json.dumps(saved, indent=2) + "\n", encoding="utf-8")
                return (f"Saved and validated {target}", "message status-message status-ok",
                        blank_card, blank_back, blank_fit, *blank_figures,
                        preflight, blank_comparison, html.Code(str(target)), "")
            if action == "validate-analysis":
                has_stop = any(item["level"] == "stop" for item in input_checks)
                has_warning = any(item["level"] == "warning" for item in input_checks)
                if has_stop:
                    validation_message = "Configuration loaded, but preflight found blocking concerns."
                    validation_class = "message status-message status-stop"
                elif has_warning:
                    validation_message = "Configuration valid with preflight warnings."
                    validation_class = "message status-message status-warning"
                else:
                    validation_message = "Configuration valid. All referenced files and settings passed validation."
                    validation_class = "message status-message status-ok"
                warning_text = "\n".join(
                    f"{item.category.__name__}: {item.message}" for item in input_warnings
                )
                return (validation_message, validation_class, blank_card, blank_back, blank_fit,
                        *blank_figures, preflight, blank_comparison,
                        "No analysis was run.", warning_text)

            if action == "compare-fit-methods":
                with warnings.catch_warnings(record=True) as caught_warnings:
                    warnings.simplefilter("always")
                    comparison = _compare_fit_methods(config)
                warning_text = "\n".join(
                    f"{item.category.__name__}: {item.message}"
                    for item in caught_warnings
                )
                return (
                    "Nominal fit-method comparison completed.",
                    "message status-message status-ok",
                    blank_card, blank_back, blank_fit, *blank_figures,
                    preflight, _render_comparison(html, comparison, document["fit"]["method"]),
                    "No result files were generated by the comparison.", warning_text,
                )

            with warnings.catch_warnings(record=True) as caught_warnings:
                warnings.simplefilter("always")
                output = run_analysis(config)
            summary = result_summary(
                output.result, output.data,
                document["experiment"]["irradiation_wavelength_nm"],
            )
            figures = _interactive_figures(
                go, make_subplots, output.result, output.data,
                document["plots"]["absorbance_residual_percentile"],
            )
            rp = _yield_card(html, "R → P", summary, "R_to_P")
            pr = _yield_card(html, "P → R", summary, "P_to_R")
            fit = [html.Span("Fit"), html.Strong(_method_label(output.result.fit_method)),
                   html.Small(_fit_note(summary))]
            files = [html.Code(str(path)) for path in output.files]
            health_checks = _fit_health_checks(output.result, output.data, document)
            all_checks = input_checks + health_checks
            preflight = _render_checks(html, all_checks)
            warning_text = "\n".join(
                f"{item.category.__name__}: {item.message}"
                for item in [*input_warnings, *caught_warnings]
            )
            has_stop = any(item["level"] == "stop" for item in all_checks)
            has_warning = caught_warnings or any(
                item["level"] == "warning" for item in all_checks
            )
            if has_stop:
                message = "Analysis completed, but the fit is not reliable. Review fit health."
                message_class = "message status-message status-stop"
            elif has_warning:
                message = "Analysis completed with warnings. Review preflight and fit health."
                message_class = "message status-message status-warning"
            else:
                message = "Analysis completed successfully."
                message_class = "message status-message status-ok"
            return (message, message_class,
                    rp, pr, fit, *figures, preflight, blank_comparison,
                    files, warning_text)
        except Exception as error:
            message = f"{type(error).__name__}: {error}"
            return ("Validation or analysis stopped. See Python errors.",
                    "message status-message status-stop", blank_card, blank_back, blank_fit,
                    *blank_figures,
                    _render_checks(html, [{"level": "stop", "title": "Stopped", "body": message}]),
                    blank_comparison, "No files were generated.", message)

    return app


def _configuration(values):
    formats = {}
    for name, _ in (*SPECTRAL_INPUTS, TIMESTAMP_INPUT):
        kind = values.get(f"format_{name}")
        specification = {"type": kind}
        if kind == "generic_delimited":
            specification.update(
                delimiter=values.get(f"delimiter_{name}") or ",",
                header=_on(values, f"header_{name}"),
            )
        formats[name] = specification
    return {
        "schema_version": 1,
        "analysis": {
            "id": values.get("analysis_id") or "",
            "reactant_name": values.get("reactant_name") or "",
            "product_name": values.get("product_name") or "",
            "expected_pss_reactant_percent": values.get("expected_pss_reactant_percent"),
        },
        "inputs": {
            **{name: values.get(name) or "" for name, _ in (*SPECTRAL_INPUTS, TIMESTAMP_INPUT)},
            "formats": formats,
        },
        "experiment": {
            "volume_ul": _volume_ul(values),
            "power_mw": values.get("power_mw"),
            "power_error_mw": values.get("power_error_mw"),
            "thermal_back_reaction_s_1": values.get("thermal_rate"),
            "irradiation_wavelength_nm": values.get("irradiation_wavelength_nm"),
            "path_length_cm": values.get("path_length_cm"),
        },
        "processing": {
            "wavelength_range_nm": [values.get("wavelength_low"), values.get("wavelength_high")],
            "led_smoothing": {"method": "savgol", "window_points": values.get("led_window"),
                              "polynomial_order": values.get("led_order")},
            "led_baseline": {"enabled": _on(values, "led_baseline"),
                             "exclusion_fwhm_multiplier": values.get("baseline_multiplier")},
        },
        "fit": {
            "method": values.get("fit_method"),
            "emission_threshold_fraction": values.get("emission_threshold"),
            "regularization_strength": values.get("regularization_strength"),
            "absorbance_baseline_order": values.get("absorbance_baseline_order"),
            "robust_loss_scale": values.get("robust_loss_scale"),
            "initial_quantum_yields": {"R_to_P": values.get("initial_rp"),
                                        "P_to_R": values.get("initial_pr")},
            "quantum_yield_bounds": {"minimum": values.get("yield_min"),
                                      "maximum": values.get("yield_max")},
        },
        "plots": {"absorbance_residual_percentile": values.get("residual_percentile")},
        "uncertainty": {"epsilon": {"method": values.get("epsilon_method"),
                                      "error_metric": values.get("epsilon_metric")}},
        "outputs": {
            "directory": values.get("output_directory") or "results",
            "stem": values.get("output_stem") or "",
            "write_text": _on(values, "write_text"),
            "write_figures": _on(values, "write_figures"),
            "write_json": _on(values, "write_json"),
            "write_config": _on(values, "write_config"),
            "write_detailed_data": _on(values, "write_detailed"),
            "overwrite": _on(values, "overwrite"),
        },
    }


def _form_values(document):
    inputs, experiment = document["inputs"], document["experiment"]
    processing, fit = document["processing"], document["fit"]
    output = document["outputs"]
    values = {
        "analysis_id": document["analysis"]["id"],
        "reactant_name": document["analysis"]["reactant_name"],
        "product_name": document["analysis"]["product_name"],
        "output_stem": output["stem"],
        "volume_value": experiment["volume_ul"],
        "volume_unit": "ul",
        "power_mw": experiment["power_mw"],
        "power_error_mw": experiment["power_error_mw"],
        "thermal_rate": experiment["thermal_back_reaction_s_1"],
        "irradiation_wavelength_nm": experiment["irradiation_wavelength_nm"],
        "path_length_cm": experiment["path_length_cm"],
        "wavelength_low": processing["wavelength_range_nm"][0],
        "wavelength_high": processing["wavelength_range_nm"][1],
        "led_window": processing["led_smoothing"]["window_points"],
        "led_order": processing["led_smoothing"]["polynomial_order"],
        "led_baseline": _toggle_value(processing["led_baseline"]["enabled"]),
        "baseline_multiplier": processing["led_baseline"]["exclusion_fwhm_multiplier"],
        "fit_method": fit["method"],
        "emission_threshold": fit.get("emission_threshold_fraction", 0.01),
        "regularization_strength": fit.get("regularization_strength", 1),
        "absorbance_baseline_order": fit.get("absorbance_baseline_order", 1),
        "robust_loss_scale": fit.get("robust_loss_scale", 0.02),
        "expected_pss_reactant_percent": document["analysis"].get(
            "expected_pss_reactant_percent"
        ),
        "initial_rp": fit["initial_quantum_yields"]["R_to_P"],
        "initial_pr": fit["initial_quantum_yields"]["P_to_R"],
        "yield_min": fit["quantum_yield_bounds"]["minimum"],
        "yield_max": fit["quantum_yield_bounds"]["maximum"],
        "residual_percentile": document["plots"]["absorbance_residual_percentile"],
        "epsilon_method": document.get("uncertainty", {}).get("epsilon", {}).get("method", "none"),
        "epsilon_metric": document.get("uncertainty", {}).get("epsilon", {}).get("error_metric", "sd"),
        "output_directory": output["directory"],
        "write_text": _toggle_value(output["write_text"]),
        "write_figures": _toggle_value(output["write_figures"]),
        "write_json": _toggle_value(output["write_json"]),
        "write_config": _toggle_value(output["write_config"]),
        "write_detailed": _toggle_value(output.get("write_detailed_data", False)),
        "overwrite": _toggle_value(output["overwrite"]),
        "relative_paths": ["on"],
    }
    formats = inputs.get("formats", {})
    for name, _ in (*SPECTRAL_INPUTS, TIMESTAMP_INPUT):
        values[name] = inputs[name]
        spec = formats.get(name, "ahk_csv" if name == "timestamps" else "spectragryph_tsv")
        spec = {"type": spec} if isinstance(spec, str) else spec
        values[f"format_{name}"] = spec.get("type")
        values[f"delimiter_{name}"] = spec.get("delimiter", ",")
        values[f"header_{name}"] = _toggle_value(spec.get("header", True))
    return values


def _portable_document(document, folder):
    copied = json.loads(json.dumps(document))
    for name, _ in (*SPECTRAL_INPUTS, TIMESTAMP_INPUT):
        path = Path(copied["inputs"][name]).expanduser()
        if path.is_absolute():
            try:
                copied["inputs"][name] = str(path.relative_to(folder))
            except ValueError:
                try:
                    copied["inputs"][name] = str(Path(_relative_path(path, folder)))
                except ValueError:
                    copied["inputs"][name] = str(path)
    return copied


def _relative_path(path, folder):
    import os
    return os.path.relpath(path, folder)


def _volume_ul(values):
    value = values.get("volume_value")
    if value is None:
        return None
    return float(value) * (1000.0 if values.get("volume_unit") == "ml" else 1.0)


def _physical_parameter_preview(values):
    try:
        volume_ul = _volume_ul(values)
        if volume_ul is None:
            return "Enter the sample volume.", "message status-message status-warning"
        text = f"Interpreted volume: {volume_ul:g} µL = {volume_ul / 1000:g} mL."
        rate = values.get("thermal_rate")
        if rate is not None and float(rate) > 0:
            half_life = np.log(2) / float(rate)
            text += f" Thermal half-life: {half_life:g} s."
        warning = volume_ul < 50
        if warning:
            text += " Confirm that the volume was not entered in mL."
        return text, ("message status-message status-warning" if warning
                      else "message status-message status-ok")
    except (TypeError, ValueError, ZeroDivisionError) as error:
        return f"Cannot interpret physical parameters: {error}", "message status-message status-stop"


def _input_checks(config):
    values = config.values
    experiment, processing = values["experiment"], values["processing"]
    checks = []
    volume_ul = float(experiment["volume_ul"])
    checks.append(_check(
        "warning" if volume_ul < 50 else "ok",
        "Sample volume",
        f"{volume_ul:g} µL = {volume_ul / 1000:g} mL" +
        ("; confirm that an mL value was not entered as µL" if volume_ul < 50 else ""),
    ))

    wavelengths, absorbance = load_spectra(
        config.input_path("measurement_spectra"),
        input_format(config, "measurement_spectra"),
    )
    timestamps = load_timestamps(
        config.input_path("timestamps"), input_format(config, "timestamps")
    )
    count_ok = absorbance.shape[1] == len(timestamps)
    checks.append(_check(
        "ok" if count_ok else "stop", "Spectra and timestamps",
        f"{absorbance.shape[1]} spectra and {len(timestamps)} timestamps" +
        (" are aligned" if count_ok else " do not match"),
    ))
    if len(timestamps) > 1:
        intervals = np.diff(timestamps)
        median_interval = float(np.median(intervals))
        duration = float(timestamps[-1] - timestamps[0])
        checks.append(_check(
            "ok", "Acquisition timeline",
            f"median interval {median_interval:.4g} s; duration {duration:.4g} s",
        ))
    else:
        median_interval, duration = np.nan, 0.0

    rate = float(experiment["thermal_back_reaction_s_1"])
    if rate > 0:
        half_life = float(np.log(2) / rate)
        suspicious = ((np.isfinite(median_interval) and half_life < median_interval) or
                      (duration > 0 and half_life < duration / 10))
        checks.append(_check(
            "warning" if suspicious else "ok", "Thermal back-reaction",
            f"k = {rate:g} s⁻¹; half-life {half_life:.4g} s" +
            (" is short relative to the experiment" if suspicious else ""),
        ))
    else:
        checks.append(_check("ok", "Thermal back-reaction", "disabled (k = 0 s⁻¹)"))

    epsilon_r = _load_reference(config, "reactant_absorptivity")
    epsilon_p = _load_reference(config, "product_absorptivity")
    led_wavelengths, _ = load_spectrum(
        config.input_path("led_emission"), input_format(config, "led_emission")
    )
    requested_low, requested_high = processing["wavelength_range_nm"]
    common_low = max(float(wavelengths[0]), float(epsilon_r[0][0]),
                     float(epsilon_p[0][0]), float(requested_low))
    common_high = min(float(wavelengths[-1]), float(epsilon_r[0][-1]),
                      float(epsilon_p[0][-1]), float(requested_high))
    checks.append(_check(
        "ok" if common_high > common_low else "stop", "Wavelength overlap",
        f"usable measurement/reference range {common_low:.1f}–{common_high:.1f} nm",
    ))
    irradiation = float(experiment["irradiation_wavelength_nm"])
    led_ok = float(led_wavelengths[0]) <= irradiation <= float(led_wavelengths[-1])
    checks.append(_check(
        "ok" if led_ok else "stop", "LED coverage",
        f"irradiation wavelength {irradiation:g} nm; LED file spans "
        f"{led_wavelengths[0]:.1f}–{led_wavelengths[-1]:.1f} nm",
    ))
    if common_high > common_low:
        grid = wavelengths[(wavelengths >= common_low) & (wavelengths <= common_high)]
        reference_matrix = np.column_stack((
            np.interp(grid, epsilon_r[0], epsilon_r[1]),
            np.interp(grid, epsilon_p[0], epsilon_p[1]),
        ))
        condition = float(np.linalg.cond(reference_matrix))
        level = "stop" if condition > 1e6 else "warning" if condition > 1e3 else "ok"
        checks.append(_check(
            level, "Reference-spectrum separation",
            f"two-spectrum condition number {condition:.3g}; lower is better",
        ))
    return checks


def _load_reference(config, name):
    path = config.input_path(name)
    nominal = load_epsilon_nominal(path)
    if nominal is not None:
        return nominal
    return load_spectrum(path, input_format(config, name))


def _fit_health_checks(result, data, document):
    checks = []
    measured = result.concentration_fit.concentrations
    fitted = result.yield_fit.concentrations
    totals = measured.sum(axis=1)
    initial_total = float(totals[0])
    initial_product = (float(measured[0, 1]) / initial_total if initial_total > 0 else 0.0)
    checks.append(_check(
        "warning" if initial_product > 0.02 else "ok", "Initial product",
        f"{initial_product:.2%} of total concentration" +
        ("; verify references and baseline, while retaining it if physically real"
         if initial_product > 0.02 else ""),
    ))

    measured_fraction = result.concentration_fit.fractions[:, 0]
    fitted_total = fitted.sum(axis=1)
    fitted_fraction = np.divide(fitted[:, 0], fitted_total, out=np.zeros_like(fitted_total),
                                where=fitted_total != 0)
    fraction_residual = measured_fraction - fitted_fraction
    fraction_rmse = float(np.sqrt(np.mean(fraction_residual ** 2)))
    fraction_max = float(np.max(np.abs(fraction_residual)))
    fraction_level = "stop" if fraction_rmse > 0.05 else "warning" if fraction_rmse > 0.02 else "ok"
    checks.append(_check(
        fraction_level, "Kinetic fraction residual",
        f"RMSE {fraction_rmse:.4g}; maximum |data − fit| {fraction_max:.4g}",
    ))

    epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
    fitted_absorbance = fitted @ epsilon * data.path_length_cm
    if result.yield_fit.absorbance_correction is not None:
        fitted_absorbance += result.yield_fit.absorbance_correction
    absorbance_residual = result.absorbance.T - fitted_absorbance
    absorbance_rmse = float(np.sqrt(np.mean(absorbance_residual ** 2)))
    absorbance_scale = max(float(np.max(np.abs(result.absorbance))), np.finfo(float).eps)
    relative_rmse = absorbance_rmse / absorbance_scale
    absorbance_level = "stop" if relative_rmse > 0.05 else "warning" if relative_rmse > 0.02 else "ok"
    checks.append(_check(
        absorbance_level, "Full-spectrum reconstruction",
        f"absorbance RMSE {absorbance_rmse:.4g} ({relative_rmse:.2%} of maximum absorbance)",
    ))

    total_cv = float(np.std(totals) / np.mean(totals)) if np.mean(totals) else np.inf
    total_level = "stop" if total_cv > 0.05 else "warning" if total_cv > 0.02 else "ok"
    checks.append(_check(
        total_level, "Concentration conservation",
        f"total-concentration coefficient of variation {total_cv:.3%}",
    ))

    fit = result.yield_fit
    condition = float(fit.jacobian_condition)
    optimizer_level = ("stop" if not fit.optimizer_success or condition > 1e8 else
                       "warning" if condition > 1e6 or any(fit.active_bounds) else "ok")
    bound_text = "; a quantum yield touches its bound" if any(fit.active_bounds) else ""
    checks.append(_check(
        optimizer_level, "Optimizer identifiability",
        f"Jacobian condition number {condition:.3g}{bound_text}",
    ))

    model_time = _transition_time(np.asarray(data.timestamps), 1 - fitted_fraction)
    data_time = _transition_time(np.asarray(data.timestamps),
                                 1 - measured_fraction)
    if model_time is not None and data_time is not None:
        interval = float(np.median(np.diff(data.timestamps))) if len(data.timestamps) > 1 else 0.0
        mismatch = abs(model_time - data_time) > max(2 * interval, 0.25 * max(data_time, interval))
        checks.append(_check(
            "warning" if mismatch else "ok", "Conversion timescale",
            f"95% transition: model {model_time:.4g} s; data {data_time:.4g} s",
        ))

    expected = document["analysis"].get("expected_pss_reactant_percent")
    fitted_pss = float(result.extrapolated_pss[0] / np.sum(result.extrapolated_pss) * 100)
    if expected is None:
        checks.append(_check(
            "ok", "Extrapolated PSS", f"reactant {fitted_pss:.2f}% (no expected value entered)"
        ))
    else:
        difference = abs(fitted_pss - float(expected))
        checks.append(_check(
            "warning" if difference > 5 else "ok", "Expected versus fitted PSS",
            f"reactant expected {float(expected):.2f}%; fitted {fitted_pss:.2f}%",
        ))

    errors = np.asarray(result.yield_errors, float)
    if np.any(~np.isfinite(errors)) or np.any(errors >= 1):
        checks.append(_check(
            "stop", "Quantum-yield uncertainty",
            "uncertainty is at least 100 percentage points; treat the fit as unidentifiable",
        ))
    uncertainty = result.epsilon_uncertainty
    if uncertainty is not None and len(uncertainty.bound_labels):
        worst = int(np.argmax(np.max(uncertainty.bound_optimizer_errors, axis=1)))
        worst_errors = uncertainty.bound_optimizer_errors[worst] * 100
        level = "stop" if np.max(worst_errors) >= 100 else "warning" if np.max(worst_errors) >= 25 else "ok"
        checks.append(_check(
            level, "Most uncertain ε-bound combination",
            f"{uncertainty.bound_labels[worst]}: optimizer errors "
            f"{worst_errors[0]:.3g}% and {worst_errors[1]:.3g}%",
        ))
    return checks


def _transition_time(times, values, fraction=0.95):
    values = np.asarray(values, float)
    change = float(values[-1] - values[0])
    if not np.isfinite(change) or abs(change) <= np.finfo(float).eps:
        return None
    target = values[0] + fraction * change
    indices = np.flatnonzero(values >= target) if change > 0 else np.flatnonzero(values <= target)
    return float(times[indices[0]]) if len(indices) else None


def _compare_fit_methods(config):
    methods = ("concentrations", "regularized_concentrations", "ode_absorbance")
    rows = []
    with TemporaryDirectory(prefix="autoqy-method-comparison-") as temporary:
        for method in methods:
            values = deepcopy(config.values)
            values["fit"]["method"] = method
            values.setdefault("uncertainty", {}).setdefault("epsilon", {})["method"] = "none"
            for name in ("write_text", "write_figures", "write_json", "write_config",
                         "write_detailed_data"):
                values["outputs"][name] = False
            output = run_analysis(
                AnalysisConfig(values, config.base_directory, config.source), temporary
            )
            result, data = output.result, output.data
            measured_fraction = result.concentration_fit.fractions[:, 0]
            fitted_total = result.yield_fit.concentrations.sum(axis=1)
            fitted_fraction = np.divide(
                result.yield_fit.concentrations[:, 0], fitted_total,
                out=np.zeros_like(fitted_total), where=fitted_total != 0,
            )
            epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
            fitted_absorbance = result.yield_fit.concentrations @ epsilon * data.path_length_cm
            if result.yield_fit.absorbance_correction is not None:
                fitted_absorbance += result.yield_fit.absorbance_correction
            rows.append({
                "method": method,
                "values": result.yield_fit.values * 100,
                "errors": result.yield_errors * 100,
                "fraction_rmse": float(np.sqrt(np.mean(
                    (measured_fraction - fitted_fraction) ** 2
                ))),
                "absorbance_rmse": float(np.sqrt(np.mean(
                    (result.absorbance.T - fitted_absorbance) ** 2
                ))),
                "active_bounds": any(result.yield_fit.active_bounds),
            })
    return rows


def _render_comparison(html, rows, selected_method):
    headings = ("Method", "Φ R→P", "Φ P→R", "Fraction RMSE", "Absorbance RMSE", "Health")
    body = []
    for row in rows:
        health = "Bound hit" if row["active_bounds"] else "Nominal"
        body.append(html.Tr(className="comparison-selected" if row["method"] == selected_method else "",
                            children=[
            html.Td(_method_label(row["method"])),
            html.Td(f"{row['values'][0]:.3g} ± {row['errors'][0]:.2g}%"),
            html.Td(f"{row['values'][1]:.3g} ± {row['errors'][1]:.2g}%"),
            html.Td(f"{row['fraction_rmse']:.4g}"),
            html.Td(f"{row['absorbance_rmse']:.4g}"),
            html.Td(health),
        ]))
    return html.Div([
        html.P("Nominal ε comparison; no result files were written.", className="helper-text"),
        html.Div(className="comparison-table-wrap", children=[
            html.Table(className="comparison-table", children=[
                html.Thead(html.Tr([html.Th(name) for name in headings])),
                html.Tbody(body),
            ])
        ]),
    ])


def _yield_card(html, label, summary, name):
    value = float(summary["quantum_yield_percent"][name])
    error = float(summary["quantum_yield_error_percent"][name])
    if not np.isfinite(error) or error >= 100:
        return [html.Span(label), html.Strong(f"Nominal {value:.3g}%"),
                html.Small("Unidentifiable uncertainty · inspect fit health")]
    formatted = summary["quantum_yield_formatted_percent"][name]
    return [html.Span(label),
            html.Strong(f"{formatted['value']} ± {formatted['error']}%"),
            html.Small("Quantum yield")]


def _check(level, title, body):
    return {"level": level, "title": title, "body": body}


def _render_checks(html, checks):
    return [html.Div(className=f"diagnostic-item diagnostic-{item['level']}", children=[
        html.Strong(item["title"]), html.Span(item["body"]),
    ]) for item in checks]


def _interactive_figures(go, make_subplots, result, data, residual_percentile):
    times = np.asarray(data.timestamps)
    measured = result.concentration_fit.concentrations
    fitted = result.yield_fit.concentrations
    measured_fraction = result.concentration_fit.fractions[:, 0]
    fitted_total = fitted.sum(axis=1)
    fitted_fraction = np.divide(fitted[:, 0], fitted_total, out=np.zeros_like(fitted_total),
                                where=fitted_total != 0)
    fraction_residual = measured_fraction - fitted_fraction
    uncertainty = result.epsilon_uncertainty
    blue, orange = "#2d6f8e", "#d67b36"

    concentration = go.Figure()
    if uncertainty is not None:
        for index, (name, colour) in enumerate((("Reactant", blue), ("Product", orange))):
            _add_band(go, concentration, times,
                      uncertainty.concentration_fit_minimum[:, index],
                      uncertainty.concentration_fit_maximum[:, index], colour,
                      f"{name} ε-bound fit range")
    for index, (name, colour) in enumerate((("Reactant", blue), ("Product", orange))):
        error_y = None
        if uncertainty is not None:
            error_y = {
                "type": "data", "symmetric": False, "visible": True,
                "array": uncertainty.concentration_data_maximum[:, index] - measured[:, index],
                "arrayminus": measured[:, index] - uncertainty.concentration_data_minimum[:, index],
                "color": _hex_rgba(colour, 0.38), "thickness": 1, "width": 2,
            }
        concentration.add_trace(go.Scatter(
            x=times, y=measured[:, index], mode="markers", name=f"{name} data",
            marker={"symbol": "circle-open", "size": 7, "color": colour},
            error_y=error_y,
        ))
        concentration.add_trace(go.Scatter(
            x=times, y=fitted[:, index], mode="lines", name=f"{name} fit",
            line={"color": colour, "width": 2},
        ))
    _style_figure(concentration, "Concentration fit", "Irradiation time (s)",
                  "Concentration (mol/L)")

    fraction = go.Figure()
    if uncertainty is not None:
        _add_band(go, fraction, times, uncertainty.fraction_residual_minimum,
                  uncertainty.fraction_residual_maximum, blue, "ε-bound range")
    fraction.add_trace(go.Scatter(x=times, y=fraction_residual, mode="lines+markers",
                                  name="Nominal residual", line={"color": blue}))
    fraction.add_hline(y=0, line={"color": "#20242c", "width": 1})
    _style_figure(fraction, "Reactant fraction residual", "Irradiation time (s)",
                  "Fraction data − fit")

    spectra = go.Figure()
    measured_absorbance = result.absorbance.T
    count = len(times)
    indices = np.unique(np.linspace(0, count - 1, min(count, 80), dtype=int))
    for index in indices:
        ratio = 0 if count == 1 else index / (count - 1)
        colour = _interpolate_colour((45, 111, 142), (214, 123, 54), ratio)
        spectra.add_trace(go.Scatter(
            x=result.wavelengths, y=measured_absorbance[index], mode="lines",
            name=f"{times[index]:g} s", showlegend=False,
            line={"color": colour, "width": 1},
            hovertemplate=f"time={times[index]:g} s<br>λ=%{{x:.1f}} nm<br>A=%{{y:.5g}}<extra></extra>",
        ))
    _style_figure(spectra, "Measured absorption spectra over time", "Wavelength (nm)",
                  "Absorbance")
    spectra.update_layout(coloraxis={"colorscale": [[0, blue], [1, orange]],
                                     "cmin": float(times.min()), "cmax": float(times.max()),
                                     "colorbar": {"title": "Time (s)"}})
    spectra.add_trace(go.Scatter(x=[None], y=[None], mode="markers", showlegend=False,
                                 marker={"color": [float(times.min())], "coloraxis": "coloraxis"},
                                 hoverinfo="skip"))

    epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
    fitted_absorbance = fitted @ epsilon * data.path_length_cm
    if result.yield_fit.absorbance_correction is not None:
        fitted_absorbance += result.yield_fit.absorbance_correction

    input_diagnostic = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.12,
        specs=[[{"secondary_y": True}], [{}]],
        row_heights=[0.46, 0.54],
        subplot_titles=("Reference absorptivity and fitted LED", "Measured and reconstructed spectra"),
    )
    input_diagnostic.add_trace(go.Scatter(
        x=result.wavelengths, y=result.epsilon_r, mode="lines",
        name="Reactant ε", line={"color": blue, "width": 2},
    ), row=1, col=1, secondary_y=False)
    input_diagnostic.add_trace(go.Scatter(
        x=result.wavelengths, y=result.epsilon_p, mode="lines",
        name="Product ε", line={"color": orange, "width": 2},
    ), row=1, col=1, secondary_y=False)
    led_processed = process_led(
        *data.led, data.baseline_correct_led, data.led_smoothing_window,
        data.led_polynomial_order, data.baseline_exclusion_fwhm_multiplier,
    )
    led_aligned = np.interp(result.wavelengths, data.led[0], led_processed)
    led_scale = max(float(np.max(led_aligned)), np.finfo(float).eps)
    input_diagnostic.add_trace(go.Scatter(
        x=result.wavelengths, y=led_aligned / led_scale, mode="lines",
        name="Processed LED (normalized)", line={"color": "#8a5a9b", "width": 1.8},
    ), row=1, col=1, secondary_y=True)
    for index, label, colour in ((0, "Initial", blue), (-1, "Final", orange)):
        input_diagnostic.add_trace(go.Scatter(
            x=result.wavelengths, y=measured_absorbance[index], mode="lines",
            name=f"{label} measured", line={"color": colour, "width": 2},
        ), row=2, col=1)
        input_diagnostic.add_trace(go.Scatter(
            x=result.wavelengths, y=fitted_absorbance[index], mode="lines",
            name=f"{label} reconstructed",
            line={"color": colour, "width": 1.8, "dash": "dash"},
        ), row=2, col=1)
    _style_figure(input_diagnostic, "Input and reconstruction compatibility", "", "")
    input_diagnostic.update_yaxes(title_text="ε (M⁻¹ cm⁻¹)", row=1, col=1,
                                  secondary_y=False)
    input_diagnostic.update_yaxes(title_text="LED (normalized)", row=1, col=1,
                                  secondary_y=True)
    input_diagnostic.update_yaxes(title_text="Absorbance", row=2, col=1)
    input_diagnostic.update_xaxes(title_text="Wavelength (nm)", row=2, col=1)

    absorbance_residual = measured_absorbance - fitted_absorbance
    limit = float(np.percentile(np.abs(absorbance_residual), residual_percentile))
    limit = limit or float(np.finfo(float).eps)
    heatmap = go.Figure(go.Heatmap(
        x=result.wavelengths, y=times, z=absorbance_residual,
        colorscale="RdBu", reversescale=True, zmin=-limit, zmax=limit,
        colorbar={"title": "Data − fit"},
        hovertemplate="time=%{y:.3g} s<br>λ=%{x:.1f} nm<br>residual=%{z:.5g}<extra></extra>",
    ))
    title = "Absorbance residuals"
    if uncertainty is not None:
        title += (f" · ε-bound RMSE {uncertainty.absorbance_residual_rmse_minimum:.3g}–"
                  f"{uncertainty.absorbance_residual_rmse_maximum:.3g}")
    _style_figure(heatmap, title, "Wavelength (nm)", "Irradiation time (s)")
    heatmap.update_yaxes(autorange="reversed")
    return concentration, fraction, spectra, input_diagnostic, heatmap


def _add_band(go, figure, x, lower, upper, colour, name):
    rgba = _hex_rgba(colour, 0.15)
    figure.add_trace(go.Scatter(x=x, y=lower, mode="lines", line={"width": 0},
                                hoverinfo="skip", showlegend=False, legendgroup=name))
    figure.add_trace(go.Scatter(x=x, y=upper, mode="lines", line={"width": 0},
                                fill="tonexty", fillcolor=rgba, name=name,
                                legendgroup=name))


def _style_figure(figure, title, x_title, y_title):
    figure.update_layout(
        template="plotly_white", height=740,
        title={"text": title, "x": 0.02, "y": 0.98,
               "yanchor": "top", "yref": "container"},
        margin={"l": 72, "r": 34, "t": 205, "b": 60}, hovermode="closest",
        legend={"orientation": "h", "y": 1.06, "yanchor": "bottom", "x": 0,
                "xanchor": "left", "bgcolor": "rgba(255,255,255,.9)"},
    )
    figure.update_xaxes(title=x_title, showgrid=False)
    figure.update_yaxes(title=y_title, gridcolor="rgba(32,36,44,.10)")
    return figure


def _empty_figure(go, message):
    figure = go.Figure()
    figure.add_annotation(text=message, showarrow=False, font={"color": "#6c7280"})
    return _style_figure(figure, "Analysis result", "", "")


def _method_label(method):
    return {
        "concentrations": "Concentrations",
        "regularized_concentrations": "Regularized",
        "ode_absorbance": "ODE absorbance",
        "emission": "Emission",
    }[method]


def _fit_note(summary):
    uncertainty = summary.get("epsilon_uncertainty")
    return (f"{uncertainty['bound_combination_count']} ε-bound combinations"
            if uncertainty else "Power + optimizer uncertainty")


def _on(values, name):
    return "on" in (values.get(name) or [])


def _toggle_value(enabled):
    return ["on"] if enabled else []


def _interpolate_colour(start, stop, ratio):
    rgb = tuple(round(a + (b - a) * ratio) for a, b in zip(start, stop))
    return f"rgb{rgb}"


def _hex_rgba(colour, alpha):
    colour = colour.lstrip("#")
    rgb = tuple(int(colour[index:index + 2], 16) for index in (0, 2, 4))
    return f"rgba({rgb[0]},{rgb[1]},{rgb[2]},{alpha})"


def _file_filter(name):
    if name == "timestamps":
        return "Timestamp files|*.csv;*.txt;*.tsv|All files|*.*"
    return "Spectral files|*.dat;*.txt;*.tsv;*.csv|All files|*.*"


def _choose_file(initial_directory=None, file_filter="All files|*.*"):
    initial = _powershell_quote(str(initial_directory or ""))
    file_filter = _powershell_quote(file_filter)
    script = (
        "Add-Type -AssemblyName System.Windows.Forms; "
        + _foreground_owner_script() +
        "$dialog = New-Object System.Windows.Forms.OpenFileDialog; "
        f"$dialog.Filter = '{file_filter}'; "
        f"if ('{initial}' -and (Test-Path -LiteralPath '{initial}')) {{ $dialog.InitialDirectory = '{initial}' }}; "
        "$result = $dialog.ShowDialog($owner); "
        "if ($result -eq [System.Windows.Forms.DialogResult]::OK) { Write-Output $dialog.FileName }; "
        "$owner.Close(); $owner.Dispose()"
    )
    selected = _run_powershell_dialog(script)
    return selected[0] if selected else None


def _choose_folder(initial_directory=None):
    initial = _powershell_quote(str(initial_directory or ""))
    script = (
        "Add-Type -AssemblyName System.Windows.Forms; "
        + _foreground_owner_script() +
        "$dialog = New-Object System.Windows.Forms.FolderBrowserDialog; "
        f"if ('{initial}' -and (Test-Path -LiteralPath '{initial}')) {{ $dialog.SelectedPath = '{initial}' }}; "
        "$result = $dialog.ShowDialog($owner); "
        "if ($result -eq [System.Windows.Forms.DialogResult]::OK) { Write-Output $dialog.SelectedPath }; "
        "$owner.Close(); $owner.Dispose()"
    )
    selected = _run_powershell_dialog(script)
    return selected[0] if selected else None


def _choose_save_json(initial_directory, filename):
    initial = _powershell_quote(str(initial_directory))
    filename = _powershell_quote(filename)
    script = (
        "Add-Type -AssemblyName System.Windows.Forms; "
        + _foreground_owner_script() +
        "$dialog = New-Object System.Windows.Forms.SaveFileDialog; "
        "$dialog.Filter = 'JSON files|*.json'; $dialog.DefaultExt = 'json'; "
        "$dialog.AddExtension = $true; $dialog.OverwritePrompt = $true; "
        f"$dialog.InitialDirectory = '{initial}'; $dialog.FileName = '{filename}'; "
        "$result = $dialog.ShowDialog($owner); "
        "if ($result -eq [System.Windows.Forms.DialogResult]::OK) { Write-Output $dialog.FileName }; "
        "$owner.Close(); $owner.Dispose()"
    )
    selected = _run_powershell_dialog(script)
    return selected[0] if selected else None


def _run_powershell_dialog(script):
    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    completed = subprocess.run(
        ["powershell.exe", "-NoProfile", "-STA", "-Command", script],
        check=True, capture_output=True, text=True, creationflags=flags,
    )
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def _foreground_owner_script():
    return (
        "$owner = New-Object System.Windows.Forms.Form; "
        "$owner.TopMost = $true; $owner.ShowInTaskbar = $false; "
        "$owner.StartPosition = 'CenterScreen'; $owner.Width = 1; $owner.Height = 1; "
        "$owner.Opacity = 0; $owner.Show(); $owner.Activate(); "
    )


def _powershell_quote(value):
    return str(value).replace("'", "''")


def _port_open(host, port):
    with socket.socket() as connection:
        connection.settimeout(0.2)
        return connection.connect_ex((host, port)) == 0


def _start_spectral_treatment(host="127.0.0.1", port=8051, timeout=15):
    """Start Spectral Treatment and return only after its local port is ready."""
    if _port_open(host, port):
        return False
    inherited_paths = [path for path in sys.path if path]
    code = (
        "import sys; "
        f"sys.path[:0] = {inherited_paths!r}; "
        "from autoqy_core.tools.smoother_gui import run_server; "
        f"run_server({host!r}, {int(port)}, False)"
    )
    process = subprocess.Popen(
        [sys.executable, "-c", code], cwd=str(Path(__file__).parents[2]),
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    deadline = time.monotonic() + float(timeout)
    while time.monotonic() < deadline:
        if _port_open(host, port):
            return True
        if process.poll() is not None:
            raise RuntimeError(
                f"the Spectral Treatment process exited with code {process.returncode}"
            )
        time.sleep(0.2)
    process.terminate()
    raise TimeoutError(f"Spectral Treatment did not become ready on port {port}")


def run_server(host="127.0.0.1", port=8052, open_browser=True):
    from werkzeug.serving import make_server

    app = create_app()
    server = make_server(host, port, app.server, threaded=True)
    window_state, window_lock = app.server.config["AUTOQY_WINDOW_STATE"]

    def close_after_browser():
        while True:
            time.sleep(0.5)
            with window_lock:
                requested = window_state["close_requested_at"]
            if requested is not None and time.monotonic() - requested >= 3:
                server.shutdown()
                return

    Thread(target=close_after_browser, daemon=True).start()
    if open_browser:
        Timer(1, lambda: webbrowser.open(f"http://{host}:{port}")).start()
    server.serve_forever()


def main(argv=None):
    parser = ArgumentParser(prog="autoqy-analysis-gui")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8052, type=int)
    parser.add_argument("--no-open", action="store_true")
    args = parser.parse_args(argv)
    run_server(args.host, args.port, not args.no_open)


if __name__ == "__main__":
    main()
