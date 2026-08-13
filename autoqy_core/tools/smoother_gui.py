"""Browser GUI for treating spectra and calculating molar absorptivity."""

from argparse import ArgumentParser
import base64
from pathlib import Path
from threading import Timer
import webbrowser

import numpy as np

from ..epsilon import (EpsilonResult, NMRSubtractionResult,
                       calculate_epsilon_statistics, export_epsilon_tsv,
                       export_nmr_subtraction_tsv, load_epsilon_tsv,
                       nonnegative_error_bounds, reconstruct_product_from_nmr)
from ..smoother import (SpectralDataset, analyze_svd, baseline_spectra,
                        load_spectral_bytes, savgol_window_points, select_wavelengths,
                        smooth_reconstruction)


def create_app():
    try:
        from dash import ALL, Dash, Input, Output, State, ctx, dcc, html, no_update
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as error:
        raise RuntimeError("The spectral treatment GUI requires the 'power-gui' optional dependencies") from error

    assets = Path(__file__).parents[1] / "assets"
    app = Dash(__name__, assets_folder=str(assets), suppress_callback_exceptions=True)
    app.title = "AutoQY Spectral treatment"
    app.layout = html.Div(className="app-shell", children=[
        html.Header(className="app-header", children=[html.Div([
            html.P("AUTOQY CORE", className="eyebrow"),
            html.H1("Spectral treatment"),
            html.P("Prepare wavelength-resolved data and calculate molar absorptivity.",
                   className="subtitle"),
        ]), html.Div("Local session", className="local-badge")]),
        html.Main(className="workspace epsilon-workspace", children=[
            html.Aside(className="control-column smoother-controls", children=[
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("1 · Data", className="step-label"),
                                  "Spectral data"]),
                    dcc.Upload(
                        id="upload-spectra", className="upload-box", multiple=True,
                        children=html.Div([
                            html.Span("Drop or choose one or more spectral files"),
                            html.Small("SpectraGryph .dat, Avantes .Abs8, TSV, or CSV"),
                        ]),
                    ),
                    html.Label("Text input format"),
                    dcc.Dropdown(
                        id="input-format", value="spectragryph", clearable=False,
                        options=[
                            {"label": "SpectraGryph matrix (default)", "value": "spectragryph"},
                            {"label": "Avantes AvaSoft 8 (.Abs8)", "value": "avantes_abs8"},
                            {"label": "TSV: wavelength + spectra", "value": "tsv"},
                            {"label": "CSV: wavelength + spectra", "value": "csv"},
                        ],
                    ),
                    html.Button("Clear all spectra", id="clear-dataset",
                                className="button button-secondary", disabled=True),
                    html.Div(id="load-message", className="message"),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("2 · Range", className="step-label"),
                                  "Wavelengths"]),
                    html.Div([
                        dcc.Input(id="wavelength-low", type="number", placeholder="Start (nm)", disabled=True),
                        dcc.Input(id="wavelength-high", type="number", placeholder="End (nm)", disabled=True),
                    ], className="input-row"),
                    html.Small("The preview, epsilon calculation, and export use this range."),
                    html.Details(open=False, className="nested-tool", children=[
                        html.Summary("Preprocess spectra"),
                        html.Div(className="interaction-bar", children=[
                            dcc.Checklist(
                                id="baseline-enabled", value=[], className="toggle-control",
                                options=[{"label": "Baseline", "value": "on"}],
                            ),
                            dcc.RadioItems(
                                id="smoothing-method", value="off",
                                className="segmented-control two-options",
                                options=[
                                    {"label": "Raw", "value": "off"},
                                    {"label": "SavGol", "value": "savgol"},
                                ],
                            ),
                        ]),
                        html.Label("Baseline interval (nm)"),
                        html.Div([
                            dcc.Input(id="baseline-low", type="number", placeholder="Start"),
                            dcc.Input(id="baseline-high", type="number", placeholder="End"),
                        ], className="input-row"),
                        html.Details(open=False, className="parameter-details", children=[
                            html.Summary("Smoothing parameters"),
                            html.Label("Savitzky–Golay: window (nm) / polynomial order"),
                            html.Div([
                                dcc.Input(id="savgol-window", type="number", value=5,
                                          min=0, step="any"),
                                dcc.Input(id="savgol-order", type="number", value=3,
                                          min=0, step=1),
                            ], className="input-row"),
                        ]),
                        html.Div(className="svd-control-row", children=[
                            dcc.Checklist(
                                id="svd-enabled", value=[], className="toggle-control",
                                options=[{"label": "SVD", "value": "on"}],
                            ),
                            dcc.Dropdown(id="svd-rank", clearable=False, disabled=True,
                                         placeholder="Load data for rank suggestion"),
                        ]),
                        html.Div(id="svd-message", className="message"),
                        html.P("Use SVD only for ordered time-based datasets. Do not use it on "
                               "independent repeat measurements: it mixes columns and can suppress "
                               "the real between-measurement variation.", className="svd-warning"),
                        html.P("Baseline and Savitzky–Golay are applied to each spectrum "
                               "independently; uploaded values remain unchanged.",
                               className="warning-copy"),
                        html.Div(id="smoothing-message", className="message"),
                    ]),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("3 · Beer–Lambert", className="step-label"),
                                  "Final concentrations"]),
                    html.P("Enter the final concentration of every measured solution directly "
                           "in mol/L. Path length is required in centimetres.",
                           className="helper-text"),
                    html.Div(id="concentration-parameters"),
                    html.Div(id="concentration-message", className="message"),
                ]),
                html.Section(className="panel export-panel", children=[
                    html.P("4 · Output", className="step-label"),
                    html.H2("Export absorptivity dataset"),
                    html.Button("Export epsilon TSV", id="export-epsilon",
                                className="button button-primary", disabled=True),
                    html.Small("Includes processed absorbance, each epsilon spectrum, mean, SD, and SEM."),
                    dcc.Download(id="download-epsilon"),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("5 · Optional", className="step-label"),
                                  "NMR-guided PSS subtraction"]),
                    html.P("Load one UV–Vis dataset containing the full irradiation sequence. "
                           "The first spectrum is used as pure reactant and the last spectrum "
                           "as the final PSS. Subtraction is performed in the shared normalized "
                           "space as (PSS − x · reactant) / (1 − x).",
                           className="helper-text"),
                    dcc.Upload(
                        id="nmr-upload", className="upload-box compact-upload", multiple=False,
                        children=html.Div([
                            html.Span("Drop one dataset containing reactant and PSS"),
                            html.Small("First spectrum = reactant · last spectrum = PSS"),
                        ]),
                    ),
                    html.Label("Text input format"),
                    dcc.Dropdown(
                        id="nmr-input-format", value="spectragryph", clearable=False,
                        options=[
                            {"label": "SpectraGryph matrix (default)", "value": "spectragryph"},
                            {"label": "TSV", "value": "tsv"},
                            {"label": "CSV", "value": "csv"},
                        ],
                    ),
                    html.Details(open=True, className="nested-tool", children=[
                        html.Summary("Preprocess reactant and PSS"),
                        html.Div(className="interaction-bar", children=[
                            dcc.Checklist(
                                id="nmr-baseline-enabled", value=[],
                                className="toggle-control",
                                options=[{"label": "Baseline", "value": "on"}],
                            ),
                            dcc.RadioItems(
                                id="nmr-smoothing-method", value="off",
                                className="segmented-control two-options",
                                options=[
                                    {"label": "Raw", "value": "off"},
                                    {"label": "SavGol", "value": "savgol"},
                                ],
                            ),
                        ]),
                        html.Label("Baseline interval (nm)"),
                        html.Div([
                            dcc.Input(id="nmr-baseline-low", type="number", placeholder="Start"),
                            dcc.Input(id="nmr-baseline-high", type="number", placeholder="End"),
                        ], className="input-row"),
                        html.Details(open=False, className="parameter-details", children=[
                            html.Summary("Smoothing parameters"),
                            html.Label("Savitzky–Golay: window (nm) / polynomial order"),
                            html.Div([
                                dcc.Input(id="nmr-savgol-window", type="number", value=5,
                                          min=0, step="any"),
                                dcc.Input(id="nmr-savgol-order", type="number", value=3,
                                          min=0, step=1),
                            ], className="input-row"),
                        ]),
                        html.Div(id="nmr-processing-message", className="message"),
                    ]),
                    html.Div(className="input-row", children=[
                        html.Div([html.Label("Reactant in final PSS (%)"),
                                  dcc.Input(id="nmr-reactant-percent", type="number",
                                            value=10, min=0, max=99.999, step="any")]),
                        html.Div([html.Label("NMR error (%)"),
                                  dcc.Input(id="nmr-error-percent", type="number",
                                            value=1, min=0, step="any")]),
                    ]),
                    html.Button("Clear NMR spectra", id="clear-nmr",
                                className="button button-secondary", disabled=True),
                    html.Div(id="nmr-load-message", className="message"),
                    html.Div(id="nmr-result-message", className="warning-copy"),
                    html.Div(className="export-mode-control", children=[
                        dcc.Checklist(
                            id="nmr-export-raw", value=[], className="toggle-control",
                            options=[{"label": "Keep negative values in product ε", "value": "on"}],
                        ),
                        html.Small("Default: the primary product ε column is clipped at zero. "
                                   "The raw audit column is always preserved."),
                    ]),
                    html.Button("Export NMR-derived epsilon TSV", id="export-nmr",
                                className="button button-accent", disabled=True),
                    dcc.Download(id="download-nmr"),
                ]),
            ]),
            html.Div(className="plot-column", children=[
                html.Section(className="plot-panel", children=[
                    dcc.Graph(id="epsilon-plot", figure=_empty(go, "Load spectral data to begin")),
                ]),
                html.Section(id="nmr-plot-panel", className="plot-panel", style={"display": "none"},
                             children=[dcc.Graph(id="nmr-plot", figure=_empty(go, "Load NMR spectra"))]),
                html.Section(className="panel result-summary", children=[
                    html.H2("Result"), html.Div(id="result-message", className="helper-text"),
                ]),
                html.Details(open=False, className="panel tool-details error-panel", children=[
                    html.Summary("Python errors"),
                    html.Pre(id="load-error"), html.Pre(id="preview-error"),
                    html.Pre(id="svd-error"), html.Pre(id="nmr-error"),
                    html.Pre(id="export-error"),
                    html.Small("An empty panel means no Python error is active."),
                ]),
            ]),
        ]),
        dcc.Store(id="dataset-store"),
        dcc.Store(id="epsilon-store"),
        dcc.Store(id="nmr-spectra-store"),
        dcc.Store(id="nmr-result-store"),
    ])

    @app.callback(
        Output("dataset-store", "data"), Output("load-message", "children"),
        Output("upload-spectra", "contents"), Output("upload-spectra", "filename"),
        Output("load-error", "children"),
        Input("upload-spectra", "contents"), Input("clear-dataset", "n_clicks"),
        State("upload-spectra", "filename"), State("input-format", "value"),
        prevent_initial_call=True,
    )
    def load(contents, _, filenames, format_name):
        if ctx.triggered_id == "clear-dataset":
            return None, "All spectra cleared.", None, None, ""
        try:
            if not contents:
                return None, "Choose one or more files to begin.", no_update, no_update, ""
            contents = contents if isinstance(contents, list) else [contents]
            filenames = filenames if isinstance(filenames, list) else [filenames]
            restored = _try_load_autoqy_export(contents, filenames)
            if restored is not None:
                result, labels = restored
                dataset = SpectralDataset(
                    result.wavelengths, np.arange(len(labels), dtype=float),
                    result.absorbance, "autoqy_epsilon", 0,
                )
                return (
                    _pack(dataset, labels, filenames, result.concentrations_m,
                          result.path_lengths_cm),
                    f"Restored {len(labels)} processed absorbance spectrum/spectra, "
                    "concentrations, and path lengths from an AutoQY epsilon TSV.",
                    no_update, no_update, "",
                )
            loaded = []
            for content, filename in zip(contents, filenames):
                payload = base64.b64decode(content.split(",", 1)[1])
                selected_format = ("avantes_abs8" if Path(filename or "").suffix.lower() == ".abs8"
                                   else format_name)
                loaded.append((load_spectral_bytes(payload, selected_format), filename))
            dataset, labels, resampled = _combine_loaded(loaded)
            missing = sum(item.interpolated_values for item, _ in loaded)
            notes = []
            if missing:
                notes.append(f"interpolated {missing} non-finite detector value(s)")
            if resampled:
                notes.append(f"resampled {resampled} spectrum/spectra to the first common grid")
            note = f" ({'; '.join(notes)})" if notes else ""
            return (
                _pack(dataset, labels, filenames),
                f"Loaded {len(filenames)} file(s), {len(labels)} spectrum/spectra, "
                f"and {len(dataset.wavelengths)} common wavelengths{note}.",
                no_update, no_update, "",
            )
        except Exception as error:
            return (None, "Could not load the selected files.", no_update, no_update,
                    f"Load error: {type(error).__name__}: {error}")

    @app.callback(
        Output("wavelength-low", "value"), Output("wavelength-high", "value"),
        Output("wavelength-low", "disabled"), Output("wavelength-high", "disabled"),
        Output("clear-dataset", "disabled"), Output("concentration-parameters", "children"),
        Input("dataset-store", "data"),
    )
    def dataset_controls(data):
        if not data:
            return None, None, True, True, True, []
        wavelengths = np.asarray(data["wavelengths"], float)
        return (float(wavelengths.min()), float(wavelengths.max()), False, False, False,
                _parameter_cards(html, dcc, data["labels"],
                                 data.get("concentrations"), data.get("path_lengths")))

    @app.callback(
        Output("svd-rank", "options"), Output("svd-rank", "value"),
        Output("svd-rank", "disabled"), Output("svd-message", "children"),
        Output("svd-error", "children"),
        Input("dataset-store", "data"), Input("wavelength-low", "value"),
        Input("wavelength-high", "value"), Input("baseline-enabled", "value"),
        Input("baseline-low", "value"), Input("baseline-high", "value"),
        Input("smoothing-method", "value"), Input("savgol-window", "value"),
        Input("savgol-order", "value"), Input("svd-enabled", "value"),
    )
    def configure_svd(data, wavelength_low, wavelength_high, baseline_enabled,
                      baseline_low, baseline_high, method, sg_width, sg_order,
                      svd_enabled):
        if not data:
            return [], None, True, "SVD is off.", ""
        try:
            dataset, _, processed, _ = _prepare_processing(
                data, wavelength_low, wavelength_high, baseline_enabled,
                baseline_low, baseline_high, method, sg_width, sg_order,
            )
            analysis = analyze_svd(SpectralDataset(
                dataset.wavelengths, dataset.coordinates, processed,
                dataset.source_format, dataset.interpolated_values,
            ))
            maximum = min(10, len(analysis.singular_values))
            options = [{"label": f"{rank} component{'s' if rank != 1 else ''}",
                        "value": rank} for rank in range(1, maximum + 1)]
            proposed = min(max(analysis.proposed_rank, 2), maximum)
            retained = 100 * analysis.cumulative_weights[proposed - 1]
            message = (f"Proposed: {proposed} component(s), retaining {retained:.6g}% of "
                       "squared singular-value weight. Inspect the result yourself.")
            return options, proposed, "on" not in (svd_enabled or []), message, ""
        except Exception as error:
            return [], None, True, "SVD rank unavailable.", f"SVD error: {type(error).__name__}: {error}"

    @app.callback(
        Output("epsilon-plot", "figure"), Output("result-message", "children"),
        Output("concentration-message", "children"), Output("smoothing-message", "children"),
        Output("epsilon-store", "data"), Output("export-epsilon", "disabled"),
        Output("preview-error", "children"),
        Input("dataset-store", "data"), Input("wavelength-low", "value"),
        Input("wavelength-high", "value"), Input("baseline-enabled", "value"),
        Input("baseline-low", "value"), Input("baseline-high", "value"),
        Input("smoothing-method", "value"), Input("savgol-window", "value"),
        Input("savgol-order", "value"), Input("svd-enabled", "value"),
        Input("svd-rank", "value"),
        Input({"type": "direct-concentration", "index": ALL}, "value"),
        Input({"type": "path-length", "index": ALL}, "value"),
    )
    def preview(data, wavelength_low, wavelength_high, baseline_enabled,
                baseline_low, baseline_high, method, sg_width, sg_order,
                svd_enabled, svd_rank, concentrations, path_lengths):
        if not data:
            return (_empty(go, "Load spectral data to begin"),
                    "No result yet.", "", "Smoothing is off.", None, True, "")
        try:
            dataset, original, processed, smoothing_message = _prepare_processing(
                data, wavelength_low, wavelength_high, baseline_enabled,
                baseline_low, baseline_high, method, sg_width, sg_order,
                svd_enabled, svd_rank,
            )
            concentration_data = _read_concentrations(
                len(data["labels"]), concentrations, path_lengths
            )
            if concentration_data is None:
                message = ("Enter the final concentration and path length for every "
                           "spectrum to calculate molar absorptivity.")
                return (_absorbance_figure(go, dataset, original, processed, data["labels"],
                                           method, svd_enabled, svd_rank),
                        "Absorbance preview; epsilon is waiting for concentration inputs.",
                        message, smoothing_message, None, True, "")
            concentrations, paths = concentration_data
            result = calculate_epsilon_statistics(
                dataset.wavelengths, processed, concentrations, paths
            )
            count = len(concentrations)
            if count > 1:
                result_message = (f"Mean molar absorptivity from {count} independent spectra; "
                                  "the shaded region is ±1 sample SD at each wavelength. "
                                  "SEM is also included in the export.")
            else:
                result_message = ("One epsilon spectrum calculated. At least two independent "
                                  "spectra are required to estimate wavelength-resolved SD.")
            negative_mean = int(np.count_nonzero(result.mean < 0))
            if count > 1:
                raw_lower = result.mean - np.abs(result.standard_deviation)
                negative_error = int(np.count_nonzero(raw_lower < 0))
            else:
                negative_error = 0
            if negative_mean or negative_error:
                result_message += (f" Warning: {negative_mean} mean value(s) and "
                                   f"{negative_error} lower error bound(s) are negative. "
                                   "Negative means remain visible; plotted/exported error "
                                   "bounds are constrained to zero.")
            concentration_message = "Concentrations: " + ", ".join(
                f"{value:.6g} M" for value in concentrations
            )
            return (
                _epsilon_figure(go, make_subplots, dataset, original, result,
                                data["labels"], method, svd_enabled, svd_rank),
                result_message, concentration_message, smoothing_message,
                _pack_epsilon(result, data["labels"]), False, "",
            )
        except Exception as error:
            return (_empty(go, "Preview unavailable; open Python errors below."),
                    "No valid result.", "", "", None, True,
                    f"Preview error: {type(error).__name__}: {error}")

    @app.callback(
        Output("download-epsilon", "data"), Output("export-error", "children"),
        Input("export-epsilon", "n_clicks"), State("epsilon-store", "data"),
        prevent_initial_call=True,
    )
    def download(_, data):
        try:
            if not data:
                raise ValueError("Calculate epsilon before exporting")
            result, labels = _unpack_epsilon(data)
            return dict(content=export_epsilon_tsv(result, labels),
                        filename="epsilon-spectra.tsv"), ""
        except Exception as error:
            return no_update, f"Export error: {type(error).__name__}: {error}"

    @app.callback(
        Output("nmr-spectra-store", "data"), Output("nmr-load-message", "children"),
        Output("nmr-upload", "contents"), Output("nmr-upload", "filename"),
        Output("clear-nmr", "disabled"), Output("nmr-error", "children"),
        Input("nmr-upload", "contents"), Input("clear-nmr", "n_clicks"),
        State("nmr-upload", "filename"), State("nmr-input-format", "value"),
        prevent_initial_call=True,
    )
    def load_nmr(contents, _, filenames, format_name):
        if ctx.triggered_id == "clear-nmr":
            return None, "NMR spectra cleared.", None, None, True, ""
        try:
            if not contents:
                return None, "Choose one dataset containing at least two spectra.", no_update, no_update, True, ""
            if isinstance(contents, list):
                if len(contents) != 1:
                    raise ValueError("Select one dataset containing both reactant and PSS")
                contents = contents[0]
            if isinstance(filenames, list):
                filenames = filenames[0]
            payload = base64.b64decode(contents.split(",", 1)[1])
            dataset = load_spectral_bytes(payload, format_name)
            if dataset.absorbance.shape[1] < 2:
                raise ValueError("The NMR dataset must contain at least two spectra")
            return ({
                "filename": filenames,
                "reactant_wavelengths": dataset.wavelengths.tolist(),
                "reactant_absorbance": dataset.absorbance[:, 0].tolist(),
                "pss_wavelengths": dataset.wavelengths.tolist(),
                "pss_absorbance": dataset.absorbance[:, -1].tolist(),
            }, f"Loaded {dataset.absorbance.shape[1]} spectra from {filenames}; "
               "using the first as reactant and the last as PSS.",
               no_update, no_update, False, "")
        except Exception as error:
            return (None, "Could not load NMR subtraction files.", no_update, no_update, True,
                    f"NMR load error: {type(error).__name__}: {error}")

    @app.callback(
        Output("nmr-plot", "figure"), Output("nmr-plot-panel", "style"),
        Output("nmr-result-message", "children"), Output("nmr-result-store", "data"),
        Output("export-nmr", "disabled"), Output("nmr-result-message", "className"),
        Output("nmr-processing-message", "children"),
        Input("nmr-spectra-store", "data"), Input("epsilon-store", "data"),
        Input("nmr-reactant-percent", "value"), Input("nmr-error-percent", "value"),
        Input("nmr-baseline-enabled", "value"),
        Input("nmr-baseline-low", "value"), Input("nmr-baseline-high", "value"),
        Input("nmr-smoothing-method", "value"),
        Input("nmr-savgol-window", "value"), Input("nmr-savgol-order", "value"),
    )
    def calculate_nmr(nmr_data, epsilon_data, reactant_percent, error_percent,
                      baseline_enabled, baseline_low, baseline_high, smoothing_method,
                      sg_width, sg_order):
        if not nmr_data:
            return (_empty(go, "Load reactant and PSS spectra"), {"display": "none"},
                    "", None, True, "status-message", "Preprocessing is off.")
        if not epsilon_data:
            return (_empty(go, "Calculate or reload reactant epsilon first"), {"display": "block"},
                    "Reactant epsilon is required before NMR subtraction.", None, True,
                    "status-message status-warning", "Preprocessing is waiting for reactant epsilon.")
        try:
            epsilon_result, _ = _unpack_epsilon(epsilon_data)
            low = max(float(np.min(nmr_data["reactant_wavelengths"])),
                      float(np.min(nmr_data["pss_wavelengths"])))
            high = min(float(np.max(nmr_data["reactant_wavelengths"])),
                       float(np.max(nmr_data["pss_wavelengths"])))
            mask = ((epsilon_result.wavelengths >= low) &
                    (epsilon_result.wavelengths <= high))
            if np.count_nonzero(mask) < 2:
                raise ValueError("NMR and reactant epsilon datasets do not share a wavelength range")
            epsilon_result = _subset_epsilon(epsilon_result, mask)
            target = epsilon_result.wavelengths
            raw_reactant = _interpolate_to_axis(
                nmr_data["reactant_wavelengths"], nmr_data["reactant_absorbance"], target
            )
            raw_pss = _interpolate_to_axis(
                nmr_data["pss_wavelengths"], nmr_data["pss_absorbance"], target
            )
            reactant, pss, processing_message = _process_nmr_pair(
                target, raw_reactant, raw_pss, baseline_enabled,
                baseline_low, baseline_high, smoothing_method,
                sg_width, sg_order,
            )
            result = reconstruct_product_from_nmr(
                target, reactant, pss, epsilon_result,
                reactant_percent, error_percent,
            )
            minimum = float(np.min(result.product))
            figure = _nmr_figure(
                go, make_subplots, result, raw_reactant, raw_pss,
                preprocessing_changed=not (
                    np.array_equal(raw_reactant, reactant) and np.array_equal(raw_pss, pss)
                ),
            )
            if minimum < -500:
                message = (f"HARD STOP: reconstructed product ε reaches {minimum:.4g} M⁻¹ cm⁻¹, "
                           "below the −500 tolerance. Review baseline, smoothing, PSS composition, "
                           "and spectral alignment. Export is disabled.")
                return (figure, {"display": "block"}, message, None, True,
                        "status-message status-stop", processing_message)
            warnings = []
            if result.negative_product_points:
                warnings.append(
                    f"{result.negative_product_points} product point(s) lie between "
                    f"{minimum:.4g} and 0 M⁻¹ cm⁻¹"
                )
            if result.negative_bound_points:
                warnings.append(
                    f"{result.negative_bound_points} raw uncertainty bound(s) cross zero"
                )
            if warnings:
                message = ("Soft warning: " + "; ".join(warnings) +
                           ". These remain within the −500 tolerance. Negative means stay "
                           "visible; uncertainty bounds are constrained to zero.")
                status = "status-message status-warning"
            else:
                message = ("Product epsilon reconstructed from normalized "
                           "(PSS − x · reactant) / (1 − x). The band combines the "
                           "reactant epsilon scale SD and selected NMR error by min/max propagation.")
                status = "status-message status-ok"
            return (figure, {"display": "block"}, message, _pack_nmr(result), False,
                    status, processing_message)
        except Exception as error:
            return (_empty(go, "NMR reconstruction unavailable"), {"display": "block"},
                    f"NMR reconstruction error: {type(error).__name__}: {error}", None, True,
                    "status-message status-stop", "Preprocessing could not be applied.")

    @app.callback(
        Output("download-nmr", "data"),
        Input("export-nmr", "n_clicks"), State("nmr-result-store", "data"),
        State("nmr-export-raw", "value"),
        prevent_initial_call=True,
    )
    def download_nmr(_, data, preserve_negative):
        if not data:
            return no_update
        return dict(content=export_nmr_subtraction_tsv(
                        _unpack_nmr(data), "on" in (preserve_negative or [])),
                    filename="nmr-derived-product-epsilon.tsv")

    return app


def _parameter_cards(html, dcc, labels, concentrations=None, path_lengths=None):
    cards = []
    for index, label in enumerate(labels):
        cards.append(html.Div(className="spectrum-card", children=[
            html.H3(label),
            html.Div(className="parameter-grid", children=[
                html.Div([html.Label("Final concentration (mol/L)"),
                          dcc.Input(
                              id={"type": "direct-concentration", "index": index},
                              type="number", min=0, step="any", placeholder="mol/L",
                              value=(concentrations[index]
                                     if concentrations is not None else None),
                          )]),
                html.Div([html.Label("Path length (cm)"),
                          dcc.Input(id={"type": "path-length", "index": index},
                                    type="number", min=0, step="any",
                                    value=(path_lengths[index]
                                           if path_lengths is not None else 1.0))]),
            ]),
        ]))
    return cards


def _combine_loaded(loaded):
    if not loaded:
        raise ValueError("No spectral files were supplied")
    low = max(float(np.min(dataset.wavelengths)) for dataset, _ in loaded)
    high = min(float(np.max(dataset.wavelengths)) for dataset, _ in loaded)
    if low >= high:
        raise ValueError("The selected files do not share a wavelength range")
    first = np.sort(np.unique(np.asarray(loaded[0][0].wavelengths, float)))
    target = first[(first >= low) & (first <= high)]
    if len(target) < 2:
        raise ValueError("The files share fewer than two wavelength values")
    columns, labels = [], []
    resampled = 0
    for dataset, filename in loaded:
        order = np.argsort(dataset.wavelengths)
        wavelengths, unique_index = np.unique(
            np.asarray(dataset.wavelengths, float)[order], return_index=True
        )
        values = np.asarray(dataset.absorbance, float)[order][unique_index]
        exact = len(wavelengths) == len(target) and np.allclose(wavelengths, target)
        stem = Path(filename or "spectrum").stem
        for index in range(values.shape[1]):
            columns.append(values[:, index] if exact else np.interp(
                target, wavelengths, values[:, index]
            ))
            coordinate = dataset.coordinates[index]
            labels.append(stem if values.shape[1] == 1 else f"{stem} [{coordinate:g}]")
            resampled += int(not exact)
    combined = SpectralDataset(target, np.arange(len(columns), dtype=float),
                               np.column_stack(columns), "combined", 0)
    return combined, _display_unique(labels), resampled


def _read_concentrations(count, concentrations, path_lengths):
    if any(len(values or []) != count for values in (concentrations, path_lengths)):
        return None
    parsed_concentrations, parsed_paths = [], []
    for concentration, path_length in zip(concentrations, path_lengths):
        if concentration is None or path_length is None:
            return None
        parsed_concentrations.append(float(concentration))
        parsed_paths.append(float(path_length))
    return np.asarray(parsed_concentrations), np.asarray(parsed_paths)


def _prepare_processing(data, wavelength_low, wavelength_high, baseline_enabled,
                        baseline_low, baseline_high, method, sg_width, sg_order,
                        svd_enabled=None, svd_rank=None):
    raw = _unpack(data)
    original = select_wavelengths(raw, _wavelength_interval(wavelength_low, wavelength_high))
    baselined = baseline_spectra(raw, _interval(baseline_enabled, baseline_low, baseline_high))
    baseline_dataset = SpectralDataset(raw.wavelengths, raw.coordinates, baselined,
                                       raw.source_format, raw.interpolated_values)
    selected = select_wavelengths(
        baseline_dataset, _wavelength_interval(wavelength_low, wavelength_high)
    )
    sg_points = (savgol_window_points(selected.wavelengths, sg_width, sg_order)
                 if method == "savgol" else 11)
    processed = smooth_reconstruction(
        selected.absorbance, method, savgol_window=sg_points,
        savgol_order=sg_order,
    )
    rms = float(np.sqrt(np.mean((processed - selected.absorbance) ** 2)))
    if method == "savgol":
        spacing = float(np.median(np.abs(np.diff(selected.wavelengths))))
        message = (f"Savitzky–Golay: {sg_points} points "
                   f"(~{sg_points * spacing:.3g} nm); smoothing RMS {rms:.4g} absorbance.")
    else:
        message = "Smoothing is off."
    if "on" in (baseline_enabled or []):
        message = "Baseline correction applied independently. " + message
    if "on" in (svd_enabled or []):
        if svd_rank is None:
            raise ValueError("Select an SVD component count")
        analysis = analyze_svd(SpectralDataset(
            selected.wavelengths, selected.coordinates, processed,
            selected.source_format, selected.interpolated_values,
        ))
        processed = analysis.reconstruct(int(svd_rank))
        retained = 100 * analysis.cumulative_weights[int(svd_rank) - 1]
        message += (f" SVD rank {int(svd_rank)} retains {retained:.6g}% of squared "
                    "singular-value weight.")
    else:
        message += " SVD is off."
    return selected, original.absorbance, processed, message


def _interval(enabled, low, high):
    if "on" not in (enabled or []):
        return None
    if low is None or high is None:
        raise ValueError("Enter both baseline interval limits")
    return float(low), float(high)


def _wavelength_interval(low, high):
    if low is None or high is None:
        raise ValueError("Enter both wavelength range limits")
    return float(low), float(high)


def _pack(dataset, labels, filenames, concentrations=None, path_lengths=None):
    packed = {
        "filenames": list(filenames), "format": dataset.source_format,
        "labels": labels, "wavelengths": dataset.wavelengths.tolist(),
        "coordinates": dataset.coordinates.tolist(),
        "absorbance": dataset.absorbance.tolist(),
        "interpolated_values": dataset.interpolated_values,
    }
    if concentrations is not None:
        packed["concentrations"] = np.asarray(concentrations, float).tolist()
    if path_lengths is not None:
        packed["path_lengths"] = np.asarray(path_lengths, float).tolist()
    return packed


def _unpack(data):
    return SpectralDataset(
        np.asarray(data["wavelengths"], float), np.asarray(data["coordinates"], float),
        np.asarray(data["absorbance"], float), data["format"],
        data.get("interpolated_values", 0),
    )


def _pack_epsilon(result, labels):
    return {
        "labels": labels, "wavelengths": result.wavelengths.tolist(),
        "absorbance": result.absorbance.tolist(),
        "concentrations": result.concentrations_m.tolist(),
        "path_lengths": result.path_lengths_cm.tolist(),
        "individual": result.individual.tolist(), "mean": result.mean.tolist(),
        "sd": result.standard_deviation.tolist(), "sem": result.standard_error.tolist(),
    }


def _unpack_epsilon(data):
    result = EpsilonResult(
        np.asarray(data["wavelengths"]), np.asarray(data["absorbance"]),
        np.asarray(data["concentrations"]), np.asarray(data["path_lengths"]),
        np.asarray(data["individual"]), np.asarray(data["mean"]),
        np.asarray(data["sd"]), np.asarray(data["sem"]),
    )
    return result, data["labels"]


def _try_load_autoqy_export(contents, filenames):
    if len(contents) != 1 or Path(filenames[0] or "").suffix.lower() != ".tsv":
        return None
    payload = base64.b64decode(contents[0].split(",", 1)[1])
    try:
        text = payload.decode("utf-8-sig")
        return load_epsilon_tsv(text)
    except (UnicodeDecodeError, ValueError):
        return None


def _interpolate_to_axis(wavelengths, values, target):
    wavelengths = np.asarray(wavelengths, float)
    values = np.asarray(values, float)
    target = np.asarray(target, float)
    order = np.argsort(wavelengths)
    wavelengths, unique = np.unique(wavelengths[order], return_index=True)
    values = values[order][unique]
    if target[0] < wavelengths[0] or target[-1] > wavelengths[-1]:
        raise ValueError("NMR spectra must span the full reactant epsilon wavelength range")
    return np.interp(target, wavelengths, values)


def _process_nmr_pair(wavelengths, reactant, pss, baseline_enabled,
                      baseline_low, baseline_high, smoothing_method,
                      sg_width, sg_order):
    values = np.column_stack([reactant, pss])
    dataset = SpectralDataset(
        np.asarray(wavelengths, float), np.array([0.0, 1.0]), values, "nmr_pair", 0
    )
    interval = _interval(baseline_enabled, baseline_low, baseline_high)
    processed = baseline_spectra(dataset, interval)
    sg_points = (savgol_window_points(dataset.wavelengths, sg_width, sg_order)
                 if smoothing_method == "savgol" else 11)
    processed = smooth_reconstruction(
        processed, smoothing_method, savgol_window=sg_points,
        savgol_order=sg_order,
    )
    actions = []
    if interval is not None:
        actions.append(f"baseline {interval[0]:g}–{interval[1]:g} nm")
    if smoothing_method == "savgol":
        actions.append(f"Savitzky–Golay ({sg_points} points, order {int(sg_order)})")
    rms = float(np.sqrt(np.mean((processed - values) ** 2)))
    message = ("Applied " + " + ".join(actions) + f" · total change RMS {rms:.4g} absorbance."
               if actions else "NMR preprocessing is off; raw first/last spectra are used.")
    return processed[:, 0], processed[:, 1], message


def _subset_epsilon(result, mask):
    return EpsilonResult(
        result.wavelengths[mask], result.absorbance[mask],
        result.concentrations_m, result.path_lengths_cm,
        result.individual[mask], result.mean[mask],
        result.standard_deviation[mask], result.standard_error[mask],
    )


def _pack_nmr(result):
    return {field: np.asarray(getattr(result, field)).tolist()
            for field in ("wavelengths", "normalized_reactant", "normalized_pss",
                          "reconstructed_reactant", "reconstructed_pss", "product",
                          "product_lower", "product_upper", "product_error_lower",
                          "product_error_upper")} | {
        "negative_product_points": result.negative_product_points,
        "negative_bound_points": result.negative_bound_points,
    }


def _unpack_nmr(data):
    return NMRSubtractionResult(
        *(np.asarray(data[field], float) for field in
          ("wavelengths", "normalized_reactant", "normalized_pss",
           "reconstructed_reactant", "reconstructed_pss", "product",
           "product_lower", "product_upper", "product_error_lower",
           "product_error_upper")),
        int(data["negative_product_points"]), int(data["negative_bound_points"]),
    )


def _display_unique(labels):
    totals = {}
    result = []
    for label in labels:
        totals[label] = totals.get(label, 0) + 1
        result.append(label if totals[label] == 1 else f"{label} ({totals[label]})")
    return result


def _absorbance_figure(go, dataset, original, processed, labels, method,
                       svd_enabled=None, svd_rank=None):
    figure = go.Figure()
    colors = _colors()
    changed = not np.allclose(original, processed)
    for index, label in enumerate(labels):
        color = colors[index % len(colors)]
        if changed:
            figure.add_trace(go.Scatter(
                x=dataset.wavelengths, y=original[:, index], mode="lines",
                line={"color": "rgba(90,96,108,.22)", "width": 1},
                showlegend=index == 0, name="Uploaded absorbance", hoverinfo="skip",
            ))
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=processed[:, index], mode="lines",
            line={"color": color, "width": 1.5}, name=label,
        ))
    figure.update_yaxes(title_text="Absorbance")
    figure.update_xaxes(title_text="Wavelength (nm)")
    figure.update_layout(title={"text": _processing_title(
        method, svd_enabled, svd_rank), "x": 0.02})
    return _style(figure, 520)


def _epsilon_figure(go, make_subplots, dataset, original, result, labels, method,
                    svd_enabled=None, svd_rank=None):
    figure = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        row_heights=[0.34, 0.66],
        subplot_titles=("Processed absorbance spectra", "Molar absorptivity"),
    )
    colors = _colors()
    for index, label in enumerate(labels):
        color = colors[index % len(colors)]
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=result.absorbance[:, index], mode="lines",
            line={"color": color, "width": 1.3}, name=label,
            legendgroup=f"spectrum-{index}",
        ), row=1, col=1)
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=result.individual[:, index], mode="lines",
            line={"color": color, "width": 1}, opacity=0.42,
            name=f"ε · {label}", legendgroup=f"spectrum-{index}", showlegend=False,
        ), row=2, col=1)
    if result.individual.shape[1] > 1:
        lower, upper = nonnegative_error_bounds(result.mean, result.standard_deviation)
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=lower, mode="lines",
            line={"width": 0}, hoverinfo="skip", showlegend=False,
            legendgroup="epsilon-error",
        ), row=2, col=1)
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=upper, mode="lines", fill="tonexty",
            fillcolor="rgba(214,123,54,.24)", line={"width": 0},
            name="Mean ± sample SD", legendgroup="epsilon-error",
        ), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=dataset.wavelengths, y=result.mean, mode="lines",
        line={"color": "#173b4c", "width": 2.4}, name="Mean ε",
    ), row=2, col=1)
    figure.update_yaxes(title_text="Absorbance", row=1, col=1)
    figure.update_yaxes(title_text="ε (M⁻¹ cm⁻¹)", row=2, col=1)
    figure.update_xaxes(title_text="Wavelength (nm)", row=2, col=1)
    figure.update_layout(title={"text": _processing_title(
        method, svd_enabled, svd_rank), "x": 0.02})
    return _style(figure, 690)


def _nmr_figure(go, make_subplots, result, raw_reactant=None, raw_pss=None,
                preprocessing_changed=False):
    figure = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        row_heights=[0.34, 0.66],
        subplot_titles=("Normalized UV–Vis inputs", "NMR-derived product absorptivity"),
    )
    if preprocessing_changed:
        raw_peak = float(np.max(raw_reactant))
        if raw_peak > 0:
            figure.add_trace(go.Scatter(
                x=result.wavelengths, y=np.asarray(raw_reactant) / raw_peak,
                mode="lines", name="Reactant raw", hoverinfo="skip",
                line={"color": "rgba(45,111,142,.25)", "width": 1}), row=1, col=1)
            figure.add_trace(go.Scatter(
                x=result.wavelengths, y=np.asarray(raw_pss) / raw_peak,
                mode="lines", name="PSS raw", hoverinfo="skip",
                line={"color": "rgba(214,123,54,.25)", "width": 1}), row=1, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.normalized_reactant, mode="lines",
        name="Reactant processed", line={"color": "#2d6f8e"}), row=1, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.normalized_pss, mode="lines",
        name="PSS processed", line={"color": "#d67b36"}), row=1, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.product_lower, mode="lines",
        line={"width": 0}, hoverinfo="skip", showlegend=False,
        legendgroup="product-error"), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.product_upper, mode="lines", fill="tonexty",
        fillcolor="rgba(42,157,143,.24)", line={"width": 0},
        name="Propagation envelope", legendgroup="product-error"), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.reconstructed_reactant, mode="lines",
        name="Reactant ε", line={"color": "#2d6f8e", "dash": "dot"}), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.reconstructed_pss, mode="lines",
        name="PSS reconstructed ε", line={"color": "#d67b36", "dash": "dash"}), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.product, mode="lines",
        name="Product ε", line={"color": "#173b4c", "width": 2.4}), row=2, col=1)
    figure.add_hline(y=0, line={"color": "#b85c4a", "dash": "dash", "width": 1}, row=2, col=1)
    figure.update_yaxes(title_text="Normalized absorbance", row=1, col=1)
    figure.update_yaxes(title_text="ε (M⁻¹ cm⁻¹)", row=2, col=1)
    figure.update_xaxes(title_text="Wavelength (nm)", row=2, col=1)
    figure.update_layout(title={"text": "Weighted PSS subtraction", "x": 0.02})
    return _style(figure, 650)


def _processing_title(method, svd_enabled=None, svd_rank=None):
    smoothing = {"off": "raw spectra", "savgol": "Savitzky–Golay"}[method]
    svd = (f" · rank-{int(svd_rank)} SVD"
           if "on" in (svd_enabled or []) and svd_rank is not None else " · SVD off")
    return f"Spectral treatment · {smoothing}{svd}"


def _colors():
    return ["#2d6f8e", "#d67b36", "#2a9d8f", "#8a5a9b", "#b85c4a",
            "#527a4f", "#725c42", "#467aa8"]


def _empty(go, message):
    figure = go.Figure()
    figure.add_annotation(text=message, showarrow=False)
    return _style(figure, 430)


def _style(figure, height):
    figure.update_layout(
        template="plotly_white", height=height,
        margin=dict(l=68, r=24, t=72, b=52), hovermode="closest",
        legend={"orientation": "v", "y": 0.99, "yanchor": "top",
                "x": 0.99, "xanchor": "right", "bgcolor": "rgba(255,255,255,.92)",
                "bordercolor": "#cfd1d4", "borderwidth": 1},
    )
    return figure


def run_server(host="127.0.0.1", port=8051, open_browser=True):
    app = create_app()
    if open_browser:
        Timer(1, lambda: webbrowser.open(f"http://{host}:{port}")).start()
    app.run(host=host, port=port, debug=False)


def main(argv=None):
    parser = ArgumentParser(prog="autoqy-smoother-gui")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8051, type=int)
    parser.add_argument("--no-open", action="store_true")
    args = parser.parse_args(argv)
    run_server(args.host, args.port, not args.no_open)


if __name__ == "__main__":
    main()
