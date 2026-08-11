"""Browser GUI for inspecting and applying SVD spectral denoising."""

from argparse import ArgumentParser
import base64
from pathlib import Path
from threading import Timer
import webbrowser

import numpy as np

from ..smoother import (SpectralDataset, analyze_svd, baseline_spectra,
                        export_smoothed_text, load_spectral_text,
                        savgol_window_points, select_wavelengths,
                        smooth_reconstruction)


WARNING = ("The proposed rank is advisory. Inspect the component spectra and weights: "
           "baseline drift, correlated noise, and weak chemical signals can each form components.")


def create_app():
    try:
        from dash import Dash, Input, Output, State, ctx, dcc, html, no_update
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as error:
        raise RuntimeError("The smoother GUI requires the 'power-gui' optional dependencies") from error

    assets = Path(__file__).parents[1] / "assets"
    app = Dash(__name__, assets_folder=str(assets))
    app.title = "AutoQY Spectral smoother"
    app.layout = html.Div(className="app-shell", children=[
        html.Header(className="app-header", children=[html.Div([
            html.P("AUTOQY CORE", className="eyebrow"), html.H1("Spectral smoother"),
            html.P("Denoise wavelength-resolved datasets and export an inspected reconstruction."),
        ]), html.Div("Local session", className="local-badge")]),
        html.Main(className="workspace", children=[
            html.Aside(className="control-column smoother-controls", children=[
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("1 · Data", className="step-label"), "Dataset"]),
                    dcc.Upload(id="upload-spectra", className="upload-box", multiple=False,
                               children=html.Div([html.Span("Drop or choose a spectral dataset"),
                                                  html.Small("SpectraGryph .dat, TSV, or CSV")])),
                    html.Label("Input format"),
                    dcc.Dropdown(id="input-format", value="spectragryph", clearable=False,
                                 options=[{"label": "SpectraGryph matrix (default)", "value": "spectragryph"},
                                          {"label": "TSV: wavelength + spectra", "value": "tsv"},
                                          {"label": "CSV: wavelength + spectra", "value": "csv"}]),
                    html.Button("Clear dataset", id="clear-dataset", className="button button-secondary", disabled=True),
                    html.Div(id="load-message", className="message"),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("2 · Wavelengths", className="step-label"), "Analysis range"]),
                    html.Div([dcc.Input(id="wavelength-low", type="number", placeholder="Start (nm)", disabled=True),
                              dcc.Input(id="wavelength-high", type="number", placeholder="End (nm)", disabled=True)],
                             className="input-row"),
                    html.Small("Only this wavelength range is decomposed, previewed, and exported."),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("3 · Optional", className="step-label"), "Baseline correction"]),
                    dcc.Checklist(id="baseline-enabled", options=[{"label": " Apply baseline", "value": "on"}], value=[]),
                    html.Div([dcc.Input(id="baseline-low", type="number", placeholder="Start (nm)"),
                              dcc.Input(id="baseline-high", type="number", placeholder="End (nm)")], className="input-row"),
                    html.Small("Each spectrum is shifted by its own mean within this wavelength interval."),
                ]),
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("4 · Primary", className="step-label"), "Spectral smoothing"]),
                    dcc.Dropdown(id="smoothing-method", value="savgol", clearable=False,
                                 options=[{"label": "Off", "value": "off"},
                                          {"label": "Savitzky–Golay", "value": "savgol"},
                                          {"label": "Whittaker–Eilers", "value": "whittaker"}]),
                    html.Label("Savitzky–Golay: window (nm) / polynomial order"),
                    html.Div([dcc.Input(id="savgol-window", type="number", value=5, min=0, step="any"),
                              dcc.Input(id="savgol-order", type="number", value=3, min=0, step=1)],
                             className="input-row"),
                    html.Label("Whittaker–Eilers: λ / difference order"),
                    html.Div([dcc.Input(id="whittaker-strength", type="number", value=1000, min=0, step="any"),
                              dcc.Input(id="whittaker-order", type="number", value=2, min=1, step=1)],
                             className="input-row"),
                    html.Div(id="smoothing-message", className="message"),
                    html.P("Applied before SVD. Original uploaded data remains unchanged.",
                           className="warning-copy"),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("5 · Optional", className="step-label"), "SVD reduction"]),
                    dcc.Dropdown(id="rank", clearable=False),
                    html.Div(id="rank-message", className="message"), html.P(WARNING, className="warning-copy"),
                ]),
                html.Section(className="panel export-panel", children=[
                    html.P("6 · Output", className="step-label"), html.H2("Export dataset"),
                    html.Button("Export smoothed dataset", id="export-smooth", className="button button-primary", disabled=True),
                    dcc.Download(id="download-smooth"),
                ]),
            ]),
            html.Div(className="plot-column", children=[
                html.Section(className="plot-panel", children=[dcc.Graph(id="reconstruction-plot", figure=_empty(go, "Load a dataset to begin"))]),
                html.Section(className="panel error-panel", children=[
                    html.H2("Python errors"),
                    html.Pre(id="load-error"), html.Pre(id="analysis-error"),
                    html.Pre(id="preview-error"), html.Pre(id="export-error"),
                    html.Small("Callback and parsing errors appear here. An empty window means no Python error is active."),
                ]),
            ]),
        ]),
        dcc.Store(id="dataset-store"),
    ])

    @app.callback(Output("dataset-store", "data"), Output("load-message", "children"),
                  Output("upload-spectra", "contents"), Output("upload-spectra", "filename"),
                  Output("load-error", "children"),
                  Input("upload-spectra", "contents"), Input("clear-dataset", "n_clicks"),
                  State("upload-spectra", "filename"),
                  State("input-format", "value"), prevent_initial_call=True)
    def load(contents, _, filename, format_name):
        if ctx.triggered_id == "clear-dataset":
            return None, "Dataset cleared.", None, None, ""
        try:
            if not contents:
                return None, "Choose a dataset to begin.", no_update, no_update, ""
            text = base64.b64decode(contents.split(",", 1)[1]).decode("utf-8-sig")
            dataset = load_spectral_text(text, format_name)
            note = (f" Interpolated {dataset.interpolated_values} missing absorbance value(s) along wavelength."
                    if dataset.interpolated_values else "")
            return (_pack(dataset, filename),
                    f"Loaded {len(dataset.coordinates)} spectra × {len(dataset.wavelengths)} wavelengths.{note}",
                    no_update, no_update, "")
        except Exception as error:
            message = f"Load error: {type(error).__name__}: {error}"
            return None, "Could not load dataset.", no_update, no_update, message

    @app.callback(Output("wavelength-low", "value"), Output("wavelength-high", "value"),
                  Output("wavelength-low", "disabled"), Output("wavelength-high", "disabled"),
                  Output("clear-dataset", "disabled"),
                  Input("dataset-store", "data"))
    def wavelength_controls(data):
        if not data:
            return None, None, True, True, True
        wavelengths = np.asarray(data["wavelengths"])
        low, high = float(np.min(wavelengths)), float(np.max(wavelengths))
        return low, high, False, False, False

    @app.callback(Output("rank", "options"), Output("rank", "value"),
                  Output("export-smooth", "disabled"),
                  Output("analysis-error", "children"), Output("smoothing-message", "children"),
                  Input("dataset-store", "data"), Input("wavelength-low", "value"),
                  Input("wavelength-high", "value"),
                  Input("baseline-enabled", "value"),
                  Input("baseline-low", "value"), Input("baseline-high", "value"),
                  Input("smoothing-method", "value"), Input("savgol-window", "value"),
                  Input("savgol-order", "value"), Input("whittaker-strength", "value"),
                  Input("whittaker-order", "value"))
    def analyze(data, wavelength_low, wavelength_high, enabled, low, high,
                method, sg_width, sg_order, wh_strength, wh_order):
        if not data:
            return [], None, True, "", ""
        try:
            _, _, _, result, smoothing_message = _prepare_processing(
                data, wavelength_low, wavelength_high, enabled, low, high,
                method, sg_width, sg_order, wh_strength, wh_order,
            )
            options = [{"label": "Off (spectral smoothing only)", "value": "off"}]
            options += [{"label": f"{rank} component{'s' if rank != 1 else ''}", "value": rank}
                        for rank in range(1, min(10, len(result.singular_values)) + 1)]
            selected = min(max(result.proposed_rank, 2), 10, len(result.singular_values))
            return options, selected, False, "", smoothing_message
        except Exception as error:
            message = f"Analysis error: {type(error).__name__}: {error}"
            return [], None, True, message, ""

    @app.callback(Output("reconstruction-plot", "figure"), Output("preview-error", "children"),
                  Output("rank-message", "children"),
                  Input("dataset-store", "data"), Input("wavelength-low", "value"),
                  Input("wavelength-high", "value"), Input("baseline-enabled", "value"), Input("baseline-low", "value"),
                  Input("baseline-high", "value"), Input("rank", "value"),
                  Input("smoothing-method", "value"), Input("savgol-window", "value"),
                  Input("savgol-order", "value"), Input("whittaker-strength", "value"),
                  Input("whittaker-order", "value"))
    def preview(data, wavelength_low, wavelength_high, enabled, low, high, rank,
                method, sg_window, sg_order, wh_strength, wh_order):
        if not data or rank is None:
            return _empty(go, "Select a rank to preview"), "", ""
        try:
            dataset, baseline_values, spectral_values, result, _ = _prepare_processing(
                data, wavelength_low, wavelength_high, enabled, low, high,
                method, sg_window, sg_order, wh_strength, wh_order,
            )
            final = spectral_values if rank == "off" else result.reconstruct(int(rank))
            if rank == "off":
                rank_message = "SVD reduction is off: 100% of the spectrally smoothed matrix is retained."
            else:
                retained = 100 * result.cumulative_weights[int(rank) - 1]
                rank_message = (f"{rank} component(s) retain {retained:.6g}% of squared "
                                "singular-value weight. Change the selection to compare interactively.")
            return _reconstruction_figure(
                go, make_subplots, dataset, baseline_values, spectral_values,
                final, rank, method
            ), "", rank_message
        except Exception as error:
            return (_empty(go, "Preview unavailable; see Python errors."),
                    f"Preview error: {type(error).__name__}: {error}", "")

    @app.callback(Output("download-smooth", "data"), Output("export-error", "children"),
                  Input("export-smooth", "n_clicks"),
                  State("dataset-store", "data"), State("baseline-enabled", "value"),
                  State("baseline-low", "value"), State("baseline-high", "value"),
                  State("rank", "value"), State("wavelength-low", "value"),
                  State("wavelength-high", "value"), State("smoothing-method", "value"),
                  State("savgol-window", "value"), State("savgol-order", "value"),
                  State("whittaker-strength", "value"), State("whittaker-order", "value"),
                  prevent_initial_call=True)
    def download(_, data, enabled, low, high, rank, wavelength_low, wavelength_high,
                 method, sg_window, sg_order, wh_strength, wh_order):
        try:
            dataset, _, spectral_values, result, _ = _prepare_processing(
                data, wavelength_low, wavelength_high, enabled, low, high,
                method, sg_window, sg_order, wh_strength, wh_order,
            )
            final = spectral_values if rank == "off" else result.reconstruct(int(rank))
            text = export_smoothed_text(dataset, final)
            original = Path(data["filename"] or "spectra.dat")
            suffix = "" if method == "off" else f"-{method}"
            svd_suffix = "" if rank == "off" else f"-svd-rank{rank}"
            return dict(content=text, filename=f"{original.stem}{suffix}{svd_suffix}{original.suffix or '.dat'}"), ""
        except Exception as error:
            return no_update, f"Export error: {type(error).__name__}: {error}"
    return app


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


def _prepare_processing(data, wavelength_low, wavelength_high, baseline_enabled,
                        baseline_low, baseline_high, method, sg_width, sg_order,
                        wh_strength, wh_order):
    dataset = select_wavelengths(
        _unpack(data), _wavelength_interval(wavelength_low, wavelength_high)
    )
    baseline_values = baseline_spectra(
        dataset, _interval(baseline_enabled, baseline_low, baseline_high)
    )
    sg_points = (savgol_window_points(dataset.wavelengths, sg_width, sg_order)
                 if method == "savgol" else 11)
    spectral_values = smooth_reconstruction(
        baseline_values, method, savgol_window=sg_points, savgol_order=sg_order,
        whittaker_strength=wh_strength, whittaker_order=wh_order,
    )
    processed_dataset = SpectralDataset(
        dataset.wavelengths, dataset.coordinates, spectral_values,
        dataset.source_format, dataset.interpolated_values,
    )
    result = analyze_svd(processed_dataset)
    rms = float(np.sqrt(np.mean((spectral_values - baseline_values) ** 2)))
    if method == "savgol":
        spacing = float(np.median(np.abs(np.diff(dataset.wavelengths))))
        message = (f"Savitzky–Golay uses {sg_points} points "
                   f"(~{sg_points * spacing:.3g} nm); smoothing RMS {rms:.4g} absorbance.")
    elif method == "whittaker":
        message = f"Whittaker–Eilers smoothing RMS: {rms:.4g} absorbance."
    else:
        message = "Primary spectral smoothing is off."
    return dataset, baseline_values, spectral_values, result, message


def _pack(dataset, filename):
    return {"filename": filename, "format": dataset.source_format,
            "wavelengths": dataset.wavelengths.tolist(), "coordinates": dataset.coordinates.tolist(),
            "absorbance": dataset.absorbance.tolist(), "interpolated_values": dataset.interpolated_values}


def _unpack(data):
    return SpectralDataset(np.array(data["wavelengths"]), np.array(data["coordinates"]),
                           np.array(data["absorbance"]), data["format"], data.get("interpolated_values", 0))


def _reconstruction_figure(go, make_subplots, dataset, baseline_values,
                           spectral_values, final, rank, method="off"):
    figure = make_subplots(rows=2, cols=1, shared_xaxes=True,
                           row_heights=[0.72, 0.28], vertical_spacing=0.08,
                           subplot_titles=("Processed spectra and final result", "Removed residual"))
    chosen = np.linspace(0, baseline_values.shape[1] - 1,
                         min(10, baseline_values.shape[1]), dtype=int)
    for position, index in enumerate(chosen):
        first = position == 0
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=baseline_values[:, index], mode="lines",
            name="Original/baselined", legendgroup="input", showlegend=first,
            line={"color": "rgba(80,90,110,.24)"}, hoverinfo="skip"), row=1, col=1)
        if method != "off":
            figure.add_trace(go.Scatter(
                x=dataset.wavelengths, y=spectral_values[:, index], mode="lines",
                name="Spectrally smoothed", legendgroup="spectral", showlegend=first,
                line={"color": "#2a9d8f", "dash": "dot"}), row=1, col=1)
        if rank != "off":
            figure.add_trace(go.Scatter(
                x=dataset.wavelengths, y=final[:, index], mode="lines",
                name="Final SVD", legendgroup="final", showlegend=first,
                line={"color": "#e76f51"}), row=1, col=1)
        figure.add_trace(go.Scatter(
            x=dataset.wavelengths, y=baseline_values[:, index] - final[:, index],
            mode="lines", showlegend=False, line={"color": "rgba(45,111,142,.5)"}),
            row=2, col=1)
    method_label = {"off": "No spectral smoothing", "savgol": "Savitzky–Golay",
                    "whittaker": "Whittaker–Eilers"}[method]
    rank_label = "SVD off" if rank == "off" else f"rank-{rank} SVD"
    removed_rms = float(np.sqrt(np.mean((baseline_values - final) ** 2)))
    figure.update_layout(
        title={"text": f"{method_label} · {rank_label}", "x": 0.02,
               "xanchor": "left", "font": {"size": 18}},
        legend={"orientation": "v", "y": 0.98, "yanchor": "top",
                "x": 0.99, "xanchor": "right", "bgcolor": "rgba(255,255,255,.94)",
                "bordercolor": "#cfd1d4", "borderwidth": 1},
        annotations=[*figure.layout.annotations,
                     dict(text=f"Total removed RMS: {removed_rms:.4g}", x=0.01, y=0.29,
                          xref="paper", yref="paper", showarrow=False)],
    )
    figure.update_yaxes(title_text="Absorbance", row=1, col=1)
    figure.update_yaxes(title_text="Removed", row=2, col=1)
    figure.update_xaxes(title_text="Wavelength (nm)", row=2, col=1)
    return _style(figure, height=520)


def _empty(go, message):
    figure = go.Figure(); figure.add_annotation(text=message, showarrow=False)
    return _style(figure, height=430)


def _style(figure, height):
    figure.update_layout(template="plotly_white", height=height, margin=dict(l=55, r=25, t=65, b=50), hovermode="closest")
    return figure


def run_server(host="127.0.0.1", port=8051, open_browser=True):
    app = create_app()
    if open_browser:
        Timer(1, lambda: webbrowser.open(f"http://{host}:{port}")).start()
    app.run(host=host, port=port, debug=False)


def main(argv=None):
    parser = ArgumentParser(prog="autoqy-smoother-gui")
    parser.add_argument("--host", default="127.0.0.1"); parser.add_argument("--port", default=8051, type=int)
    parser.add_argument("--no-open", action="store_true"); args = parser.parse_args(argv)
    run_server(args.host, args.port, not args.no_open)


if __name__ == "__main__":
    main()
