"""Local browser interface for standalone optical-power processing."""

from argparse import ArgumentParser
from base64 import b64decode
from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
import re
from threading import Timer
from uuid import uuid4
import webbrowser

import numpy as np

from .output import format_value_uncertainty
from .power import (REGION_NAMES, PowerTraceResult, baseline_power_trace,
                    combine_power_traces, load_thorlabs_opm_text,
                    process_power_trace)


REGION_LABELS = (
    "Background", "Open beam", "Background",
    "Background", "Cuvette", "Background",
)
REGION_COLOURS = ("#b85c4a", "#2d6f8e", "#b85c4a",
                  "#b85c4a", "#2d6f8e", "#b85c4a")


def create_app():
    """Create the Dash application; web dependencies are imported lazily."""
    try:
        from dash import (Dash, Input, Output, State, ctx, dcc, html,
                          no_update)
        from dash.exceptions import PreventUpdate
    except ImportError as error:
        raise RuntimeError(
            "The power GUI requires the optional web dependencies. "
            "Install them with: pip install -e .[power-gui]"
        ) from error

    assets = Path(__file__).with_name("assets")
    app = Dash(__name__, title="AutoQY Power", assets_folder=str(assets))
    app.layout = html.Div(className="app-shell", children=[
        dcc.Store(id="power-store", data={"datasets": []}),
        dcc.Download(id="download-results"),
        html.Header(className="app-header", children=[
            html.Div([
                html.P("AUTOQY CORE", className="eyebrow"),
                html.H1("Power treatment"),
                html.P("Baseline optical-power traces and calculate the power at the cuvette.",
                       className="subtitle"),
            ]),
            html.Div("Local session", className="local-badge"),
        ]),
        html.Main(className="workspace", children=[
            html.Aside(className="control-column", children=[
                html.Section(className="panel", children=[
                    html.P("1 · Data", className="step-label"),
                    html.H2("Measurements"),
                    dcc.Upload(
                        id="upload-power", multiple=True,
                        accept=".csv,text/csv",
                        className="upload-box",
                        children=html.Div([
                            html.Span("Choose Thorlabs CSV files"),
                            html.Small("Up to three measurements"),
                        ]),
                    ),
                    html.Label("Active measurement", htmlFor="dataset-select"),
                    dcc.Dropdown(id="dataset-select", placeholder="Upload a CSV file",
                                 clearable=False),
                    html.Button("Remove active measurement", id="delete-measurement",
                                className="button button-secondary"),
                ]),
                html.Section(className="panel", children=[
                    html.P("2 · Baseline", className="step-label"),
                    html.H2("Correction"),
                    html.P("Drag the twelve vertical boundaries on the raw trace. "
                           "The shaded areas are the six regions used by the calculation.",
                           className="helper-text"),
                    html.Label("Polynomial degree", htmlFor="baseline-degree"),
                    dcc.Input(id="baseline-degree", type="number", min=0, step=1,
                              value=3),
                    html.Button("Apply baseline correction", id="apply-baseline",
                                className="button button-primary"),
                    html.Div(id="baseline-message", className="message"),
                ]),
                html.Section(className="panel", children=[
                    html.P("3 · Result", className="step-label"),
                    html.H2("Power at cuvette"),
                    html.Label("Reported uncertainty", htmlFor="uncertainty-output"),
                    dcc.RadioItems(
                        id="uncertainty-output", value="repeatability_sd",
                        options=[
                            {"label": " Repeatability SD", "value": "repeatability_sd"},
                            {"label": " Standard error", "value": "standard_error"},
                        ],
                        className="radio-group",
                    ),
                    html.Button("Calculate power", id="calculate-power",
                                className="button button-accent"),
                    html.Div(id="action-message", className="message"),
                    html.Button("Export JSON", id="export-results",
                                className="button button-secondary"),
                ]),
            ]),
            html.Div(className="plot-column", children=[
                html.Section(className="result-strip", children=[
                    html.Div(id="active-result", className="result-card"),
                    html.Div(id="average-result", className="result-card result-card-accent"),
                ]),
                html.Section(className="plot-panel", children=[
                    dcc.Tabs(id="plot-tabs", value="raw", children=[
                        dcc.Tab(label="Raw trace", value="raw", children=[
                            dcc.Graph(
                                id="raw-plot", figure=_empty_figure("Upload a power trace"),
                                config={"displaylogo": False, "editable": True,
                                        "edits": {"shapePosition": True}},
                            ),
                        ]),
                        dcc.Tab(label="Open-beam baseline", value="open", children=[
                            dcc.Graph(id="open-baseline-plot",
                                      figure=_empty_figure("Apply baseline correction"),
                                      config={"displaylogo": False}),
                        ]),
                        dcc.Tab(label="Cuvette baseline", value="cuvette", children=[
                            dcc.Graph(id="cuvette-baseline-plot",
                                      figure=_empty_figure("Apply baseline correction"),
                                      config={"displaylogo": False}),
                        ]),
                    ]),
                ]),
                html.Section(className="panel repeat-panel", children=[
                    html.H2("Measurements in this session"),
                    html.Div(id="repeat-results"),
                ]),
            ]),
        ]),
    ])

    @app.callback(
        Output("power-store", "data"), Output("action-message", "children"),
        Input("upload-power", "contents"), Input("raw-plot", "relayoutData"),
        Input("calculate-power", "n_clicks"), Input("delete-measurement", "n_clicks"),
        State("upload-power", "filename"), State("power-store", "data"),
        State("dataset-select", "value"), State("baseline-degree", "value"),
        prevent_initial_call=True,
    )
    def reduce_session(contents, relayout, _calculate, _delete, filenames,
                       store, selected, degree):
        store = deepcopy(store or {"datasets": []})
        trigger = ctx.triggered_id
        try:
            if trigger == "upload-power":
                if not contents:
                    raise PreventUpdate
                contents = contents if isinstance(contents, list) else [contents]
                filenames = filenames if isinstance(filenames, list) else [filenames]
                remaining = 3 - len(store["datasets"])
                if remaining <= 0:
                    return no_update, "Three measurements are already loaded."
                loaded = []
                for content, filename in list(zip(contents, filenames))[:remaining]:
                    values = _decode_upload(content, filename)
                    if len(values) < 12:
                        raise ValueError(f"{filename} contains fewer than 12 samples")
                    loaded.append({
                        "id": str(uuid4()), "name": filename,
                        "values_mw": values.tolist(),
                        "regions": _initial_regions(len(values)), "result": None,
                    })
                store["datasets"].extend(loaded)
                return store, f"Loaded {len(loaded)} measurement(s)."

            dataset = _selected_dataset(store, selected)
            if trigger == "raw-plot":
                positions = _positions_from_relayout(relayout, dataset["regions"])
                if positions is None:
                    raise PreventUpdate
                dataset["regions"] = _positions_to_regions(
                    positions, len(dataset["values_mw"])
                )
                dataset["result"] = None
                return store, "Regions updated; calculate again to refresh the result."

            if trigger == "calculate-power":
                degree = _degree(degree)
                result = process_power_trace(
                    dataset["values_mw"], dataset["regions"], degree, dataset["name"]
                )
                dataset["result"] = asdict(result)
                dataset["baseline_degree"] = degree
                return store, "Power calculated from the selected regions."

            if trigger == "delete-measurement":
                store["datasets"] = [item for item in store["datasets"]
                                     if item["id"] != selected]
                return store, "Measurement removed from this local session."
        except PreventUpdate:
            raise
        except Exception as error:
            return no_update, str(error)
        raise PreventUpdate

    @app.callback(
        Output("dataset-select", "options"), Output("dataset-select", "value"),
        Input("power-store", "data"), State("dataset-select", "value"),
    )
    def update_dataset_selector(store, current):
        datasets = (store or {}).get("datasets", [])
        options = [{"label": item["name"], "value": item["id"]} for item in datasets]
        identifiers = {item["id"] for item in datasets}
        return options, current if current in identifiers else (options[0]["value"] if options else None)

    @app.callback(
        Output("raw-plot", "figure"),
        Input("power-store", "data"), Input("dataset-select", "value"),
    )
    def update_raw_plot(store, selected):
        try:
            return _raw_figure(_selected_dataset(store, selected))
        except ValueError:
            return _empty_figure("Upload a power trace")

    @app.callback(
        Output("open-baseline-plot", "figure"),
        Output("cuvette-baseline-plot", "figure"),
        Output("plot-tabs", "value"), Output("baseline-message", "children"),
        Input("apply-baseline", "n_clicks"),
        State("power-store", "data"), State("dataset-select", "value"),
        State("baseline-degree", "value"), prevent_initial_call=True,
    )
    def apply_baseline(_clicks, store, selected, degree):
        try:
            dataset = _selected_dataset(store, selected)
            detail = baseline_power_trace(
                dataset["values_mw"], dataset["regions"], _degree(degree)
            )
            return (_baseline_figure(detail, "open", dataset["name"]),
                    _baseline_figure(detail, "cuvette", dataset["name"]),
                    "open", "Baseline correction applied.")
        except Exception as error:
            empty = _empty_figure("Check the selected regions")
            return empty, empty, "raw", str(error)

    @app.callback(
        Output("active-result", "children"), Output("average-result", "children"),
        Output("repeat-results", "children"),
        Input("power-store", "data"), Input("dataset-select", "value"),
        Input("uncertainty-output", "value"),
    )
    def show_results(store, selected, uncertainty_output):
        datasets = (store or {}).get("datasets", [])
        selected_item = next((item for item in datasets if item["id"] == selected), None)
        active = _single_result_component(html, selected_item, uncertainty_output)
        completed = [PowerTraceResult(**item["result"]) for item in datasets
                     if item.get("result")]
        if completed:
            combined = combine_power_traces(completed, uncertainty_output)
            value, error = format_value_uncertainty(combined.power_mw,
                                                    combined.power_error_mw)
            average = [html.Span("Session average"), html.Strong(f"{value} ± {error} mW"),
                       html.Small(f"{len(completed)} calculated measurement(s)")]
        else:
            average = [html.Span("Session average"), html.Strong("—"),
                       html.Small("Calculate at least one measurement")]
        rows = [_repeat_row(html, item, uncertainty_output) for item in datasets]
        return active, average, rows or html.P("No measurements loaded.", className="empty-copy")

    @app.callback(
        Output("download-results", "data"), Input("export-results", "n_clicks"),
        State("power-store", "data"), State("uncertainty-output", "value"),
        State("baseline-degree", "value"), prevent_initial_call=True,
    )
    def export_results(_clicks, store, uncertainty_output, degree):
        datasets = (store or {}).get("datasets", [])
        completed = [PowerTraceResult(**item["result"]) for item in datasets
                     if item.get("result")]
        if not completed:
            raise PreventUpdate
        combined = combine_power_traces(completed, uncertainty_output)
        document = {
            "schema_version": 1,
            "format": "thorlabs_opm_csv",
            "baseline_polynomial_degree": _degree(degree),
            "uncertainty_output": uncertainty_output,
            "measurements": [
                {"path": item["name"], "regions": item["regions"]}
                for item in datasets
            ],
            "result": combined.as_dict(),
        }
        return dcc.send_string(json.dumps(document, indent=2) + "\n",
                               "PowerProcessing_results.json")

    return app


def _decode_upload(content, filename):
    try:
        encoded = content.split(",", 1)[1]
        raw = b64decode(encoded)
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            text = raw.decode("latin-1")
    except Exception as error:
        raise ValueError(f"Could not decode {filename}") from error
    return load_thorlabs_opm_text(text, filename)


def _degree(value):
    if isinstance(value, bool) or value is None or int(value) != value or value < 0:
        raise ValueError("Polynomial degree must be a nonnegative integer")
    return int(value)


def _selected_dataset(store, identifier):
    for dataset in (store or {}).get("datasets", []):
        if dataset["id"] == identifier:
            return dataset
    raise ValueError("Select a measurement first")


def _initial_regions(length):
    positions = np.linspace(0, length, 13)[:-1]
    return _positions_to_regions(positions, length)


def _region_positions(regions):
    return [boundary for name in REGION_NAMES for boundary in regions[name]]


def _positions_to_regions(positions, length):
    positions = sorted(float(position) for position in positions)
    if len(positions) != 12:
        raise ValueError("Exactly twelve region boundaries are required")
    result = {}
    for index, name in enumerate(REGION_NAMES):
        start = max(0, min(length - 1, int(round(positions[index * 2]))))
        end = max(start + 1, min(length, int(round(positions[index * 2 + 1]))))
        result[name] = [start, end]
    return result


def _positions_from_relayout(relayout, regions):
    if not relayout:
        return None
    positions = _region_positions(regions)
    if "shapes" in relayout:
        boundaries = [shape["x0"] for shape in relayout["shapes"]
                      if str(shape.get("name", "")).startswith("boundary-")]
        return boundaries if len(boundaries) == 12 else None
    changed = False
    for key, value in relayout.items():
        match = re.fullmatch(r"shapes\[(\d+)]\.x[01]", key)
        if match and int(match.group(1)) < 12:
            positions[int(match.group(1))] = value
            changed = True
    return positions if changed else None


def _empty_figure(message):
    try:
        import plotly.graph_objects as go
    except ImportError:
        return {}
    figure = go.Figure()
    figure.add_annotation(text=message, x=0.5, y=0.5, xref="paper", yref="paper",
                          showarrow=False, font={"color": "#6c7280", "size": 16})
    return _style_figure(figure)


def _raw_figure(dataset):
    import plotly.graph_objects as go

    values = np.asarray(dataset["values_mw"])
    x = np.arange(len(values))
    figure = go.Figure(go.Scatter(x=x, y=values, mode="lines", name="Power",
                                  line={"color": "#20242c", "width": 1.5}))
    positions = _region_positions(dataset["regions"])
    for index, (position, colour) in enumerate(zip(positions,
                                                   np.repeat(REGION_COLOURS, 2))):
        figure.add_shape(
            type="line", x0=position, x1=position, y0=0, y1=1,
            yref="paper", line={"color": colour, "width": 2}, opacity=0.65,
            editable=True, name=f"boundary-{index}",
        )
    for name, label, colour in zip(REGION_NAMES, REGION_LABELS, REGION_COLOURS):
        start, end = dataset["regions"][name]
        figure.add_vrect(x0=start, x1=end, fillcolor=colour, opacity=0.10,
                         line_width=0, layer="below", editable=False, name=f"region-{name}")
        figure.add_annotation(x=(start + end) / 2, y=1.04, xref="x", yref="paper",
                              text=label, showarrow=False,
                              font={"size": 10, "color": colour})
    figure.update_layout(title=dataset["name"], dragmode="pan")
    figure.update_xaxes(title="Sample index", range=[0, max(len(values) - 1, 1)])
    figure.update_yaxes(title="Power (mW)")
    return _style_figure(figure)


def _baseline_figure(detail, case, filename):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    if case == "open":
        baseline = detail.open_baseline_mw
        corrected = detail.open_corrected_mw
        first, signal, last = (detail.regions[name] for name in REGION_NAMES[:3])
        title = "Open beam · no jacket and no cuvette"
    else:
        baseline = detail.cuvette_baseline_mw
        corrected = detail.cuvette_corrected_mw
        first, signal, last = (detail.regions[name] for name in REGION_NAMES[3:])
        title = "With jacket and cuvette with solvent"
    start, end = first[0], last[1]
    x = np.arange(len(detail.values_mw))[start:end]
    figure = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.14,
                           subplot_titles=("Original trace and fitted baseline",
                                           "Baseline-corrected trace"))
    figure.add_trace(go.Scatter(x=x, y=detail.values_mw[start:end], mode="lines",
                                name="Power", line={"color": "#20242c"}), row=1, col=1)
    figure.add_trace(go.Scatter(x=x, y=baseline[start:end], mode="lines",
                                name="Baseline", line={"color": "#b85c4a", "dash": "dash"}),
                     row=1, col=1)
    figure.add_trace(go.Scatter(x=x, y=corrected[start:end], mode="lines",
                                name="Corrected power", line={"color": "#d67b36"}),
                     row=2, col=1)
    figure.add_hline(y=0, line={"color": "#8a8f99", "width": 1}, row=2, col=1)
    figure.add_vrect(x0=signal[0], x1=signal[1], fillcolor="#2d6f8e", opacity=0.10,
                     line_width=0, row="all", col=1)
    figure.update_layout(title=f"{filename}<br><sup>{title}</sup>", height=650)
    figure.update_xaxes(title="Sample index", row=2, col=1)
    figure.update_yaxes(title="Power (mW)", row=1, col=1)
    figure.update_yaxes(title="Power (mW)", row=2, col=1)
    return _style_figure(figure)


def _style_figure(figure):
    figure.update_layout(
        template="plotly_white", margin={"l": 70, "r": 30, "t": 80, "b": 60},
        paper_bgcolor="#ffffff", plot_bgcolor="#fbfaf7",
        font={"family": "Inter, Segoe UI, sans-serif", "color": "#20242c"},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02,
                "xanchor": "right", "x": 1},
    )
    figure.update_xaxes(showgrid=False, zeroline=False)
    figure.update_yaxes(gridcolor="#e6e2db", zeroline=False)
    return figure


def _single_result_component(html, dataset, uncertainty_output):
    if not dataset or not dataset.get("result"):
        return [html.Span("Active measurement"), html.Strong("—"),
                html.Small("Adjust regions, then calculate")]
    trace = PowerTraceResult(**dataset["result"])
    error = (trace.repeatability_sd_mw if uncertainty_output == "repeatability_sd"
             else trace.standard_error_mw)
    value, formatted_error = format_value_uncertainty(trace.power_mw, error)
    return [html.Span(dataset["name"]), html.Strong(f"{value} ± {formatted_error} mW"),
            html.Small(f"Open {trace.open_beam_power_mw:.4g} · Cuvette "
                       f"{trace.cuvette_power_mw:.4g} mW")]


def _repeat_row(html, dataset, uncertainty_output):
    if dataset.get("result"):
        trace = PowerTraceResult(**dataset["result"])
        error = (trace.repeatability_sd_mw if uncertainty_output == "repeatability_sd"
                 else trace.standard_error_mw)
        value, formatted_error = format_value_uncertainty(trace.power_mw, error)
        result = f"{value} ± {formatted_error} mW"
        state = "Calculated"
    else:
        result, state = "—", "Not calculated"
    return html.Div(className="repeat-row", children=[
        html.Div([html.Strong(dataset["name"]), html.Small(state)]), html.Span(result)
    ])


def run_server(host="127.0.0.1", port=8050, open_browser=True):
    """Run the local application and optionally open it in the default browser."""
    app = create_app()
    if open_browser:
        Timer(1, lambda: webbrowser.open(f"http://{host}:{port}")).start()
    app.run(host=host, port=port, debug=False)


def main(argv=None):
    parser = ArgumentParser(prog="autoqy-power-gui")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8050, type=int)
    parser.add_argument("--no-open", action="store_true")
    args = parser.parse_args(argv)
    run_server(args.host, args.port, not args.no_open)


if __name__ == "__main__":
    main()
