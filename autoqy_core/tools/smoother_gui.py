"""Browser GUI for treating spectra and calculating molar absorptivity."""

from argparse import ArgumentParser
import base64
from pathlib import Path
import subprocess
from threading import Lock
import time

import numpy as np

from ..epsilon import (EpsilonResult, NMRSubtractionResult,
                       calculate_epsilon_statistics, export_epsilon_csv,
                       export_nmr_subtraction_csv, load_epsilon_table,
                       nonnegative_error_bounds, reconstruct_product_from_nmr)
from ..gui_window import serve_gui
from ..plot_style import (ANALYSIS_TRACE_PALETTE, PLOT_BLUE, PLOT_ORANGE,
                          PLOT_PURPLE)
from ..smoother import (SpectralDataset, analyze_svd, baseline_spectra,
                        export_smoothed_text, load_spectral_bytes,
                        savgol_window_points, select_wavelengths,
                        smooth_reconstruction)
from ..version import get_project_version


def create_app():
    try:
        from dash import ALL, Dash, Input, Output, State, ctx, dcc, html, no_update
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as error:
        raise RuntimeError("The spectral treatment GUI requires the 'gui' optional dependencies") from error

    assets = Path(__file__).parents[1] / "assets"
    project_version = get_project_version()
    app = Dash(__name__, assets_folder=str(assets), suppress_callback_exceptions=True)
    app.title = f"AutoQY Spectral treatment · {project_version}"
    window_state = {"close_requested_at": None, "last_heartbeat_at": None}
    window_lock = Lock()
    app.server.config["AUTOQY_WINDOW_STATE"] = (window_state, window_lock)

    @app.server.post("/_autoqy_heartbeat")
    def autoqy_heartbeat():
        with window_lock:
            window_state["close_requested_at"] = None
            window_state["last_heartbeat_at"] = time.monotonic()
        return "", 204

    @app.server.post("/_autoqy_window_closed")
    def autoqy_window_closed():
        with window_lock:
            window_state["close_requested_at"] = time.monotonic()
        return "", 204

    app.index_string = """<!DOCTYPE html>
<html><head>{%metas%}<title>{%title%}</title>{%favicon%}{%css%}</head>
<body>{%app_entry%}<footer>{%config%}{%scripts%}{%renderer%}</footer>
<script>
(() => {
  let heartbeatFailures = 0;
  let closingForStoppedServer = false;
  const heartbeat = async () => {
    try {
      const response = await fetch('/_autoqy_heartbeat', {
        method: 'POST', keepalive: true, cache: 'no-store'
      });
      if (!response.ok) throw new Error('heartbeat failed');
      heartbeatFailures = 0;
    } catch (_) {
      heartbeatFailures += 1;
      if (heartbeatFailures >= 3 && !closingForStoppedServer) {
        closingForStoppedServer = true;
        window.close();
        window.setTimeout(() => {
          document.body.innerHTML = '<main style="font-family:Segoe UI,sans-serif;padding:2rem">'
            + '<h1>AutoQY GUI stopped</h1><p>The terminal has been closed. Close this webpage and restart the AutoQY GUI to continue.</p></main>';
        }, 150);
      }
    }
  };
  heartbeat();
  const heartbeatTimer = window.setInterval(heartbeat, 1000);
  window.addEventListener('pageshow', heartbeat);
  window.addEventListener('pagehide', () => {
    window.clearInterval(heartbeatTimer);
    navigator.sendBeacon('/_autoqy_window_closed');
  });
  document.addEventListener('click', (event) => {
    const popup = event.target.closest('.info-popup');
    if (popup) {
      event.preventDefault();
      event.stopPropagation();
      popup.focus();
    } else if (document.activeElement?.classList.contains('info-popup')) {
      document.activeElement.blur();
    }
  }, true);
})();
window.autoqyDownloadPlot = (graphId, figure, format, includeHeader) => {
  const graph = document.querySelector(`#${graphId} .js-plotly-plot`);
  if (!graph || !window.Plotly) return 'Plot export is unavailable.';
  const title = typeof figure?.layout?.title === 'string'
    ? figure.layout.title : figure?.layout?.title?.text;
  const filename = String(title || 'autoqy-spectral-treatment')
    .replace(/<[^>]+>/g, '')
    .normalize('NFKD')
    .replace(/[^a-zA-Z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '')
    .toLowerCase() || 'autoqy-spectral-treatment';
  const width = Math.max(900, Math.round(graph._fullLayout?.width || 1200));
  let height = Math.max(600, Math.round(graph._fullLayout?.height || 800));
  const keepHeader = Array.isArray(includeHeader) && includeHeader.includes('on');
  const exportData = JSON.parse(JSON.stringify(figure?.data || []));
  const exportLayout = JSON.parse(JSON.stringify(figure?.layout || {}));
  if (keepHeader) {
    exportLayout.showlegend = true;
  } else {
    exportLayout.title = null;
    exportLayout.showlegend = false;
    exportLayout.margin = Object.assign({}, exportLayout.margin || {}, {t: 58});
    height = Math.max(480, height - 127);
  }
  exportLayout.width = width;
  exportLayout.height = height;

  const extension = `.${format}`;
  const suggestedName = `${filename}${extension}`;
  const mimeType = format === 'svg' ? 'image/svg+xml' : 'image/png';
  let fileHandlePromise = null;
  if (window.showSaveFilePicker) {
    try {
      fileHandlePromise = window.showSaveFilePicker({
        suggestedName,
        types: [{
          description: `${format.toUpperCase()} image`,
          accept: {[mimeType]: [extension]}
        }]
      });
    } catch (_) {
      return `Could not open the ${format.toUpperCase()} save dialog.`;
    }
  }

  const holder = document.createElement('div');
  holder.style.cssText = `position:fixed;left:-100000px;top:0;width:${width}px;height:${height}px`;
  document.body.appendChild(holder);
  const renderImage = () => window.Plotly.newPlot(
    holder, exportData, exportLayout, {staticPlot: true, displayModeBar: false}
  ).then(() => window.Plotly.toImage(holder, {
    format, width, height, scale: format === 'png' ? 2 : 1
  })).finally(() => {
    window.Plotly.purge(holder);
    holder.remove();
  });

  const saveImage = (fileHandle) => renderImage().then((dataUrl) => {
    if (!fileHandle) {
      const link = document.createElement('a');
      link.href = dataUrl;
      link.download = suggestedName;
      link.click();
      return `${format.toUpperCase()} download started.`;
    }
    return window.fetch(dataUrl)
      .then((response) => response.blob())
      .then((blob) => fileHandle.createWritable()
        .then((writer) => writer.write(blob).then(() => writer.close())))
      .then(() => `${format.toUpperCase()} saved.`);
  });

  return (fileHandlePromise || Promise.resolve(null))
    .then(saveImage)
    .catch((error) => error?.name === 'AbortError'
      ? 'Save cancelled.'
      : `Could not create the ${format.toUpperCase()} image.`);
};
</script></body></html>"""

    def info_popup(text):
        return html.Span(className="info-popup", tabIndex=0, children=[
            html.Span("i", className="info-popup-icon", **{"aria-hidden": "true"}),
            html.Span(text, className="info-popup-content"),
        ])

    app.layout = html.Div(className="app-shell", children=[
        html.Header(className="app-header", children=[html.Div([
            html.P("AUTOQY", className="eyebrow"),
            html.H1("Spectral treatment"),
            html.P("Prepare wavelength-resolved data and calculate molar absorptivity.",
                   className="subtitle"),
        ]), html.Div(f"Version {project_version}", className="local-badge")]),
        html.Main(className="workspace epsilon-workspace", children=[
            html.Aside(className="control-column smoother-controls", children=[
                html.Details(open=True, className="panel tool-details", children=[
                    html.Summary([html.Span("1 · Data", className="step-label"),
                                  "Spectral data", info_popup(
                                      "Load one or more spectra. File types are detected automatically; "
                                      "folder loading also makes the source folder the default export location."
                                  )]),
                    dcc.Upload(
                        id="upload-spectra", className="upload-box", multiple=True,
                        children=html.Div([
                            html.Span("Drop or choose one or more spectral files"),
                            html.Small("SpectraGryph .dat, Avantes .Abs8, TSV, or CSV"),
                        ]),
                    ),
                    html.Button("Open files from folder", id="open-local-spectra",
                                className="button button-secondary"),
                    html.Button("Clear all spectra", id="clear-dataset",
                                className="button button-secondary", disabled=True),
                    dcc.Loading(
                        html.Div(id="load-message", className="message"),
                        type="circle",
                    ),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("2 · Range", className="step-label"),
                                  "Wavelengths", info_popup(
                                      "The selected range is used for the preview, molar-absorptivity "
                                      "calculation, uncertainty statistics, and exported CSV."
                                  )]),
                    html.Div([
                        dcc.Input(id="wavelength-low", type="number", placeholder="Start (nm)", disabled=True),
                        dcc.Input(id="wavelength-high", type="number", placeholder="End (nm)", disabled=True),
                    ], className="input-row"),
                    html.Details(open=False, className="nested-tool", children=[
                        html.Summary(["Preprocess spectra", info_popup(
                            "Baseline and Savitzky–Golay operate on each spectrum independently. "
                            "SVD mixes columns and should only be used for ordered time-series data, "
                            "never for independent concentration replicates. Uploaded values are unchanged."
                        )]),
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
                            html.Summary(["Smoothing parameters", info_popup(
                                "The Savitzky–Golay window is entered in nanometres and converted "
                                "to a valid odd number of detector points; polynomial order must be lower."
                            )]),
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
                        html.Div(id="smoothing-message", className="message"),
                    ]),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("3 · Beer–Lambert", className="step-label"),
                                  "Concentrations", info_popup(
                                      "Enter each measured solution concentration directly in mol/L "
                                      "and its optical path length in centimetres."
                                  )]),
                    html.Div(id="concentration-parameters"),
                    html.Div(id="concentration-message", className="message"),
                ]),
                html.Details(open=False, className="panel tool-details export-panel", children=[
                    html.Summary([html.Span("4 · Output", className="step-label"),
                                  "Export processed dataset", info_popup(
                                      "Processed absorbance can be saved without concentrations. "
                                      "When every concentration and path length is provided, the CSV "
                                      "also contains individual ε spectra and their statistics."
                                  )]),
                    html.Label("Save folder"),
                    html.Div(className="input-row", children=[
                        dcc.Input(id="save-folder", type="text",
                                  placeholder="Choose a save folder"),
                        html.Button("Choose folder", id="choose-save-folder",
                                    className="button button-secondary"),
                    ]),
                    html.Label("CSV file name (optional)"),
                    dcc.Input(id="save-filename", type="text",
                              placeholder="Default name"),
                    html.Small("Leave blank to use processed-absorbance.csv or "
                               "epsilon-spectra-reactant.csv."),
                    html.Div(className="export-mode-control", children=[
                        dcc.Checklist(
                            id="export-nonnegative", value=[], className="toggle-control",
                            options=[{
                                "label": "Convert negative absorbance and ε values to 0",
                                "value": "on",
                            }],
                        ),
                        html.Small("Default: negative values are preserved. Enable this "
                                   "option to clamp them to zero in the saved CSV."),
                    ]),
                    html.Button("Save processed CSV", id="export-epsilon",
                                className="button button-primary", disabled=True),
                    html.Div(id="epsilon-save-message", className="message"),
                    dcc.ConfirmDialog(id="confirm-epsilon-overwrite"),
                ]),
                html.Details(open=False, className="panel tool-details", children=[
                    html.Summary([html.Span("5 · Optional", className="step-label"),
                                  "NMR-guided PSS subtraction", info_popup(
                                      "Load one UV–Vis dataset containing reactant and final PSS. "
                                      "After shared normalization, product ε is reconstructed as "
                                      "(PSS − x · reactant) / (1 − x), where x is the reactant fraction at PSS."
                                  )]),
                    dcc.Upload(
                        id="nmr-upload", className="upload-box compact-upload", multiple=False,
                        children=html.Div([
                            html.Span("Drop one dataset containing reactant and PSS"),
                            html.Small("First spectrum = reactant · last spectrum = PSS"),
                        ]),
                    ),
                    html.Details(open=False, className="nested-tool", children=[
                        html.Summary(["Preprocess reactant and PSS", info_popup(
                            "Apply the same baseline interval and Savitzky–Golay settings to the "
                            "reactant and PSS spectra before their normalized subtraction."
                        )]),
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
                    html.Label("CSV file names (optional)"),
                    html.Div(className="input-row", children=[
                        dcc.Input(id="nmr-reactant-filename", type="text",
                                  placeholder="Reactant file name"),
                        dcc.Input(id="nmr-product-filename", type="text",
                                  placeholder="Product file name"),
                    ]),
                    html.Small("Leave blank to use the standard reactant and product names."),
                    html.Button("Save reactant + NMR-derived ε CSVs", id="export-nmr",
                                className="button button-accent", disabled=True),
                    html.Div(id="nmr-save-message", className="message"),
                    dcc.ConfirmDialog(id="confirm-nmr-overwrite"),
                ]),
            ]),
            html.Div(className="plot-column", children=[
                html.Section(className="plot-panel", children=[
                    html.Div(className="plot-toolbar", children=[
                        html.Div(className="plot-download-actions", children=[
                            html.Button("Save PNG", id="save-epsilon-png",
                                        className="button button-secondary compact-button"),
                            html.Button("Save SVG", id="save-epsilon-svg",
                                        className="button button-secondary compact-button"),
                        ]),
                        html.Div(className="plot-toolbar-controls", children=[
                            html.Div(className="plot-quick-options", children=[
                                dcc.Checklist(
                                    id="minimal-spectrum-colors", value=[],
                                    className="toggle-control plot-option-toggle",
                                    options=[{"label": "Minimal colors", "value": "on"}],
                                ),
                                dcc.Checklist(
                                    id="include-plot-header", value=[],
                                    className="toggle-control plot-option-toggle",
                                    options=[{"label": "Save title + legend", "value": "on"}],
                                ),
                            ]),
                            html.Details(className="plot-options", children=[
                                html.Summary("Legend options"),
                                dcc.Checklist(
                                    id="show-spectrum-legend", value=["on"],
                                    className="toggle-control plot-legend-toggle",
                                    options=[{"label": "Show legend", "value": "on"}],
                                ),
                                html.Label("Legend labels (one per spectrum)"),
                                dcc.Textarea(
                                    id="spectrum-legend-labels", value="",
                                    className="legend-editor",
                                    placeholder="Load spectra to edit their legend labels",
                                ),
                                html.Small(
                                    "Use one line per spectrum. These names affect only the plot; "
                                    "CSV labels and source filenames stay unchanged."
                                ),
                            ]),
                        ]),
                    ]),
                    html.Div(id="epsilon-image-message", className="image-export-message"),
                    dcc.Graph(
                        id="epsilon-plot", figure=_empty(go, "Load spectral data to begin"),
                        config=_graph_config(),
                    ),
                ]),
                html.Section(id="nmr-plot-panel", className="plot-panel", style={"display": "none"},
                             children=[
                                 html.Div(className="plot-toolbar nmr-plot-toolbar", children=[
                                     html.Strong("NMR-guided PSS subtraction"),
                                     html.Div(className="plot-download-actions", children=[
                                         html.Button("Save PNG", id="save-nmr-png",
                                                     className="button button-secondary compact-button"),
                                         html.Button("Save SVG", id="save-nmr-svg",
                                                     className="button button-secondary compact-button"),
                                     ]),
                                 ]),
                                 html.Div(id="nmr-image-message", className="image-export-message"),
                                 dcc.Graph(
                                     id="nmr-plot", figure=_empty(go, "Load NMR spectra"),
                                     config=_graph_config(),
                                 ),
                             ]),
                html.Section(className="panel result-summary", children=[
                    html.Div(className="section-title-row", children=[
                        html.H2("Result"), info_popup(
                            "Reports whether ε can be calculated or exported and highlights "
                            "missing concentration, preprocessing, or non-negativity requirements."
                        ),
                    ]),
                    html.Div(id="result-message", className="helper-text"),
                ]),
                html.Details(open=False, className="panel tool-details error-panel", children=[
                    html.Summary(["Python errors", info_popup(
                        "This panel contains detailed parser, preprocessing, reconstruction, "
                        "or export exceptions. An empty panel means no Python error is active."
                    )]),
                    html.Pre(id="load-error"), html.Pre(id="preview-error"),
                    html.Pre(id="svd-error"), html.Pre(id="nmr-error"),
                    html.Pre(id="export-error"),
                ]),
            ]),
        ]),
        dcc.Store(id="dataset-store"),
        dcc.Store(id="epsilon-store"),
        dcc.Store(id="processed-store"),
        dcc.Store(id="nmr-spectra-store"),
        dcc.Store(id="nmr-result-store"),
        dcc.Store(id="source-folder-store"),
    ])

    def register_image_download(graph_id, png_button, svg_button, message_id):
        app.clientside_callback(
            """
            function(pngClicks, svgClicks, figure, includeHeader) {
              const context = window.dash_clientside.callback_context;
              if (!context.triggered || !context.triggered.length) {
                return window.dash_clientside.no_update;
              }
              const trigger = context.triggered[0].prop_id.split('.')[0];
              const format = trigger.endsWith('svg') ? 'svg' : 'png';
              return window.autoqyDownloadPlot(GRAPH_ID, figure, format, includeHeader);
            }
            """.replace("GRAPH_ID", repr(graph_id)),
            Output(message_id, "children"),
            Input(png_button, "n_clicks"),
            Input(svg_button, "n_clicks"),
            State(graph_id, "figure"),
            State("include-plot-header", "value"),
            prevent_initial_call=True,
        )

    register_image_download(
        "epsilon-plot", "save-epsilon-png", "save-epsilon-svg",
        "epsilon-image-message",
    )
    register_image_download(
        "nmr-plot", "save-nmr-png", "save-nmr-svg", "nmr-image-message",
    )

    @app.callback(
        Output("dataset-store", "data"), Output("load-message", "children"),
        Output("upload-spectra", "contents"), Output("upload-spectra", "filename"),
        Output("load-error", "children"), Output("source-folder-store", "data"),
        Input("upload-spectra", "contents"), Input("clear-dataset", "n_clicks"),
        Input("open-local-spectra", "n_clicks"),
        State("upload-spectra", "filename"),
        prevent_initial_call=True,
    )
    def load(contents, _, __, filenames):
        if ctx.triggered_id == "clear-dataset":
            return None, "All spectra cleared.", None, None, "", None
        try:
            if ctx.triggered_id == "open-local-spectra":
                paths = _choose_files()
                if not paths:
                    return no_update, "File selection cancelled.", no_update, no_update, "", no_update
                packed, message = _load_local_paths(paths)
                source_folder = str(Path(paths[0]).resolve().parent)
                return packed, message, None, None, "", source_folder
            if not contents:
                return None, "Choose one or more files to begin.", no_update, no_update, "", None
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
                    "concentrations, and path lengths from an AutoQY ε CSV or legacy TSV.",
                    no_update, no_update, "", None,
                )
            loaded = []
            for content, filename in zip(contents, filenames):
                payload = base64.b64decode(content.split(",", 1)[1])
                selected_format = ("avantes_abs8" if Path(filename or "").suffix.lower() == ".abs8"
                                   else "auto")
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
                no_update, no_update, "", None,
            )
        except Exception as error:
            return (None, "Could not load the selected files.", no_update, no_update,
                    f"Load error: {type(error).__name__}: {error}", no_update)

    @app.callback(
        Output("save-folder", "value"),
        Input("source-folder-store", "data"), Input("choose-save-folder", "n_clicks"),
        State("save-folder", "value"), prevent_initial_call=True,
    )
    def select_save_folder(source_folder, _, current_folder):
        if ctx.triggered_id == "source-folder-store":
            return source_folder
        selected = _choose_folder(current_folder or source_folder)
        return selected or no_update

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
        Output("spectrum-legend-labels", "value"),
        Input("dataset-store", "data"),
    )
    def populate_legend_labels(data):
        return "\n".join(data["labels"]) if data else ""

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
        Output("epsilon-store", "data"), Output("processed-store", "data"),
        Output("export-epsilon", "disabled"), Output("preview-error", "children"),
        Input("dataset-store", "data"), Input("wavelength-low", "value"),
        Input("wavelength-high", "value"), Input("baseline-enabled", "value"),
        Input("baseline-low", "value"), Input("baseline-high", "value"),
        Input("smoothing-method", "value"), Input("savgol-window", "value"),
        Input("savgol-order", "value"), Input("svd-enabled", "value"),
        Input("svd-rank", "value"),
        Input({"type": "direct-concentration", "index": ALL}, "value"),
        Input({"type": "path-length", "index": ALL}, "value"),
        Input("show-spectrum-legend", "value"),
        Input("spectrum-legend-labels", "value"),
        Input("minimal-spectrum-colors", "value"),
    )
    def preview(data, wavelength_low, wavelength_high, baseline_enabled,
                baseline_low, baseline_high, method, sg_width, sg_order,
                svd_enabled, svd_rank, concentrations, path_lengths,
                show_legend, legend_text, minimal_colors):
        if not data:
            return (_empty(go, "Load spectral data to begin"),
                    "No result yet.", "", "Smoothing is off.", None, None, True, "")
        try:
            dataset, original, processed, smoothing_message = _prepare_processing(
                data, wavelength_low, wavelength_high, baseline_enabled,
                baseline_low, baseline_high, method, sg_width, sg_order,
                svd_enabled, svd_rank,
            )
            concentration_data = _read_concentrations(
                len(data["labels"]), concentrations, path_lengths
            )
            processed_data = _pack(
                SpectralDataset(
                    dataset.wavelengths, dataset.coordinates, processed,
                    dataset.source_format, dataset.interpolated_values,
                ),
                data["labels"], data.get("filenames", []),
            )
            plot_labels = _legend_labels(data["labels"], legend_text)
            legend_visible = "on" in (show_legend or [])
            use_minimal_colors = "on" in (minimal_colors or [])
            if concentration_data is None:
                message = ("Enter the concentration and path length for every "
                           "spectrum to calculate molar absorptivity.")
                return (_absorbance_figure(
                            go, dataset, original, processed, plot_labels,
                            method, svd_enabled, svd_rank, legend_visible,
                            use_minimal_colors,
                        ),
                        "Processed absorbance is ready to export; ε is waiting for "
                        "concentration inputs.",
                        message, smoothing_message, None, processed_data, False, "")
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
                result_message = ("One ε spectrum calculated. At least two independent "
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
                _epsilon_figure(
                    go, make_subplots, dataset, original, result,
                    plot_labels, method, svd_enabled, svd_rank, legend_visible,
                    use_minimal_colors,
                ),
                result_message, concentration_message, smoothing_message,
                _pack_epsilon(result, data["labels"]), processed_data, False, "",
            )
        except Exception as error:
            return (_empty(go, "Preview unavailable; open Python errors below."),
                    "No valid result.", "", "", None, None, True,
                    f"Preview error: {type(error).__name__}: {error}")

    @app.callback(
        Output("epsilon-save-message", "children"),
        Output("confirm-epsilon-overwrite", "displayed"),
        Output("confirm-epsilon-overwrite", "message"),
        Output("export-error", "children"),
        Input("export-epsilon", "n_clicks"),
        Input("confirm-epsilon-overwrite", "submit_n_clicks"),
        State("epsilon-store", "data"), State("processed-store", "data"),
        State("save-folder", "value"), State("save-filename", "value"),
        State("export-nonnegative", "value"),
        prevent_initial_call=True,
    )
    def save_epsilon(_, __, epsilon_data, processed_data, folder, requested_filename,
                     nonnegative):
        try:
            filename, csv_text = _export_csv_payload(
                epsilon_data, processed_data, "on" in (nonnegative or [])
            )
            destination = _save_path(
                folder, _csv_filename(requested_filename, filename)
            )
            if ctx.triggered_id == "export-epsilon" and destination.exists():
                return (f"Existing file: {destination}", True,
                        f"Overwrite {destination.name}?", "")
            destination.write_text(csv_text, encoding="utf-8")
            return f"Saved {destination}", False, "", ""
        except Exception as error:
            return ("Save failed.", False, "",
                    f"Export error: {type(error).__name__}: {error}")

    @app.callback(
        Output("nmr-spectra-store", "data"), Output("nmr-load-message", "children"),
        Output("nmr-upload", "contents"), Output("nmr-upload", "filename"),
        Output("clear-nmr", "disabled"), Output("nmr-error", "children"),
        Input("nmr-upload", "contents"), Input("clear-nmr", "n_clicks"),
        State("nmr-upload", "filename"),
        prevent_initial_call=True,
    )
    def load_nmr(contents, _, filenames):
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
            selected_format = ("avantes_abs8" if Path(filenames or "").suffix.lower() == ".abs8"
                               else "auto")
            dataset = load_spectral_bytes(payload, selected_format)
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
            return (_empty(go, "Calculate or reload reactant ε first"), {"display": "block"},
                    "Reactant ε is required before NMR subtraction.", None, True,
                    "status-message status-warning", "Preprocessing is waiting for reactant ε.")
        try:
            epsilon_result, _ = _unpack_epsilon(epsilon_data)
            low = max(float(np.min(nmr_data["reactant_wavelengths"])),
                      float(np.min(nmr_data["pss_wavelengths"])))
            high = min(float(np.max(nmr_data["reactant_wavelengths"])),
                       float(np.max(nmr_data["pss_wavelengths"])))
            mask = ((epsilon_result.wavelengths >= low) &
                    (epsilon_result.wavelengths <= high))
            if np.count_nonzero(mask) < 2:
                raise ValueError("NMR and reactant ε datasets do not share a wavelength range")
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
                message = ("Product ε reconstructed from normalized "
                           "(PSS − x · reactant) / (1 − x). The band combines the "
                           "reactant ε scale SD and selected NMR error through minimum and maximum values.")
                status = "status-message status-ok"
            return (figure, {"display": "block"}, message, _pack_nmr(result), False,
                    status, processing_message)
        except Exception as error:
            return (_empty(go, "NMR reconstruction unavailable"), {"display": "block"},
                    f"NMR reconstruction error: {type(error).__name__}: {error}", None, True,
                    "status-message status-stop", "Preprocessing could not be applied.")

    @app.callback(
        Output("nmr-save-message", "children"),
        Output("confirm-nmr-overwrite", "displayed"),
        Output("confirm-nmr-overwrite", "message"),
        Input("export-nmr", "n_clicks"),
        Input("confirm-nmr-overwrite", "submit_n_clicks"),
        State("nmr-result-store", "data"), State("epsilon-store", "data"),
        State("nmr-export-raw", "value"), State("save-folder", "value"),
        State("nmr-reactant-filename", "value"),
        State("nmr-product-filename", "value"),
        prevent_initial_call=True,
    )
    def save_nmr(_, __, nmr_data, epsilon_data, preserve_negative, folder,
                 requested_reactant_filename, requested_product_filename):
        try:
            if not nmr_data or not epsilon_data:
                raise ValueError("Calculate reactant and NMR-derived ε before saving")
            reactant_path = _save_path(folder, _csv_filename(
                requested_reactant_filename, "epsilon-spectra-reactant.csv"
            ))
            product_path = _save_path(folder, _csv_filename(
                requested_product_filename, "epsilon-spectra-product.csv"
            ))
            if reactant_path == product_path:
                raise ValueError("Reactant and product file names must be different")
            existing = [path.name for path in (reactant_path, product_path) if path.exists()]
            if ctx.triggered_id == "export-nmr" and existing:
                return ("Existing file(s): " + ", ".join(existing), True,
                        "Overwrite " + " and ".join(existing) + "?",)
            epsilon_result, labels = _unpack_epsilon(epsilon_data)
            reactant_text = export_epsilon_csv(epsilon_result, labels)
            product_text = export_nmr_subtraction_csv(
                _unpack_nmr(nmr_data), "on" in (preserve_negative or [])
            )
            reactant_path.write_text(reactant_text, encoding="utf-8")
            product_path.write_text(product_text, encoding="utf-8")
            return f"Saved {reactant_path.name} and {product_path.name} in {reactant_path.parent}", False, ""
        except Exception as error:
            return f"Save error: {type(error).__name__}: {error}", False, ""

    return app


def _parameter_cards(html, dcc, labels, concentrations=None, path_lengths=None):
    cards = []
    for index, label in enumerate(labels):
        cards.append(html.Div(className="spectrum-card", children=[
            html.H3(label),
            html.Div(className="parameter-grid", children=[
                html.Div([html.Label("Concentration (mol/L)"),
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


def _load_local_paths(paths):
    paths = [Path(path) for path in paths]
    if len(paths) == 1 and paths[0].suffix.lower() in {".csv", ".tsv"}:
        try:
            result, labels = load_epsilon_table(paths[0].read_text(encoding="utf-8-sig"))
            dataset = SpectralDataset(
                result.wavelengths, np.arange(len(labels), dtype=float),
                result.absorbance, "autoqy_epsilon", 0,
            )
            return (_pack(dataset, labels, [paths[0].name], result.concentrations_m,
                          result.path_lengths_cm),
                    f"Restored {len(labels)} processed spectrum/spectra from {paths[0].name}.")
        except (UnicodeDecodeError, ValueError):
            pass
    loaded = []
    for path in paths:
        selected_format = "avantes_abs8" if path.suffix.lower() == ".abs8" else "auto"
        loaded.append((load_spectral_bytes(path.read_bytes(), selected_format), path.name))
    dataset, labels, resampled = _combine_loaded(loaded)
    missing = sum(item.interpolated_values for item, _ in loaded)
    notes = []
    if missing:
        notes.append(f"interpolated {missing} non-finite detector value(s)")
    if resampled:
        notes.append(f"resampled {resampled} spectrum/spectra to the first common grid")
    note = f" ({'; '.join(notes)})" if notes else ""
    return (_pack(dataset, labels, [path.name for path in paths]),
            f"Loaded {len(paths)} file(s), {len(labels)} spectrum/spectra, "
            f"and {len(dataset.wavelengths)} common wavelengths{note}.")


def _choose_files(initial_directory=None):
    initial = _powershell_quote(str(initial_directory or ""))
    script = (
        "Add-Type -AssemblyName System.Windows.Forms; "
        + _foreground_owner_script() +
        "$dialog = New-Object System.Windows.Forms.OpenFileDialog; "
        "$dialog.Multiselect = $true; "
        "$dialog.Filter = 'Spectral files|*.dat;*.txt;*.tsv;*.csv;*.Abs8|All files|*.*'; "
        f"if ('{initial}') {{ $dialog.InitialDirectory = '{initial}' }}; "
        "$result = $dialog.ShowDialog($owner); "
        "if ($result -eq [System.Windows.Forms.DialogResult]::OK) { "
        "$dialog.FileNames | ForEach-Object { Write-Output $_ } }; "
        "$owner.Close(); $owner.Dispose()"
    )
    return _run_powershell_dialog(script)


def _choose_folder(initial_directory=None):
    initial = _powershell_quote(str(initial_directory or ""))
    script = (
        "Add-Type -AssemblyName System.Windows.Forms; "
        + _foreground_owner_script() +
        "$dialog = New-Object System.Windows.Forms.FolderBrowserDialog; "
        f"if ('{initial}') {{ $dialog.SelectedPath = '{initial}' }}; "
        "$result = $dialog.ShowDialog($owner); "
        "if ($result -eq [System.Windows.Forms.DialogResult]::OK) { "
        "Write-Output $dialog.SelectedPath }; $owner.Close(); $owner.Dispose()"
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


def _save_path(folder, filename):
    if not folder:
        raise ValueError("Choose a save folder")
    directory = Path(folder).expanduser()
    if not directory.is_dir():
        raise ValueError(f"Save folder does not exist: {directory}")
    return directory / filename


def _csv_filename(value, default):
    value = str(value or "").strip()
    if not value:
        return default
    candidate = Path(value)
    if candidate.name != value or candidate.parent != Path("."):
        raise ValueError("Enter a file name only, without a folder path")
    if not candidate.suffix:
        return value + ".csv"
    if candidate.suffix.lower() != ".csv":
        raise ValueError("The export file name must use the .csv extension")
    return value


def _export_csv_payload(epsilon_data, processed_data, nonnegative=False):
    if epsilon_data:
        result, labels = _unpack_epsilon(epsilon_data)
        if nonnegative:
            result = EpsilonResult(
                result.wavelengths,
                np.maximum(result.absorbance, 0.0),
                result.concentrations_m,
                result.path_lengths_cm,
                np.maximum(result.individual, 0.0),
                np.maximum(result.mean, 0.0),
                result.standard_deviation,
                result.standard_error,
            )
        return "epsilon-spectra-reactant.csv", export_epsilon_csv(result, labels)
    if processed_data:
        dataset = _unpack(processed_data)
        absorbance = (np.maximum(dataset.absorbance, 0.0)
                      if nonnegative else dataset.absorbance)
        return (
            "processed-absorbance.csv",
            export_smoothed_text(dataset, absorbance, "csv"),
        )
    raise ValueError("Process spectral data before exporting")


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
    if (len(contents) != 1 or
            Path(filenames[0] or "").suffix.lower() not in {".csv", ".tsv"}):
        return None
    payload = base64.b64decode(contents[0].split(",", 1)[1])
    try:
        text = payload.decode("utf-8-sig")
        return load_epsilon_table(text)
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
        raise ValueError("NMR spectra must span the full reactant ε wavelength range")
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
                          "reconstructed_reactant", "reconstructed_pss",
                          "reactant_lower", "reactant_upper", "product",
                          "product_lower", "product_upper", "product_error_lower",
                          "product_error_upper")} | {
        "negative_product_points": result.negative_product_points,
        "negative_bound_points": result.negative_bound_points,
    }


def _unpack_nmr(data):
    return NMRSubtractionResult(
        *(np.asarray(data[field], float) for field in
          ("wavelengths", "normalized_reactant", "normalized_pss",
           "reconstructed_reactant", "reconstructed_pss",
           "reactant_lower", "reactant_upper", "product",
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


def _legend_labels(labels, text):
    """Return plot-only labels edited one-per-line, retaining missing defaults."""
    defaults = list(labels)
    if not text or not str(text).strip():
        return defaults
    edited = str(text).splitlines()
    return [
        edited[index].strip() if index < len(edited) and edited[index].strip() else label
        for index, label in enumerate(defaults)
    ]


def _absorbance_figure(go, dataset, original, processed, labels, method,
                       svd_enabled=None, svd_rank=None, show_legend=True,
                       minimal_colors=False):
    figure = go.Figure()
    colors = _spectrum_colors(len(labels), minimal_colors)
    changed = not np.allclose(original, processed)
    for index, label in enumerate(labels):
        color = colors[index]
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
    return _style(figure, 520, show_legend)


def _epsilon_figure(go, make_subplots, dataset, original, result, labels, method,
                    svd_enabled=None, svd_rank=None, show_legend=True,
                    minimal_colors=False):
    figure = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        row_heights=[0.34, 0.66],
        subplot_titles=("Processed absorbance spectra", "Molar absorptivity"),
    )
    colors = _spectrum_colors(len(labels), minimal_colors)
    for index, label in enumerate(labels):
        color = colors[index]
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
        line={"color": PLOT_PURPLE, "width": 2.4}, name="Mean ε",
    ), row=2, col=1)
    figure.update_yaxes(title_text="Absorbance", row=1, col=1)
    figure.update_yaxes(title_text="ε (M⁻¹ cm⁻¹)", row=2, col=1)
    figure.update_xaxes(title_text="Wavelength (nm)", row=2, col=1)
    figure.update_layout(title={"text": _processing_title(
        method, svd_enabled, svd_rank), "x": 0.02})
    return _style(figure, 690, show_legend)


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
        name="Reactant processed", line={"color": PLOT_BLUE}), row=1, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.normalized_pss, mode="lines",
        name="PSS processed", line={"color": PLOT_ORANGE}), row=1, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.reactant_lower, mode="lines",
        line={"width": 0}, hoverinfo="skip", showlegend=False,
        legendgroup="reactant-error"), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.reactant_upper, mode="lines", fill="tonexty",
        fillcolor="rgba(45,111,142,.20)", line={"width": 0},
        name="Reactant ε SD envelope", legendgroup="reactant-error"), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.product_lower, mode="lines",
        line={"width": 0}, hoverinfo="skip", showlegend=False,
        legendgroup="product-error"), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.product_upper, mode="lines", fill="tonexty",
        fillcolor="rgba(138,90,155,.24)", line={"width": 0},
        name="Product propagation envelope", legendgroup="product-error"), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.reconstructed_reactant, mode="lines",
        name="Reactant ε", line={"color": PLOT_BLUE, "dash": "dot"}), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.reconstructed_pss, mode="lines",
        name="PSS reconstructed ε", line={"color": PLOT_ORANGE, "dash": "dash"}), row=2, col=1)
    figure.add_trace(go.Scatter(
        x=result.wavelengths, y=result.product, mode="lines",
        name="Product ε", line={"color": PLOT_PURPLE, "width": 2.4}), row=2, col=1)
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
    return list(ANALYSIS_TRACE_PALETTE)


def _spectrum_colors(count, minimal=False):
    """Return full or initial/intermediate/final spectrum colors."""
    count = max(0, int(count))
    if not minimal:
        palette = _colors()
        return [palette[index % len(palette)] for index in range(count)]
    if count == 0:
        return []
    if count == 1:
        return [PLOT_BLUE]
    return [PLOT_BLUE] + ["rgba(108,114,128,0.38)"] * (count - 2) + [PLOT_ORANGE]


def _graph_config():
    return {
        "displaylogo": False,
        "modeBarButtonsToRemove": ["toImage"],
    }


def _empty(go, message):
    figure = go.Figure()
    figure.add_annotation(text=message, showarrow=False)
    return _style(figure, 430)


def _style(figure, height, show_legend=True):
    figure.update_layout(
        template="plotly_white", height=height + 113,
        title={"x": 0.02, "y": 0.98, "yanchor": "top", "yref": "container"},
        margin=dict(l=68, r=24, t=185, b=52), hovermode="closest",
        legend={"orientation": "h", "y": 1.06, "yanchor": "bottom",
                "x": 0, "xanchor": "left", "bgcolor": "rgba(255,255,255,.92)"},
        showlegend=show_legend,
    )
    return figure


def run_server(host="127.0.0.1", port=8051, open_browser=True):
    serve_gui(create_app, host, port, open_browser, "AutoQY Spectral Treatment")


def main(argv=None):
    parser = ArgumentParser(prog="autoqy-smoother-gui")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8051, type=int)
    parser.add_argument("--no-open", action="store_true")
    args = parser.parse_args(argv)
    run_server(args.host, args.port, not args.no_open)


if __name__ == "__main__":
    main()
