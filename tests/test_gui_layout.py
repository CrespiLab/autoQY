import tomllib
import unittest
from pathlib import Path

import numpy as np

from autoqy_core.version import get_project_version

try:
    from dash import dcc, html

    from autoqy_core.power_web import create_app as create_power_app
    from autoqy_core.tools.analysis_gui import (
        _pss_card,
        _render_nipe_window_analysis,
        _nipe_headline_warning,
        _spectra_led_figure,
        create_app as create_analysis_app,
    )
    from autoqy_core.tools.smoother_gui import create_app as create_spectral_app
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    dcc = html = go = make_subplots = None


def _components(root):
    if isinstance(root, (list, tuple)):
        for child in root:
            yield from _components(child)
        return
    yield root
    children = getattr(root, "children", None)
    if children is not None:
        yield from _components(children)


def _by_id(root, component_id):
    return next(
        component for component in _components(root)
        if getattr(component, "id", None) == component_id
    )


class ProjectVersionTests(unittest.TestCase):
    def test_version_is_read_from_pyproject(self):
        pyproject = Path(__file__).parents[1] / "pyproject.toml"
        with pyproject.open("rb") as source:
            expected = tomllib.load(source)["project"]["version"]
        self.assertEqual(get_project_version(), expected)


@unittest.skipUnless(html is not None, "Dash GUI dependencies are not installed")
class GuiLayoutTests(unittest.TestCase):
    def test_every_gui_displays_the_project_version(self):
        expected = f"Version {get_project_version()}"
        for app in (create_power_app(), create_spectral_app(), create_analysis_app()):
            visible_text = " ".join(
                component.children
                for component in _components(app.layout)
                if isinstance(getattr(component, "children", None), str)
            )
            self.assertIn(expected, visible_text)
            self.assertNotIn("AUTOQY CORE", visible_text.upper())

    def test_every_gui_closes_its_page_when_the_terminal_stops(self):
        for app in (create_power_app(), create_spectral_app(), create_analysis_app()):
            self.assertIn("window.setInterval(heartbeat, 1000)", app.index_string)
            self.assertIn("window.close()", app.index_string)
            self.assertIn("AutoQY GUI stopped", app.index_string)
            self.assertIn("restart the AutoQY GUI to continue", app.index_string)

    def test_power_gui_stops_its_server_when_the_browser_closes(self):
        app = create_power_app()
        routes = {rule.rule for rule in app.server.url_map.iter_rules()}
        self.assertIn("/_autoqy_power_heartbeat", routes)
        self.assertIn("/_autoqy_power_window_closed", routes)
        self.assertIn("AUTOQY_WINDOW_STATE", app.server.config)
        state, _ = app.server.config["AUTOQY_WINDOW_STATE"]
        self.assertIn("last_heartbeat_at", state)

    def test_spectral_loader_uses_the_analysis_spinner(self):
        app = create_spectral_app()
        loading = next(
            component for component in _components(app.layout)
            if isinstance(component, dcc.Loading)
            and any(getattr(child, "id", None) == "load-message"
                    for child in _components(component.children))
        )
        self.assertEqual(loading.type, "circle")

    def test_spectral_default_filename_copy_is_short(self):
        app = create_spectral_app()
        self.assertEqual(_by_id(app.layout, "save-filename").placeholder, "Default name")

    def test_spectral_export_buttons_precede_the_plot_options(self):
        app = create_spectral_app()
        toolbar = next(
            component for component in _components(app.layout)
            if getattr(component, "className", None) == "plot-toolbar"
            and any(getattr(child, "id", None) == "save-epsilon-png"
                    for child in _components(component))
        )
        first_ids = [
            getattr(component, "id", None)
            for component in _components(toolbar.children[0])
        ]
        self.assertIn("save-epsilon-png", first_ids)
        self.assertIn("save-epsilon-svg", first_ids)

    def test_analysis_has_no_validate_button_and_explains_json_saving(self):
        app = create_analysis_app()
        component_ids = [getattr(component, "id", None) for component in _components(app.layout)]
        self.assertNotIn("validate-analysis", component_ids)
        visible_text = " ".join(
            component.children
            for component in _components(app.layout)
            if isinstance(getattr(component, "children", None), str)
        )
        self.assertIn("does not save analysis.json", visible_text)
        self.assertNotIn("validate-analysis", " ".join(app.callback_map))

    def test_analysis_preprocessing_preview_is_available_before_the_fit(self):
        app = create_analysis_app()
        tabs = _by_id(app.layout, "analysis-plot-tabs")
        labels = [tab.label for tab in tabs.children]
        self.assertIn("Preprocessing", labels)
        self.assertIn("Endpoint reconstruction", labels)
        self.assertNotIn("References and reconstruction", labels)
        preview = _by_id(app.layout, "spectra-figure")
        self.assertEqual(preview.style["height"], "900px")
        self.assertTrue(any(
            key.startswith("spectra-figure.figure@")
            for key in app.callback_map
        ))

    def test_analysis_offers_chemical_actinometer_photon_flux(self):
        app = create_analysis_app()
        actinometer = _by_id(
            app.layout, {"type": "analysis-field", "name": "chemical_actinometer"}
        )
        photon_flux = _by_id(
            app.layout, {"type": "analysis-field", "name": "photon_flux_mol_s"}
        )
        self.assertEqual(actinometer.options[0]["value"], "on")
        self.assertIn("chemical actinometer", actinometer.options[0]["label"])
        self.assertTrue(photon_flux.disabled)

    def test_analysis_offers_specord_and_cary_binary_spectra(self):
        app = create_analysis_app()
        formats = _by_id(
            app.layout,
            {"type": "analysis-field", "name": "format_measurement_spectra"},
        )
        values = {option["value"] for option in formats.options}
        self.assertIn("specord", values)
        self.assertIn("agilent_cary", values)

        spectral_app = create_spectral_app()
        visible_text = " ".join(
            component.children
            for component in _components(spectral_app.layout)
            if isinstance(getattr(component, "children", None), str)
        )
        self.assertIn("SPECORD", visible_text)
        self.assertIn("Cary .DSW/.BSW", visible_text)

    def test_analysis_offers_nipe_as_a_fit_method(self):
        app = create_analysis_app()
        methods = _by_id(
            app.layout, {"type": "analysis-field", "name": "fit_method"}
        )
        self.assertIn("nipe", {option["value"] for option in methods.options})

    def test_analysis_has_collapsible_nipe_window_report(self):
        app = create_analysis_app()
        _by_id(app.layout, "nipe-headline-warning")
        panel = _by_id(app.layout, "nipe-window-panel")
        self.assertIsInstance(panel, html.Details)
        self.assertFalse(panel.open)
        _by_id(panel, "nipe-window-analysis")

    def test_nipe_window_report_shows_recommended_yields_and_windows(self):
        summary = {
            "ab_model_assessment": {"status": "stop"},
            "nipe": {"pre_plateau_window_analysis": {
                "plateau_detected": True,
                "plateau_time_s": 240.0,
                "analysis_end_time_s": 210.0,
                "window_point_count": 4,
                "window_duration_s": 90.0,
                "window_count": 1,
                "extrapolated_zero_exposure_yield_percent": {
                    "R_to_P": 14.870101, "P_to_R": 12.38181,
                },
                "extrapolated_standard_error_percent": {
                    "R_to_P": 0.173003, "P_to_R": 0.255767,
                },
                "full_trace_change_percent": {
                    "R_to_P": -16.6, "P_to_R": -28.6,
                },
                "windows": [{
                    "start_s": 0.0, "end_s": 90.0, "midpoint_s": 45.0,
                    "R_to_P_percent": 14.8, "P_to_R_percent": 11.1,
                    "R_to_P_standard_error_percent": 0.2,
                    "P_to_R_standard_error_percent": 0.3,
                    "jacobian_condition": 7.0,
                }],
            }},
        }
        rendered = _render_nipe_window_analysis(html, summary, "F2", "F1")
        visible_text = " ".join(
            component.children for component in _components(rendered)
            if isinstance(getattr(component, "children", None), str)
        )
        self.assertIn("14.87 ± 0.17%", visible_text)
        self.assertIn("12.4 ± 0.3%", visible_text)
        self.assertIn("Sustained flattening was detected at 240 s", visible_text)
        self.assertIn("0–90", visible_text)

    def test_nipe_full_trace_result_warns_user_to_open_window_analysis(self):
        warning = _nipe_headline_warning(html, {
            "fit_method": "nipe",
            "ab_model_assessment": {"status": "stop"},
        })
        visible_text = " ".join(
            child.children if hasattr(child, "children") else str(child)
            for child in warning.children
        )
        self.assertIn("not corrected for degradation", visible_text)
        self.assertIn("Do not report them", visible_text)
        self.assertIn("NIPE pre-plateau window analysis", visible_text)
        self.assertIn("status-stop", warning.className)
        self.assertEqual(_nipe_headline_warning(html, {"fit_method": "emission"}), "")

    def test_pss_distribution_is_visible_beside_quantum_yields(self):
        app = create_analysis_app()
        result_strip = next(
            component for component in _components(app.layout)
            if getattr(component, "className", None) == "result-strip analysis-result-strip"
        )
        result_ids = [getattr(component, "id", None) for component in result_strip.children]
        self.assertEqual(
            result_ids,
            ["result-rp", "result-pr", "result-pss", "result-fit"],
        )
        card = _pss_card(
            html,
            {"extrapolated_pss_percent": {"reactant": 23.3459, "product": 76.6541}},
            "trans", "cis",
        )
        visible_text = " ".join(
            component.children
            for component in _components(card)
            if isinstance(getattr(component, "children", None), str)
        )
        self.assertIn("trans 23.3%", visible_text)
        self.assertIn("cis 76.7%", visible_text)

    def test_analysis_preprocessing_separates_led_from_spectral_decay(self):
        wavelengths = np.array([400.0, 450.0, 500.0])
        reference_wavelengths = np.array([350.0, 400.0, 450.0, 500.0, 550.0])
        absorbance = np.array([
            [1.0, 0.8, 0.6],
            [0.9, 0.7, 0.5],
            [0.8, 0.6, 0.4],
        ])
        figure = _spectra_led_figure(
            go, make_subplots, wavelengths, absorbance,
            (reference_wavelengths, np.array([10_000.0, 100.0, 80.0, 60.0, 9_000.0])),
            (reference_wavelengths, np.array([8_000.0, 20.0, 50.0, 90.0, 7_000.0])),
            wavelengths, np.array([0.1, 1.0, 0.1]), (400.0, 500.0), 450.0,
        )
        traces = {trace.name: trace for trace in figure.data}
        self.assertEqual(traces["Reactant ε"].xaxis, "x")
        self.assertEqual(traces["Product ε"].xaxis, "x")
        self.assertEqual(traces["Processed LED (normalized)"].yaxis, "y2")
        self.assertEqual(traces["Initial spectrum"].xaxis, "x2")
        self.assertEqual(traces["Final spectrum"].xaxis, "x2")
        self.assertFalse(figure.layout.yaxis2.showgrid)
        self.assertEqual(tuple(figure.layout.yaxis.range)[0], 0)
        self.assertEqual(tuple(figure.layout.yaxis2.range)[0], 0)
        self.assertLess(tuple(figure.layout.yaxis.range)[1], 1_000)
        self.assertTrue(np.all(np.asarray(traces["Reactant ε"].x) >= 400.0))
        self.assertTrue(np.all(np.asarray(traces["Reactant ε"].x) <= 500.0))

    def test_nested_panels_have_independent_open_and_closed_symbols(self):
        css = (Path(__file__).parents[1] / "autoqy_core" / "assets" / "power_web.css").read_text(
            encoding="utf-8"
        )
        self.assertIn('.tool-details > summary::after', css)
        self.assertIn('.nested-tool > summary::after', css)
        self.assertIn('.nested-tool[open] > summary::after', css)
        self.assertNotIn('.tool-details summary::after', css)

    def test_analysis_results_and_method_names_wrap_instead_of_clipping(self):
        css = (Path(__file__).parents[1] / "autoqy_core" / "assets" / "analysis_gui.css").read_text(
            encoding="utf-8"
        )
        result_rule = css.split(
            ".analysis-result-strip .result-card strong {", 1
        )[1].split("}", 1)[0]
        species_rule = css.split(
            ".analysis-result-strip .pss-species span {", 1
        )[1].split("}", 1)[0]
        comparison_rule = css.split(
            ".comparison-table th, .comparison-table td {", 1
        )[1].split("}", 1)[0]
        for rule in (result_rule, species_rule, comparison_rule):
            self.assertIn("white-space: normal", rule)
        self.assertIn("overflow-wrap: anywhere", result_rule)
        self.assertNotIn("text-overflow: ellipsis", result_rule)


if __name__ == "__main__":
    unittest.main()
