import tomllib
import unittest
from pathlib import Path

from autoqy_core.version import get_project_version

try:
    from dash import dcc, html

    from autoqy_core.power_web import create_app as create_power_app
    from autoqy_core.tools.analysis_gui import create_app as create_analysis_app
    from autoqy_core.tools.smoother_gui import create_app as create_spectral_app
except ImportError:
    dcc = html = None


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

    def test_nested_panels_have_independent_open_and_closed_symbols(self):
        css = (Path(__file__).parents[1] / "autoqy_core" / "assets" / "power_web.css").read_text(
            encoding="utf-8"
        )
        self.assertIn('.tool-details > summary::after', css)
        self.assertIn('.nested-tool > summary::after', css)
        self.assertIn('.nested-tool[open] > summary::after', css)
        self.assertNotIn('.tool-details summary::after', css)


if __name__ == "__main__":
    unittest.main()
