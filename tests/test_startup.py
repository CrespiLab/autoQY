"""Assert deferred imports in fresh interpreters, independent of test order."""

from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch


class StartupTests(unittest.TestCase):
    def isolated(self, code):
        root = str(Path(__file__).parents[1])
        result = subprocess.run(
            [sys.executable, '-B', '-c',
             f'import sys; sys.path.insert(0, {root!r}); ' + code],
            capture_output=True, text=True, timeout=90,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_package_and_window_launcher_do_not_import_scientific_stack(self):
        self.isolated(
            'import autoqy_core; import autoqy_core.gui_window; '
            'assert not ({"numpy", "scipy", "pandas", "matplotlib"} & sys.modules.keys()); '
            'assert "run_analysis" in dir(autoqy_core)'
        )

    def test_public_exports_keep_their_original_objects(self):
        self.isolated(
            'import autoqy_core; from autoqy_core.config import AnalysisConfig; '
            'assert autoqy_core.AnalysisConfig is AnalysisConfig; '
            'assert autoqy_core.AnalysisConfig is autoqy_core.AnalysisConfig; '
            'assert all(getattr(autoqy_core, name) is not None for name in autoqy_core.__all__)'
        )

    def test_analysis_initial_page_does_not_load_solvers_or_matplotlib(self):
        self.isolated(
            'from autoqy_core.tools.analysis_gui import create_app; '
            'app = create_app(); client = app.server.test_client(); '
            'assert all(client.get(p).status_code == 200 for p in '
            '("/", "/_dash-layout", "/_dash-dependencies")); '
            'assert "autoqy_core.runner" not in sys.modules; '
            'assert "scipy" not in sys.modules; '
            'assert "matplotlib.pyplot" not in sys.modules'
        )

    def test_runner_defers_static_plotting_until_export(self):
        self.isolated(
            'from autoqy_core.runner import run_analysis; '
            'assert "matplotlib.pyplot" not in sys.modules'
        )

    def test_ready_server_has_no_fixed_browser_delay(self):
        from autoqy_core.gui_window import GuiWindowSession
        with patch('autoqy_core.gui_window.Timer') as timer:
            session = GuiWindowSession('http://127.0.0.1:8052')
            session.open_after()
            self.assertEqual(timer.call_args.args[0], 0)
            session.close()


if __name__ == '__main__':
    unittest.main()
