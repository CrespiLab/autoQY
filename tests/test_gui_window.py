import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

from autoqy_core.gui_window import (
    GuiWindowHandle,
    GuiWindowSession,
    _app_window_command,
    _local_gui_url,
    _minimize_console_window,
    _remove_profile_directory,
    open_gui_window,
)


class GuiWindowTests(unittest.TestCase):
    def test_console_minimization_is_a_noop_outside_windows(self):
        with patch("autoqy_core.gui_window.sys.platform", "linux"):
            self.assertFalse(_minimize_console_window())

    def test_windows_console_is_minimized(self):
        kernel32 = Mock()
        kernel32.GetConsoleWindow.return_value = 123
        user32 = Mock()
        with (
            patch("autoqy_core.gui_window.sys.platform", "win32"),
            patch("autoqy_core.gui_window.ctypes.WinDLL",
                  side_effect=[kernel32, user32]),
        ):
            self.assertTrue(_minimize_console_window())
        user32.ShowWindow.assert_called_once_with(123, 6)

    def test_profile_cleanup_accepts_an_already_removed_directory(self):
        with TemporaryDirectory() as temporary:
            missing = Path(temporary) / "missing"
            self.assertTrue(_remove_profile_directory(missing, attempts=1, delay=0))

    def test_local_url_uses_a_reachable_loopback_address(self):
        self.assertEqual(_local_gui_url("0.0.0.0", 8050), "http://127.0.0.1:8050")
        self.assertEqual(_local_gui_url("::1", 8051), "http://[::1]:8051")

    def test_app_command_uses_a_dedicated_window_and_profile(self):
        command = _app_window_command(
            Path("C:/Program Files/Microsoft/Edge/Application/msedge.exe"),
            "http://127.0.0.1:8051", Path("C:/Temp/autoqy-profile"),
        )
        self.assertIn("--app=http://127.0.0.1:8051", command)
        self.assertIn("--new-window", command)
        self.assertIn("--user-data-dir=C:\\Temp\\autoqy-profile", command)
        self.assertIn("--disable-background-mode", command)

    @patch("autoqy_core.gui_window.webbrowser.open_new")
    @patch("autoqy_core.gui_window._find_app_browser", return_value=None)
    def test_default_browser_is_the_fallback(self, _, open_new):
        handle = open_gui_window("http://127.0.0.1:8050", "AutoQY Power")
        open_new.assert_called_once_with("http://127.0.0.1:8050")
        self.assertEqual(handle.mode, "browser")

    def test_app_launch_tracks_and_cleans_its_process_tree(self):
        with TemporaryDirectory() as temporary:
            profile = Path(temporary) / "profile"
            profile.mkdir()
            process = Mock()
            process._handle = 123
            process.wait.return_value = 0
            job = Mock()
            with (
                patch("autoqy_core.gui_window._find_app_browser",
                      return_value=Path("C:/edge.exe")),
                patch("autoqy_core.gui_window.mkdtemp", return_value=str(profile)),
                patch("autoqy_core.gui_window.subprocess.Popen", return_value=process) as popen,
                patch("autoqy_core.gui_window._WindowsJob", return_value=job),
            ):
                handle = open_gui_window(
                    "http://127.0.0.1:8052", "AutoQY Analysis"
                )
            self.assertEqual(handle.mode, "app")
            job.assign.assert_called_once_with(process)
            command = popen.call_args.args[0]
            self.assertIn("--app=http://127.0.0.1:8052", command)
            self.assertIn(f"--user-data-dir={profile}", command)
            handle.close()
            job.close.assert_called_once_with()
            process.wait.assert_called_once_with(timeout=4)
            self.assertFalse(profile.exists())

    def test_app_window_without_a_job_still_terminates_cleanly(self):
        process = Mock()
        process.poll.return_value = None
        process.wait.return_value = 0
        handle = GuiWindowHandle("app", process=process)
        handle.close()
        process.terminate.assert_called_once_with()
        process.wait.assert_called_once_with(timeout=4)

    def test_session_closes_the_opened_window(self):
        handle = Mock()
        opener = Mock(return_value=handle)
        session = GuiWindowSession(
            "http://127.0.0.1:8051", "AutoQY Spectral Treatment", opener
        )
        session._open()
        opener.assert_called_once_with(
            "http://127.0.0.1:8051", "AutoQY Spectral Treatment"
        )
        session.close()
        handle.close.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
