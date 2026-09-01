"""Launch local AutoQY GUIs in dedicated browser application windows."""

from dataclasses import dataclass
import ctypes
from ctypes import wintypes
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import mkdtemp
from threading import Lock, Thread, Timer
import time
import webbrowser


_WINDOWS_BROWSER_LOCATIONS = (
    ("PROGRAMFILES(X86)", "Microsoft/Edge/Application/msedge.exe"),
    ("PROGRAMFILES", "Microsoft/Edge/Application/msedge.exe"),
    ("LOCALAPPDATA", "Microsoft/Edge/Application/msedge.exe"),
    ("PROGRAMFILES", "Google/Chrome/Application/chrome.exe"),
    ("PROGRAMFILES(X86)", "Google/Chrome/Application/chrome.exe"),
    ("LOCALAPPDATA", "Google/Chrome/Application/chrome.exe"),
)


def _remove_profile_directory(directory, attempts=20, delay=0.2):
    """Remove a browser profile after Windows releases its short-lived locks."""
    if directory is None:
        return True
    directory = Path(directory)
    for attempt in range(attempts):
        try:
            shutil.rmtree(directory)
            return True
        except FileNotFoundError:
            return True
        except OSError:
            if attempt + 1 < attempts:
                time.sleep(delay)
    return not directory.exists()


def _find_app_browser():
    """Return Edge or Chrome when a Windows app-mode browser is available."""
    if sys.platform != "win32":
        return None
    for executable in ("msedge.exe", "chrome.exe"):
        discovered = shutil.which(executable)
        if discovered:
            return Path(discovered)
    for environment_name, relative_path in _WINDOWS_BROWSER_LOCATIONS:
        root = os.environ.get(environment_name)
        if root:
            candidate = Path(root) / relative_path
            if candidate.is_file():
                return candidate
    return None


def _app_window_command(browser, url, profile_directory):
    return [
        str(browser),
        f"--app={url}",
        "--new-window",
        f"--user-data-dir={profile_directory}",
        "--no-first-run",
        "--no-default-browser-check",
        "--disable-background-mode",
        "--disable-extensions",
        "--disable-sync",
        "--disable-features=msEdgeFirstRunExperience",
    ]


def _local_gui_url(host, port):
    browser_host = "127.0.0.1" if host in {"0.0.0.0", "::"} else host
    if ":" in browser_host and not browser_host.startswith("["):
        browser_host = f"[{browser_host}]"
    return f"http://{browser_host}:{int(port)}"


class _WindowsJob:
    """Keep the isolated browser process tree attached to the GUI terminal."""

    _KILL_ON_JOB_CLOSE = 0x00002000
    _EXTENDED_LIMIT_INFORMATION = 9

    class _BasicLimitInformation(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class _IoCounters(ctypes.Structure):
        _fields_ = [
            ("ReadOperationCount", ctypes.c_ulonglong),
            ("WriteOperationCount", ctypes.c_ulonglong),
            ("OtherOperationCount", ctypes.c_ulonglong),
            ("ReadTransferCount", ctypes.c_ulonglong),
            ("WriteTransferCount", ctypes.c_ulonglong),
            ("OtherTransferCount", ctypes.c_ulonglong),
        ]

    class _ExtendedLimitInformation(ctypes.Structure):
        pass

    _ExtendedLimitInformation._fields_ = [
        ("BasicLimitInformation", _BasicLimitInformation),
        ("IoInfo", _IoCounters),
        ("ProcessMemoryLimit", ctypes.c_size_t),
        ("JobMemoryLimit", ctypes.c_size_t),
        ("PeakProcessMemoryUsed", ctypes.c_size_t),
        ("PeakJobMemoryUsed", ctypes.c_size_t),
    ]

    def __init__(self):
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.CreateJobObjectW.restype = wintypes.HANDLE
        kernel32.SetInformationJobObject.argtypes = (
            wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD
        )
        kernel32.AssignProcessToJobObject.argtypes = (wintypes.HANDLE, wintypes.HANDLE)
        kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
        handle = kernel32.CreateJobObjectW(None, None)
        if not handle:
            raise ctypes.WinError(ctypes.get_last_error())
        information = self._ExtendedLimitInformation()
        information.BasicLimitInformation.LimitFlags = self._KILL_ON_JOB_CLOSE
        configured = kernel32.SetInformationJobObject(
            handle, self._EXTENDED_LIMIT_INFORMATION,
            ctypes.byref(information), ctypes.sizeof(information),
        )
        if not configured:
            error = ctypes.get_last_error()
            kernel32.CloseHandle(handle)
            raise ctypes.WinError(error)
        self._kernel32 = kernel32
        self._handle = handle

    def assign(self, process):
        assigned = self._kernel32.AssignProcessToJobObject(
            self._handle, wintypes.HANDLE(process._handle)
        )
        if not assigned:
            raise ctypes.WinError(ctypes.get_last_error())

    def close(self):
        if self._handle:
            self._kernel32.CloseHandle(self._handle)
            self._handle = None


@dataclass
class GuiWindowHandle:
    """Resources associated with one opened GUI window."""

    mode: str
    process: object = None
    job: object = None
    profile_directory: Path | None = None

    def close(self):
        if self.job is not None:
            self.job.close()
            self.job = None
        elif self.process is not None and self.process.poll() is None:
            self.process.terminate()
        if self.process is not None:
            try:
                self.process.wait(timeout=4)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=2)
        if self.profile_directory is not None:
            _remove_profile_directory(self.profile_directory)
            self.profile_directory = None


def open_gui_window(url, window_name="AutoQY"):
    """Open a dedicated app window, falling back to the default browser."""
    browser = _find_app_browser()
    if browser is None:
        webbrowser.open_new(url)
        return GuiWindowHandle("browser")

    safe_name = "".join(character if character.isalnum() else "-"
                        for character in window_name).strip("-").lower() or "gui"
    profile_directory = Path(mkdtemp(prefix=f"autoqy-{safe_name}-"))
    flags = (getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) |
             getattr(subprocess, "CREATE_NO_WINDOW", 0))
    try:
        process = subprocess.Popen(
            _app_window_command(browser, url, profile_directory),
            creationflags=flags,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except OSError:
        _remove_profile_directory(profile_directory)
        webbrowser.open_new(url)
        return GuiWindowHandle("browser")

    job = None
    if sys.platform == "win32":
        try:
            job = _WindowsJob()
            job.assign(process)
        except OSError:
            if job is not None:
                job.close()
            job = None
    return GuiWindowHandle("app", process, job, profile_directory)


class GuiWindowSession:
    """Schedule and clean up one local GUI window around a blocking server."""

    def __init__(self, url, window_name="AutoQY", opener=None):
        self.url = url
        self.window_name = window_name
        self._opener = opener or open_gui_window
        self._lock = Lock()
        self._timer = None
        self._handle = None
        self._closed = False

    def open_after(self, delay=1.0):
        with self._lock:
            if self._closed:
                return
            self._timer = Timer(delay, self._open)
            self._timer.daemon = True
            self._timer.start()

    def _open(self):
        handle = self._opener(self.url, self.window_name)
        with self._lock:
            if self._closed:
                handle.close()
            else:
                self._handle = handle

    def close(self):
        with self._lock:
            self._closed = True
            timer, self._timer = self._timer, None
            handle, self._handle = self._handle, None
        if timer is not None:
            timer.cancel()
        if handle is not None:
            handle.close()


def serve_gui(create_app, host, port, open_window, window_name):
    """Serve a Dash GUI and keep its dedicated window lifecycle synchronized."""
    from werkzeug.serving import make_server

    app = create_app()
    server = make_server(host, port, app.server, threaded=True)
    window_state, window_lock = app.server.config["AUTOQY_WINDOW_STATE"]

    def close_after_window():
        while True:
            time.sleep(0.5)
            with window_lock:
                requested = window_state["close_requested_at"]
                last_heartbeat = window_state["last_heartbeat_at"]
            now = time.monotonic()
            close_requested = requested is not None and now - requested >= 3
            heartbeat_expired = last_heartbeat is not None and now - last_heartbeat >= 5
            if close_requested or heartbeat_expired:
                server.shutdown()
                return

    Thread(target=close_after_window, daemon=True).start()
    window = GuiWindowSession(_local_gui_url(host, port), window_name)
    if open_window:
        window.open_after()
    try:
        server.serve_forever()
    finally:
        window.close()
