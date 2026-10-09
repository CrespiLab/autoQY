"""Atomic, disk-backed recovery for the spectral treatment session."""

import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from threading import Lock


def recovery_path():
    root = Path(os.environ.get("LOCALAPPDATA") or Path.home() / ".local" / "share")
    return root / "AutoQY" / "spectral-treatment" / "autosave.json"


class SpectralRecovery:
    def __init__(self, path=None):
        self.path = Path(path) if path is not None else recovery_path()
        self.lock = Lock()

    def read(self):
        with self.lock:
            if not self.path.exists():
                return None
            snapshot = json.loads(self.path.read_text(encoding="utf-8"))
            if not isinstance(snapshot, dict) or snapshot.get("version") != 1:
                raise ValueError("Unsupported spectral recovery version")
            if not isinstance(snapshot.get("settings", {}), dict):
                raise ValueError("Invalid spectral recovery settings")
            for key in ("dataset", "nmr"):
                if snapshot.get(key) is not None and not isinstance(snapshot[key], dict):
                    raise ValueError(f"Invalid spectral recovery {key}")
            return snapshot

    def write(self, snapshot):
        with self.lock:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            temporary = None
            try:
                with NamedTemporaryFile(mode="w", encoding="utf-8", delete=False,
                                        dir=self.path.parent, suffix=".tmp") as stream:
                    temporary = Path(stream.name)
                    json.dump(snapshot, stream, separators=(",", ":"), allow_nan=False)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary, self.path)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
