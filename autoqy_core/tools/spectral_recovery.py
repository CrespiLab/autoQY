"""Atomic, disk-backed recovery for the spectral treatment session."""

import json
import base64
from io import BytesIO
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from threading import Lock
from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED, BadZipFile

import numpy as np


SESSION_MAGIC = b"AUTOQY-SESSION\x00\x01\n"


def validate_snapshot(snapshot):
    if not isinstance(snapshot, dict) or snapshot.get("version") != 1:
        raise ValueError("Unsupported spectral session version")
    if not isinstance(snapshot.get("settings", {}), dict):
        raise ValueError("Invalid spectral session settings")
    for key in ("dataset", "nmr", "isomer_request"):
        if snapshot.get(key) is not None and not isinstance(snapshot[key], dict):
            raise ValueError(f"Invalid spectral session {key}")
    return snapshot


def encode_session(snapshot):
    """Versioned archive of binary arrays and small JSON metadata; never pickle."""
    validate_snapshot(snapshot)
    output = BytesIO()
    output.write(SESSION_MAGIC)
    with ZipFile(output, "w", compression=ZIP_DEFLATED) as archive:
        counter = 0
        def encode(value):
            nonlocal counter
            if isinstance(value, dict):
                result = {}
                for key, item in value.items():
                    if key == "absorbance_data" and value.get("absorbance_encoding") == "zlib-base64-float64":
                        name = f"arrays/{counter}.bin"
                        counter += 1
                        archive.writestr(name, base64.b64decode(item, validate=True), compress_type=ZIP_STORED)
                        result[key] = {"__autoqy_bytes__": name}
                    else:
                        result[key] = encode(item)
                return result
            if isinstance(value, list) and value:
                try:
                    array = np.asarray(value)
                except ValueError:
                    array = None
                if array is not None and array.dtype.kind in "biuf" and np.isfinite(array).all():
                    name = f"arrays/{counter}.npy"
                    counter += 1
                    stream = BytesIO()
                    np.save(stream, array, allow_pickle=False)
                    archive.writestr(name, stream.getvalue())
                    return {"__autoqy_array__": name}
                return [encode(item) for item in value]
            return value
        archive.writestr("metadata.json", json.dumps(encode(snapshot), separators=(",", ":"), allow_nan=False))
    return output.getvalue()


def decode_session(payload):
    if not payload.startswith(SESSION_MAGIC):
        raise ValueError("Not an AutoQY binary session (.aqys)")
    try:
        with ZipFile(BytesIO(payload[len(SESSION_MAGIC):])) as archive:
            def decode(value):
                if isinstance(value, dict):
                    for tag in ("__autoqy_array__", "__autoqy_bytes__"):
                        if set(value) == {tag}:
                            name = value[tag]
                            if not isinstance(name, str) or not name.startswith("arrays/") or ".." in name:
                                raise ValueError("Invalid session array reference")
                            data = archive.read(name)
                            if tag == "__autoqy_bytes__":
                                return base64.b64encode(data).decode("ascii")
                            array = np.load(BytesIO(data), allow_pickle=False)
                            if array.dtype.kind not in "biuf" or not np.isfinite(array).all():
                                raise ValueError("Invalid session array")
                            return array.tolist()
                    return {key: decode(item) for key, item in value.items()}
                if isinstance(value, list):
                    return [decode(item) for item in value]
                return value
            return validate_snapshot(decode(json.loads(archive.read("metadata.json"))))
    except (BadZipFile, KeyError, EOFError) as error:
        raise ValueError("Corrupt or incomplete AutoQY session") from error


def recovery_path():
    root = Path(os.environ.get("LOCALAPPDATA") or Path.home() / ".local" / "share")
    return root / "AutoQY" / "spectral-treatment" / "autosave.aqys"


class SpectralRecovery:
    def __init__(self, path=None):
        self.path = Path(path) if path is not None else recovery_path()
        self.lock = Lock()

    def read(self):
        with self.lock:
            source = self.path
            if not source.exists() and source.name == "autosave.aqys":
                source = source.with_name("autosave.json")
            if not source.exists():
                return None
            payload = source.read_bytes()
            return (decode_session(payload) if payload.startswith(SESSION_MAGIC)
                    else validate_snapshot(json.loads(payload.decode("utf-8"))))

    def write(self, snapshot):
        with self.lock:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            temporary = None
            try:
                with NamedTemporaryFile(mode="wb", delete=False,
                                        dir=self.path.parent, suffix=".tmp") as stream:
                    temporary = Path(stream.name)
                    stream.write(encode_session(snapshot))
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary, self.path)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
