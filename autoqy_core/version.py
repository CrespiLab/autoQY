"""Project version shared by the browser GUIs."""

from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version as installed_version
from pathlib import Path
import tomllib


@lru_cache(maxsize=1)
def get_project_version() -> str:
    """Read the build version from pyproject.toml, with a wheel-install fallback."""
    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    if pyproject.is_file():
        with pyproject.open("rb") as source:
            value = tomllib.load(source)["project"]["version"]
        if not isinstance(value, str) or not value.strip():
            raise ValueError("pyproject.toml contains no valid project.version")
        return value.strip()

    try:
        return installed_version("autoqy-core")
    except PackageNotFoundError as error:
        raise RuntimeError("AutoQY version metadata is unavailable") from error
