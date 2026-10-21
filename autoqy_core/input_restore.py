"""Recover embedded analysis inputs as ordinary CSV tables."""

from io import StringIO
from itertools import count
from pathlib import Path

import pandas as pd

from .config import ANALYSIS_INPUTS, input_format
from .io import load_spectra


def restore_embedded_inputs(config, directory):
    """Write embedded inputs only; retain all table columns and never overwrite files.

    Conversion finishes before any output is written. Linked files and unrelated
    analysis settings are not needed. Returned paths are keyed by input role.
    """
    tables = {}
    for name in ANALYSIS_INPUTS:
        if not isinstance(config.values.get("inputs", {}).get(name), dict):
            continue
        try:
            path = config.input_path(name)  # Also verifies the embedded checksum.
            spec = input_format(config, name)
            spec = {"type": spec} if isinstance(spec, str) else dict(spec)
            kind = spec["type"]
            if kind in {"specord", "agilent_cary"} and name != "timestamps":
                wavelengths, spectra = load_spectra(path, spec)
                frame = pd.DataFrame(spectra, columns=[
                    f"Spectrum_{i + 1}" for i in range(spectra.shape[1])
                ])
                frame.insert(0, "Wavelength_nm", wavelengths)
            elif kind in {"generic_delimited", "spectragryph_tsv", "ahk_csv", "simple_csv"}:
                options = spec if kind == "generic_delimited" else {}
                separator = ("\t" if kind == "spectragryph_tsv"
                             else options.get("delimiter", ","))
                if separator == "tab":
                    separator = "\t"
                frame = pd.read_csv(
                    path, sep=separator,
                    header=0 if options.get("header", True) else None,
                    skiprows=options.get("skip_rows", 0),
                    decimal=options.get("decimal", "."),
                    float_precision="round_trip", keep_default_na=False,
                )
                if not options.get("header", True):
                    frame.columns = [f"Column_{i + 1}" for i in range(len(frame.columns))]
            else:
                raise ValueError(f"Unsupported input format: {kind}")
            # Keep every original column, including ε errors, replicate counts,
            # asymmetric ranges, and AHK events, rather than just the fit inputs.
            stream = StringIO()
            frame.to_csv(stream, index=False, lineterminator="\n")
            tables[name] = stream.getvalue()
        except Exception as error:
            filename = config.values["inputs"][name].get("filename", name)
            raise ValueError(f"Cannot restore {name} ({filename}): {error}") from error
    if not tables:
        raise ValueError("This JSON has no input data stored inside it. "
                         "Load a JSON saved with Include all input data in JSON enabled.")

    directory = Path(directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, text in tables.items():
        for number in count(1):
            suffix = "" if number == 1 else f"_{number}"
            path = directory / f"{name}{suffix}.csv"
            try:
                stream = path.open("x", encoding="utf-8", newline="")
            except FileExistsError:
                continue
            with stream:
                stream.write(text)
            paths[name] = path
            break
    return paths
