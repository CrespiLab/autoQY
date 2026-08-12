"""Load AutoQY input files without GUI dependencies."""

from pathlib import Path

import numpy as np
import pandas as pd


AVANTES_HEADER_BYTES = 328


def load_avantes_abs8(path):
    """Load wavelength and absorbance arrays from an AvaSoft 8 Abs8 file."""
    return load_avantes_abs8_bytes(Path(path).read_bytes())


def load_avantes_abs8_bytes(data):
    """Decode a single-channel AVS84 absorbance record from bytes."""
    data = bytes(data)
    if len(data) < AVANTES_HEADER_BYTES or data[:5] != b"AVS84":
        raise ValueError("Not an Avantes AvaSoft 8 AVS84 file")
    if data[11] != 2:
        raise ValueError(
            f"Expected Avantes absorbance mode (2), found measurement mode {data[11]}"
        )
    pixels = int.from_bytes(data[91:93], "little") + 1
    required = AVANTES_HEADER_BYTES + 4 * pixels * np.dtype("<f4").itemsize
    if pixels < 2 or len(data) < required:
        raise ValueError("Avantes Abs8 file is truncated or has an invalid pixel count")

    arrays = np.frombuffer(
        data, dtype="<f4", count=4 * pixels, offset=AVANTES_HEADER_BYTES
    ).reshape(4, pixels).astype(float)
    wavelengths, sample, dark, reference = arrays
    if not np.isfinite(arrays).all():
        raise ValueError("Avantes Abs8 file contains non-finite instrument data")
    if np.any(np.diff(wavelengths) <= 0):
        raise ValueError("Avantes Abs8 wavelengths must increase")

    numerator = sample - dark
    denominator = reference - dark
    valid = (numerator > 0) & (denominator > 0)
    absorbance = np.full(pixels, np.nan)
    absorbance[valid] = -np.log10(numerator[valid] / denominator[valid])
    if np.count_nonzero(valid) < 2:
        raise ValueError("Avantes Abs8 file has fewer than two valid absorbance points")
    return wavelengths, absorbance


def load_spectra(path, format_spec="spectragryph_tsv"):
    spec = _format_spec(format_spec)
    if spec["type"] == "spectragryph_tsv":
        data = pd.read_csv(path, sep="\t").drop(
            columns="Wavenumbers [1/cm]", errors="ignore"
        )
        wavelength_column, value_columns = data.columns[0], list(data.columns[1:])
    elif spec["type"] == "generic_delimited":
        data = _read_delimited(path, spec)
        wavelength_column = _column(data, spec.get("wavelength_column", 0))
        ignored = {_column(data, value) for value in spec.get("ignored_columns", [])}
        requested = spec.get("value_columns", "remaining")
        if requested == "remaining":
            value_columns = [column for column in data.columns
                             if column != wavelength_column and column not in ignored]
        else:
            value_columns = [_column(data, value) for value in requested]
    else:
        raise ValueError(f"Unsupported spectral format: {spec['type']}")

    if not value_columns:
        raise ValueError(f"No spectral value columns found in {path}")
    wavelengths = pd.to_numeric(data[wavelength_column], errors="raise").to_numpy(float)
    values = data[value_columns].apply(pd.to_numeric, errors="raise").to_numpy(float)
    if not np.isfinite(wavelengths).all() or not np.isfinite(values).all():
        raise ValueError(f"Non-finite spectral values found in {path}")
    return wavelengths, values


def load_spectrum(path, format_spec="spectragryph_tsv"):
    wavelengths, values = load_spectra(path, format_spec)
    if values.shape[1] != 1:
        raise ValueError(f"Expected one spectrum in {path}")
    return wavelengths, values[:, 0]


def load_timestamps(path, format_spec="ahk_csv"):
    spec = _format_spec(format_spec)
    if spec["type"] == "ahk_csv":
        data = pd.read_csv(path)
        event_column = _column(data, spec.get("event_column", "Event"))
        time_column = _column(data, spec.get("time_column", "ElapsedTime (s)"))
        events = data[event_column].astype(str)
        on = data.loc[events == "LEDon", time_column].to_numpy(float)
        off = data.loc[events == "LEDoff", time_column].to_numpy(float)
        measurements = int((events == "Measure").sum())
        cycles = min(len(on), len(off))
        timestamps = np.cumsum(np.r_[0, off[:cycles] - on[:cycles]])[:measurements]
    elif spec["type"] in {"generic_delimited", "simple_csv"}:
        if spec["type"] == "simple_csv":
            spec = {"type": "generic_delimited", "delimiter": ",", "header": True,
                    "time_column": 1}
        data = _read_delimited(path, spec)
        time_column = _column(data, spec.get("time_column", 1))
        timestamps = pd.to_numeric(data[time_column], errors="raise").to_numpy(float)
    else:
        raise ValueError(f"Unsupported timestamp format: {spec['type']}")
    if not len(timestamps) or not np.isfinite(timestamps).all():
        raise ValueError(f"No finite timestamps found in {path}")
    return timestamps


def _format_spec(value):
    return {"type": value} if isinstance(value, str) else dict(value)


def _read_delimited(path, spec):
    delimiter = spec.get("delimiter", ",")
    if delimiter == "tab":
        delimiter = "\t"
    header = 0 if spec.get("header", True) else None
    return pd.read_csv(Path(path), sep=delimiter, header=header,
                       skiprows=spec.get("skip_rows", 0),
                       decimal=spec.get("decimal", "."))


def _column(data, reference):
    if isinstance(reference, int):
        try:
            return data.columns[reference]
        except IndexError as error:
            raise ValueError(f"Column index {reference} is outside the input table") from error
    if reference not in data.columns:
        raise ValueError(f"Column {reference!r} was not found")
    return reference
