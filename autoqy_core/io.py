"""Load AutoQY input files without GUI dependencies."""

from pathlib import Path
import re
import struct

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


def load_cary(path):
    """Load spectra from a Varian/Agilent Cary WinUV DSW or BSW file."""
    return load_cary_bytes(Path(path).read_bytes())


def load_cary_bytes(data):
    """Decode wavelength/signal streams from a Cary WinUV binary container."""
    data = bytes(data)
    if len(data) < 18:
        raise ValueError("Cary file is too small")
    magic_length = data[0]
    magic = data[1:1 + magic_length].decode("ascii", errors="ignore")
    if not magic.startswith("Varian UV-VIS"):
        raise ValueError("Not a Varian/Agilent Cary WinUV DSW or BSW file")

    blocks = []
    offset = 0
    while offset <= len(data) - 16:
        first_wavelength, first_value, second_wavelength, second_value = struct.unpack_from(
            "<ffff", data, offset
        )
        step = second_wavelength - first_wavelength
        if (_cary_point(first_wavelength, first_value)
                and _cary_point(second_wavelength, second_value)
                and 0.001 <= abs(step) <= 30.0):
            decreasing = step < 0
            wavelengths = [first_wavelength, second_wavelength]
            values = [first_value, second_value]
            end = offset + 16
            while end + 8 <= len(data):
                wavelength, value = struct.unpack_from("<ff", data, end)
                next_step = wavelength - wavelengths[-1]
                monotonic = next_step < -0.001 if decreasing else next_step > 0.001
                if (not _cary_point(wavelength, value) or not monotonic
                        or abs(next_step) > 30.0):
                    break
                wavelengths.append(wavelength)
                values.append(value)
                end += 8
            if len(wavelengths) >= 30:
                blocks.append((offset, end, np.asarray(wavelengths, float),
                               np.asarray(values, float)))
                offset = end
                continue
        offset += 1

    unique = []
    for block in blocks:
        if unique and block[0] < unique[-1][1]:
            if len(block[2]) > len(unique[-1][2]):
                unique[-1] = block
        else:
            unique.append(block)
    if not unique:
        raise ValueError("No Cary spectral data streams were found")
    spectra = [(w[::-1], y[::-1]) if w[0] > w[-1] else (w, y)
               for _, _, w, y in unique]
    return _align_spectra(spectra, "Cary")


def load_specord(path):
    """Load spectra from an Analytik Jena SPECORD WinASPECT DAT file."""
    return load_specord_bytes(Path(path).read_bytes())


def load_specord_bytes(data):
    """Decode labelled float32 arrays from a SPECORD WinASPECT DAT file."""
    data = bytes(data)
    header_end = data.find(b"[DATA]")
    if header_end < 0 or not data.startswith(b"[GENERAL]"):
        raise ValueError("Not an Analytik Jena SPECORD WinASPECT DAT file")
    header = data[:header_end].decode("latin-1")
    if "ORIGIN=SPECORD" not in header.upper():
        raise ValueError("SPECORD origin marker is missing")
    point_match = re.search(r"(?mi)^NPOINTS=(\d+)\s*$", header)
    cycle_match = re.search(r"(?mi)^NCYCL=(\d+)\s*$", header)
    unit_match = re.search(r"(?mi)^XUNITS=([^\r\n]+)", header)
    if not point_match:
        raise ValueError("SPECORD point count is missing")
    points = int(point_match.group(1))
    cycles = int(cycle_match.group(1)) if cycle_match else 1
    if not 2 <= points <= 1_000_000 or not 1 <= cycles <= 100_000:
        raise ValueError("SPECORD point or cycle count is invalid")
    if unit_match and unit_match.group(1).strip().lower() not in {"nm", "nanometer", "nanometers"}:
        raise ValueError("SPECORD x axis must be wavelength in nm")

    wavelength_offset = _binary_marker_offset(data, b"XDATA=")
    value_offset = _binary_marker_offset(data, b"YDATA=")
    wavelength_bytes = points * np.dtype("<f4").itemsize
    value_bytes = points * cycles * np.dtype("<f4").itemsize
    if wavelength_offset + wavelength_bytes > len(data) or value_offset + value_bytes > len(data):
        raise ValueError("SPECORD data arrays are truncated")
    wavelengths = np.frombuffer(data, dtype="<f4", count=points,
                                offset=wavelength_offset).astype(float)
    values = np.frombuffer(data, dtype="<f4", count=points * cycles,
                           offset=value_offset).reshape(cycles, points).T.astype(float)
    if not np.isfinite(wavelengths).all() or not np.isfinite(values).all():
        raise ValueError("SPECORD file contains non-finite spectral values")
    differences = np.diff(wavelengths)
    if np.all(differences < 0):
        wavelengths, values = wavelengths[::-1], values[::-1]
    elif not np.all(differences > 0):
        raise ValueError("SPECORD wavelengths must be monotonic")
    return wavelengths, values


def load_spectra(path, format_spec="spectragryph_tsv"):
    spec = _format_spec(format_spec)
    if spec["type"] == "agilent_cary":
        return load_cary(path)
    if spec["type"] == "specord":
        return load_specord(path)
    if spec["type"] == "spectragryph_tsv":
        data = pd.read_csv(path, sep="\t", float_precision="round_trip").drop(
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
        data = pd.read_csv(path, float_precision="round_trip")
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
        default_column = 0 if len(data.columns) == 1 else 1
        time_column = _column(data, spec.get("time_column", default_column))
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
                       decimal=spec.get("decimal", "."),
                       float_precision="round_trip")


def _column(data, reference):
    if isinstance(reference, int):
        try:
            return data.columns[reference]
        except IndexError as error:
            raise ValueError(f"Column index {reference} is outside the input table") from error
    if reference not in data.columns:
        raise ValueError(f"Column {reference!r} was not found")
    return reference


def _cary_point(wavelength, value):
    return (np.isfinite(wavelength) and np.isfinite(value)
            and 100.0 <= wavelength <= 5000.0
            and -1_000_000.0 <= value <= 1_000_000.0)


def _align_spectra(spectra, name):
    low = max(wavelengths[0] for wavelengths, _ in spectra)
    high = min(wavelengths[-1] for wavelengths, _ in spectra)
    if high <= low:
        raise ValueError(f"{name} spectra have no common wavelength range")
    reference = min(spectra, key=lambda item: np.median(np.diff(item[0])))[0]
    common = reference[(reference >= low) & (reference <= high)]
    if len(common) < 2:
        raise ValueError(f"{name} spectra have fewer than two common wavelengths")
    values = np.column_stack([
        np.interp(common, wavelengths, signal) for wavelengths, signal in spectra
    ])
    return common, values


def _binary_marker_offset(data, marker):
    position = data.find(marker)
    if position < 0:
        raise ValueError(f"SPECORD {marker.decode('ascii').rstrip('=')} array is missing")
    position += len(marker)
    if data[position:position + 2] == b"\r\n":
        position += 2
    elif data[position:position + 1] == b"\n":
        position += 1
    return position
