"""SVD denoising for time-resolved spectral datasets."""

from dataclasses import dataclass
from io import StringIO
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.sparse import diags, eye
from scipy.sparse.linalg import factorized

from .io import load_avantes_abs8_bytes


@dataclass(frozen=True)
class SpectralDataset:
    wavelengths: np.ndarray
    coordinates: np.ndarray
    absorbance: np.ndarray  # wavelength x measurement
    source_format: str = "spectragryph"
    interpolated_values: int = 0


@dataclass(frozen=True)
class SVDResult:
    processed: np.ndarray
    left_vectors: np.ndarray
    singular_values: np.ndarray
    right_vectors: np.ndarray
    explained_weights: np.ndarray
    cumulative_weights: np.ndarray
    proposed_rank: int

    def reconstruct(self, rank):
        rank = validate_rank(rank, len(self.singular_values))
        return (self.left_vectors[:, :rank] * self.singular_values[:rank]) @ self.right_vectors[:rank]


def load_spectral_dataset(path, format_name="auto"):
    path = Path(path)
    if path.suffix.lower() == ".abs8":
        format_name = "avantes_abs8"
    return load_spectral_bytes(path.read_bytes(), format_name)


def load_spectral_bytes(data, format_name="auto"):
    """Load either an Avantes binary record or a supported text dataset."""
    if format_name == "avantes_abs8":
        wavelengths, absorbance = load_avantes_abs8_bytes(data)
        absorbance, interpolated = _fill_missing_absorbance(
            wavelengths, absorbance[:, None]
        )
        return SpectralDataset(
            wavelengths, np.array([0.0]), absorbance, format_name, interpolated
        )
    try:
        text = bytes(data).decode("utf-8-sig")
    except UnicodeDecodeError as error:
        if format_name == "auto":
            try:
                return load_spectral_bytes(data, "avantes_abs8")
            except Exception as avantes_error:
                raise ValueError(
                    "Input is neither supported UTF-8 spectral text nor an Avantes Abs8 file"
                ) from avantes_error
        raise ValueError("Input is not UTF-8 spectral text") from error
    return load_spectral_text(text, format_name)


def load_spectral_text(text, format_name="auto"):
    """Load SpectraGryph matrix, TSV, or CSV text into wavelength x measurement form."""
    if format_name == "auto":
        errors = []
        for candidate in ("spectragryph", "tsv", "csv"):
            try:
                return load_spectral_text(text, candidate)
            except (TypeError, ValueError) as error:
                errors.append(f"{candidate}: {error}")
        raise ValueError("Could not detect spectral text format (" + "; ".join(errors) + ")")
    if format_name == "spectragryph":
        try:
            table = np.loadtxt(StringIO(text))
        except ValueError:
            wavelengths, coordinates, absorbance = _load_column_table(text, "\t")
        else:
            if table.ndim != 2 or min(table.shape) < 2:
                raise ValueError("SpectraGryph data must be a rectangular numeric matrix")
            wavelengths, coordinates, absorbance = table[0, 1:], table[1:, 0], table[1:, 1:].T
    elif format_name in {"tsv", "csv"}:
        wavelengths, coordinates, absorbance = _load_column_table(
            text, "\t" if format_name == "tsv" else ","
        )
    else:
        raise ValueError(f"Unsupported smoother format: {format_name}")
    _validate_axes(wavelengths, coordinates, absorbance)
    absorbance, interpolated = _fill_missing_absorbance(wavelengths, absorbance)
    return SpectralDataset(wavelengths, coordinates, absorbance, format_name, interpolated)


def baseline_spectra(dataset, interval=None):
    """Zero each spectrum by its mean inside an inclusive wavelength interval."""
    values = np.asarray(dataset.absorbance, float)
    if interval is None:
        return values.copy()
    low, high = map(float, interval)
    if low >= high:
        raise ValueError("Baseline interval start must be below its end")
    mask = (dataset.wavelengths >= low) & (dataset.wavelengths <= high)
    if not np.any(mask):
        raise ValueError("Baseline interval contains no measured wavelengths")
    return values - np.mean(values[mask], axis=0, keepdims=True)


def select_wavelengths(dataset, interval=None):
    """Return a wavelength-restricted dataset, preserving measurement columns."""
    if interval is None:
        return dataset
    low, high = map(float, interval)
    if low >= high:
        raise ValueError("Wavelength range start must be below its end")
    mask = (dataset.wavelengths >= low) & (dataset.wavelengths <= high)
    if np.count_nonzero(mask) < 2:
        raise ValueError("Wavelength range must contain at least two measured wavelengths")
    return SpectralDataset(dataset.wavelengths[mask], dataset.coordinates,
                           dataset.absorbance[mask], dataset.source_format,
                           dataset.interpolated_values)


def analyze_svd(dataset, baseline_interval=None):
    processed = baseline_spectra(dataset, baseline_interval)
    left, singular, right = np.linalg.svd(processed, full_matrices=False)
    energy = singular ** 2
    weights = energy / energy.sum() if energy.sum() else np.zeros_like(energy)
    rank = _weight_rank(weights)
    return SVDResult(processed, left, singular, right, weights, np.cumsum(weights), rank)


def smooth_reconstruction(values, method="off", *, savgol_window=11,
                          savgol_order=3, whittaker_strength=1000.0,
                          whittaker_order=2):
    """Optionally smooth a wavelength x measurement SVD reconstruction."""
    values = np.asarray(values, float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Reconstruction must be a finite two-dimensional matrix")
    if method == "off":
        return values.copy()
    if method == "savgol":
        window = _positive_integer(savgol_window, "Savitzky-Golay window")
        order = _nonnegative_integer(savgol_order, "Savitzky-Golay polynomial order")
        if window % 2 == 0:
            raise ValueError("Savitzky-Golay window must be odd")
        if window > values.shape[0]:
            raise ValueError("Savitzky-Golay window exceeds the selected wavelengths")
        if order >= window:
            raise ValueError("Savitzky-Golay polynomial order must be below its window")
        return savgol_filter(values, window, order, axis=0, mode="interp")
    if method == "whittaker":
        strength = float(whittaker_strength)
        order = _positive_integer(whittaker_order, "Whittaker difference order")
        if not np.isfinite(strength) or strength <= 0:
            raise ValueError("Whittaker strength must be positive and finite")
        if order >= values.shape[0]:
            raise ValueError("Whittaker difference order exceeds the selected wavelengths")
        coefficients = np.array([(-1) ** (order - index) * comb(order, index)
                                 for index in range(order + 1)], float)
        difference = diags(coefficients, np.arange(order + 1),
                           shape=(values.shape[0] - order, values.shape[0]), format="csc")
        solve = factorized(eye(values.shape[0], format="csc")
                           + strength * (difference.T @ difference))
        return np.column_stack([solve(values[:, index]) for index in range(values.shape[1])])
    raise ValueError(f"Unsupported spectral smoothing method: {method}")


def savgol_window_points(wavelengths, width_nm, polynomial_order=3):
    """Convert a physical Savitzky-Golay width to a valid odd point count."""
    wavelengths = np.asarray(wavelengths, float)
    width = float(width_nm)
    order = _nonnegative_integer(polynomial_order, "Savitzky-Golay polynomial order")
    if wavelengths.ndim != 1 or len(wavelengths) < 3 or not np.isfinite(wavelengths).all():
        raise ValueError("Savitzky-Golay wavelengths must contain at least three finite values")
    spacing = float(np.median(np.abs(np.diff(wavelengths))))
    if not np.isfinite(width) or width <= 0 or not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("Savitzky-Golay width and wavelength spacing must be positive")
    points = max(int(round(width / spacing)), order + 2, 3)
    if points % 2 == 0:
        points += 1
    if points > len(wavelengths):
        points = len(wavelengths) if len(wavelengths) % 2 else len(wavelengths) - 1
    if points <= order:
        raise ValueError("Selected wavelengths are too few for the Savitzky-Golay settings")
    return points


def export_smoothed_text(dataset, reconstructed, format_name=None):
    """Serialize a reconstruction as CSV unless another format is requested."""
    format_name = format_name or "csv"
    reconstructed = np.asarray(reconstructed, float)
    if reconstructed.shape != dataset.absorbance.shape:
        raise ValueError("Reconstructed matrix shape does not match the input")
    if format_name == "spectragryph":
        table = np.block([[np.array([[0.0]]), dataset.wavelengths[None, :]],
                          [dataset.coordinates[:, None], reconstructed.T]])
        stream = StringIO()
        np.savetxt(stream, table, fmt="%.8e", delimiter="\t")
        return stream.getvalue()
    if format_name in {"tsv", "avantes_abs8"}:
        separator = "\t"
    elif format_name == "csv":
        separator = ","
    else:
        raise ValueError(f"Unsupported smoother export format: {format_name}")
    value_columns = (["Absorbance"] if format_name == "avantes_abs8"
                     else [f"{value:g}" for value in dataset.coordinates])
    columns = ["Wavelength"] + value_columns
    frame = pd.DataFrame(np.column_stack([dataset.wavelengths, reconstructed]), columns=columns)
    return frame.to_csv(
        index=False, sep=separator, float_format="%.8e", lineterminator="\n"
    )


def validate_rank(rank, maximum):
    if not isinstance(rank, (int, np.integer)) or not 1 <= rank <= maximum:
        raise ValueError(f"SVD rank must be between 1 and {maximum}")
    return int(rank)


def _weight_rank(weights, target=0.995):
    """Suggest the smallest rank retaining a declared squared-SV weight."""
    return min(int(np.searchsorted(np.cumsum(weights), target) + 1), len(weights))


def _numeric_labels(labels):
    parsed = pd.to_numeric(pd.Index(labels), errors="coerce").to_numpy(float)
    return parsed if np.isfinite(parsed).all() else np.arange(len(labels), dtype=float)


def _load_column_table(text, separator):
    first_row = pd.read_csv(
        StringIO(text), sep=separator, header=None, nrows=1, dtype=str
    )
    populated = first_row.iloc[0].dropna()
    has_header = any(not _numeric_text(value) for value in populated)
    frame = pd.read_csv(
        StringIO(text), sep=separator, header=0 if has_header else None
    )
    # SpectraGryph exports can retain many trailing delimiters after the last
    # spectrum. They are empty padding columns, not measured spectra.
    frame = frame.dropna(axis=1, how="all")
    frame = frame.drop(columns="Wavenumbers [1/cm]", errors="ignore")
    if frame.shape[1] < 2:
        raise ValueError("Delimited data requires wavelength and at least one spectrum column")
    wavelengths = pd.to_numeric(frame.iloc[:, 0], errors="raise").to_numpy(float)
    absorbance = frame.iloc[:, 1:].apply(pd.to_numeric, errors="raise").to_numpy(float)
    coordinates = (_numeric_labels(frame.columns[1:]) if has_header
                   else np.arange(frame.shape[1] - 1, dtype=float))
    return wavelengths, coordinates, absorbance


def _numeric_text(value):
    try:
        float(str(value).strip())
    except ValueError:
        return False
    return True


def _validate_axes(wavelengths, coordinates, absorbance):
    if wavelengths.ndim != 1 or coordinates.ndim != 1 or absorbance.ndim != 2:
        raise ValueError("Invalid spectral dataset dimensions")
    if absorbance.shape != (len(wavelengths), len(coordinates)):
        raise ValueError("Spectral axes do not match the absorbance matrix")
    if not len(wavelengths) or not len(coordinates):
        raise ValueError("Spectral dataset is empty")
    if not np.isfinite(wavelengths).all() or not np.isfinite(coordinates).all():
        raise ValueError("Spectral axes contain non-finite values")


def _fill_missing_absorbance(wavelengths, absorbance):
    missing = ~np.isfinite(absorbance)
    count = int(np.count_nonzero(missing))
    if not count:
        return absorbance, 0
    filled = absorbance.copy()
    for column in np.flatnonzero(np.any(missing, axis=0)):
        valid = np.isfinite(filled[:, column])
        if np.count_nonzero(valid) < 2:
            raise ValueError("A spectrum contains fewer than two finite absorbance values")
        positions = np.arange(len(wavelengths), dtype=float)
        filled[~valid, column] = np.interp(
            positions[~valid], positions[valid], filled[valid, column]
        )
    return filled, count


def _positive_integer(value, name):
    parsed = _nonnegative_integer(value, name)
    if parsed == 0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _nonnegative_integer(value, name):
    if isinstance(value, bool) or int(value) != float(value) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)
