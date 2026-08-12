"""Molar-absorptivity calculations for independently prepared spectra."""

from dataclasses import dataclass
from io import StringIO

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class EpsilonResult:
    wavelengths: np.ndarray
    absorbance: np.ndarray
    concentrations_m: np.ndarray
    path_lengths_cm: np.ndarray
    individual: np.ndarray
    mean: np.ndarray
    standard_deviation: np.ndarray
    standard_error: np.ndarray


@dataclass(frozen=True)
class NMRSubtractionResult:
    wavelengths: np.ndarray
    normalized_reactant: np.ndarray
    normalized_pss: np.ndarray
    reconstructed_reactant: np.ndarray
    reconstructed_pss: np.ndarray
    product: np.ndarray
    product_lower: np.ndarray
    product_upper: np.ndarray
    product_error_lower: np.ndarray
    product_error_upper: np.ndarray
    negative_product_points: int
    negative_bound_points: int


def concentration_from_preparation(mass_mg, volume_ml, molecular_weight_g_mol):
    """Return mol/L from solute mass, final solution volume, and molecular weight."""
    mass = _positive_float(mass_mg, "Mass")
    volume = _positive_float(volume_ml, "Solution volume")
    molecular_weight = _positive_float(molecular_weight_g_mol, "Molecular weight")
    return mass / (volume * molecular_weight)


def calculate_epsilon(wavelengths, absorbance, concentrations_m, path_lengths_cm):
    """Apply Beer-Lambert independently to each wavelength x replicate column."""
    wavelengths = np.asarray(wavelengths, float)
    absorbance = np.asarray(absorbance, float)
    concentrations = np.asarray(concentrations_m, float)
    path_lengths = np.asarray(path_lengths_cm, float)
    if wavelengths.ndim != 1 or absorbance.ndim != 2:
        raise ValueError("Wavelengths and absorbance must be one- and two-dimensional")
    if absorbance.shape[0] != len(wavelengths):
        raise ValueError("Absorbance rows must match the wavelength axis")
    if concentrations.shape != (absorbance.shape[1],):
        raise ValueError("Provide one concentration for every spectrum")
    if path_lengths.shape != (absorbance.shape[1],):
        raise ValueError("Provide one path length for every spectrum")
    if not all(np.isfinite(item).all() for item in
               (wavelengths, absorbance, concentrations, path_lengths)):
        raise ValueError("Epsilon inputs must contain only finite values")
    if np.any(concentrations <= 0) or np.any(path_lengths <= 0):
        raise ValueError("Concentrations and path lengths must be positive")
    return absorbance / (concentrations * path_lengths)[None, :]


def calculate_epsilon_statistics(wavelengths, absorbance, concentrations_m,
                                 path_lengths_cm):
    """Calculate replicate epsilon spectra and wavelength-resolved sample errors."""
    individual = calculate_epsilon(
        wavelengths, absorbance, concentrations_m, path_lengths_cm
    )
    count = individual.shape[1]
    if count > 1:
        standard_deviation = np.std(individual, axis=1, ddof=1)
        standard_error = standard_deviation / np.sqrt(count)
    else:
        standard_deviation = np.full(individual.shape[0], np.nan)
        standard_error = np.full(individual.shape[0], np.nan)
    return EpsilonResult(
        np.asarray(wavelengths, float), np.asarray(absorbance, float),
        np.asarray(concentrations_m, float), np.asarray(path_lengths_cm, float),
        individual, np.mean(individual, axis=1), standard_deviation,
        standard_error,
    )


def export_epsilon_tsv(result, labels):
    """Export processed absorbance, epsilon replicates, and wavelength errors."""
    labels = _unique_labels(labels)
    if len(labels) != result.individual.shape[1]:
        raise ValueError("Spectrum labels do not match the epsilon matrix")
    columns = {"Wavelength_nm": result.wavelengths}
    for index, label in enumerate(labels):
        columns[f"Absorbance__{label}"] = result.absorbance[:, index]
        columns[f"Concentration_M__{label}"] = np.full(
            len(result.wavelengths), result.concentrations_m[index]
        )
        columns[f"Path_length_cm__{label}"] = np.full(
            len(result.wavelengths), result.path_lengths_cm[index]
        )
        columns[f"Epsilon_M-1_cm-1__{label}"] = result.individual[:, index]
    columns["Epsilon_mean_M-1_cm-1"] = result.mean
    columns["Epsilon_SD_M-1_cm-1"] = result.standard_deviation
    columns["Epsilon_SEM_M-1_cm-1"] = result.standard_error
    lower, upper = nonnegative_error_bounds(result.mean, result.standard_deviation)
    columns["Epsilon_lower_nonnegative_M-1_cm-1"] = lower
    columns["Epsilon_upper_nonnegative_M-1_cm-1"] = upper
    stream = StringIO()
    pd.DataFrame(columns).to_csv(stream, sep="\t", index=False, float_format="%.8e")
    return stream.getvalue()


def load_epsilon_tsv(text):
    """Reload an AutoQY epsilon TSV, including its processed absorbance metadata."""
    frame = pd.read_csv(StringIO(text), sep="\t")
    required = {"Wavelength_nm", "Epsilon_mean_M-1_cm-1",
                "Epsilon_SD_M-1_cm-1", "Epsilon_SEM_M-1_cm-1"}
    if not required.issubset(frame.columns):
        raise ValueError("This is not an AutoQY epsilon TSV export")
    prefixes = ("Absorbance__", "Concentration_M__", "Path_length_cm__",
                "Epsilon_M-1_cm-1__")
    labels = [column[len(prefixes[0]):] for column in frame.columns
              if column.startswith(prefixes[0])]
    if not labels:
        raise ValueError("AutoQY epsilon TSV contains no replicate absorbance columns")
    absorbance, concentrations, paths, individual = [], [], [], []
    for label in labels:
        names = [prefix + label for prefix in prefixes]
        if not all(name in frame.columns for name in names):
            raise ValueError(f"AutoQY epsilon metadata is incomplete for {label}")
        absorbance.append(pd.to_numeric(frame[names[0]], errors="raise").to_numpy(float))
        concentration_values = pd.to_numeric(frame[names[1]], errors="raise").to_numpy(float)
        path_values = pd.to_numeric(frame[names[2]], errors="raise").to_numpy(float)
        if not np.allclose(concentration_values, concentration_values[0]):
            raise ValueError(f"Concentration metadata varies with wavelength for {label}")
        if not np.allclose(path_values, path_values[0]):
            raise ValueError(f"Path-length metadata varies with wavelength for {label}")
        concentrations.append(concentration_values[0])
        paths.append(path_values[0])
        individual.append(pd.to_numeric(frame[names[3]], errors="raise").to_numpy(float))
    result = EpsilonResult(
        pd.to_numeric(frame["Wavelength_nm"], errors="raise").to_numpy(float),
        np.column_stack(absorbance), np.asarray(concentrations), np.asarray(paths),
        np.column_stack(individual),
        pd.to_numeric(frame["Epsilon_mean_M-1_cm-1"], errors="raise").to_numpy(float),
        pd.to_numeric(frame["Epsilon_SD_M-1_cm-1"], errors="raise").to_numpy(float),
        pd.to_numeric(frame["Epsilon_SEM_M-1_cm-1"], errors="raise").to_numpy(float),
    )
    arrays = (result.wavelengths, result.absorbance, result.concentrations_m,
              result.path_lengths_cm, result.individual, result.mean)
    if not all(np.isfinite(array).all() for array in arrays):
        raise ValueError("AutoQY epsilon TSV contains non-finite required values")
    return result, labels


def nonnegative_error_bounds(mean, error):
    """Return mean ± non-negative error constrained to physical epsilon bounds."""
    mean = np.asarray(mean, float)
    error = np.abs(np.asarray(error, float))
    return np.maximum(mean - error, 0.0), np.maximum(mean + error, 0.0)


def reconstruct_product_from_nmr(wavelengths, reactant_absorbance, pss_absorbance,
                                 reactant_epsilon, reactant_percent,
                                 nmr_error_percent=1.0):
    """Infer product epsilon from a PSS spectrum and its NMR composition.

    Both absorbance spectra are divided by the reactant absorbance maximum. The
    PSS trace is then placed on the reactant-epsilon scale before the mixture is
    resolved as epsilon_PSS = f_R epsilon_R + (1-f_R) epsilon_P.
    """
    wavelengths = np.asarray(wavelengths, float)
    reactant_absorbance = np.asarray(reactant_absorbance, float)
    pss_absorbance = np.asarray(pss_absorbance, float)
    reactant_sd = np.nan_to_num(
        np.asarray(reactant_epsilon.standard_deviation, float), nan=0.0,
        posinf=0.0, neginf=0.0,
    )
    if any(array.shape != wavelengths.shape for array in
           (reactant_absorbance, pss_absorbance, reactant_epsilon.mean,
            reactant_sd)):
        raise ValueError("NMR spectra and reactant epsilon must share one wavelength axis")
    if not all(np.isfinite(array).all() for array in
               (wavelengths, reactant_absorbance, pss_absorbance,
                reactant_epsilon.mean)):
        raise ValueError("NMR reconstruction inputs must be finite")
    peak = float(np.max(reactant_absorbance))
    if peak <= 0:
        raise ValueError("The reactant absorbance maximum must be positive")
    fraction = float(reactant_percent) / 100.0
    fraction_error = abs(float(nmr_error_percent)) / 100.0
    if not 0 <= fraction < 1:
        raise ValueError("Reactant percentage in the PSS must be from 0 up to 100%")
    if not np.isfinite(fraction_error) or fraction_error < 0:
        raise ValueError("NMR percentage error must be finite and non-negative")
    fractions = np.unique(np.clip(
        [fraction - fraction_error, fraction, fraction + fraction_error],
        0.0, 1.0 - np.finfo(float).eps,
    ))
    normalized_reactant = reactant_absorbance / peak
    normalized_pss = pss_absorbance / peak
    epsilon_lower, epsilon_upper = nonnegative_error_bounds(
        reactant_epsilon.mean, reactant_sd
    )
    epsilon_scenarios = (epsilon_lower, reactant_epsilon.mean, epsilon_upper)

    def reconstruct(epsilon_curve, selected_fraction):
        scale = float(np.max(epsilon_curve))
        mixture = normalized_pss * scale
        product = (mixture - selected_fraction * epsilon_curve) / (1 - selected_fraction)
        return mixture, product

    reconstructed_pss, nominal = reconstruct(reactant_epsilon.mean, fraction)
    candidates = [reconstruct(curve, selected_fraction)[1]
                  for curve in epsilon_scenarios for selected_fraction in fractions]
    raw_lower = np.min(candidates, axis=0)
    raw_upper = np.max(candidates, axis=0)
    lower = np.maximum(raw_lower, 0.0)
    upper = np.maximum(raw_upper, 0.0)
    error_lower = np.abs(nominal - lower)
    error_upper = np.abs(upper - nominal)
    scale = float(np.max(reactant_epsilon.mean))
    return NMRSubtractionResult(
        wavelengths, normalized_reactant, normalized_pss,
        reactant_epsilon.mean, reconstructed_pss, nominal,
        lower, upper, error_lower, error_upper,
        int(np.count_nonzero(nominal < 0)),
        int(np.count_nonzero((raw_lower < 0) | (raw_upper < 0))),
    )


def export_nmr_subtraction_tsv(result, preserve_negative=False):
    """Export the NMR-derived epsilon, always retaining a raw audit column."""
    exported_product = (result.product if preserve_negative
                        else np.maximum(result.product, 0.0))
    frame = pd.DataFrame({
        "Wavelength_nm": result.wavelengths,
        "Reactant_normalized": result.normalized_reactant,
        "PSS_normalized_to_reactant": result.normalized_pss,
        "Reactant_reconstructed_M-1_cm-1": result.reconstructed_reactant,
        "PSS_reconstructed_M-1_cm-1": result.reconstructed_pss,
        "Product_epsilon_M-1_cm-1": exported_product,
        "Product_epsilon_raw_M-1_cm-1": result.product,
        "Product_lower_nonnegative_M-1_cm-1": result.product_lower,
        "Product_upper_nonnegative_M-1_cm-1": result.product_upper,
        "Product_error_lower_M-1_cm-1": result.product_error_lower,
        "Product_error_upper_M-1_cm-1": result.product_error_upper,
    })
    return frame.to_csv(sep="\t", index=False, float_format="%.8e")


def _positive_float(value, name):
    parsed = float(value)
    if not np.isfinite(parsed) or parsed <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return parsed


def _unique_labels(labels):
    used = set()
    result = []
    for position, value in enumerate(labels, 1):
        label = "".join(character if character.isalnum() or character in "-_" else "_"
                        for character in str(value)).strip("_") or f"spectrum_{position}"
        candidate = label
        number = 2
        while candidate in used:
            candidate = f"{label}_{number}"
            number += 1
        used.add(candidate)
        result.append(candidate)
    return result
