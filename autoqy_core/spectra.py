"""Spectral preprocessing and concentration recovery."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, nnls
from scipy.signal import find_peaks, peak_widths, savgol_filter


@dataclass(frozen=True)
class ConcentrationFit:
    concentrations: np.ndarray
    fractions: np.ndarray
    reconstructed_absorbance: np.ndarray
    residuals: np.ndarray


def process_led(wavelengths, intensity, baseline_correction=True, smoothing_window=12,
                polynomial_order=3, exclusion_fwhm_multiplier=10):
    processed = savgol_filter(np.asarray(intensity, float), smoothing_window, polynomial_order)
    if baseline_correction:
        peaks, _ = find_peaks(processed, height=processed.max() * 0.5)
        peak = peaks[np.argmax(processed[peaks])]
        fwhm = peak_widths(processed, [peak], rel_height=0.5)[0][0] * np.diff(wavelengths[:2])[0]
        mask = np.abs(wavelengths - wavelengths[peak]) > exclusion_fwhm_multiplier * fwhm
        processed -= np.mean(processed[mask])
    return np.maximum(processed, 0)


def interpolate_inputs(wavelengths, epsilon_r, epsilon_p, led_wavelengths, led_intensity):
    return (
        np.interp(wavelengths, epsilon_r[0], epsilon_r[1]),
        np.interp(wavelengths, epsilon_p[0], epsilon_p[1]),
        np.interp(wavelengths, led_wavelengths, led_intensity),
    )


def fit_concentrations(absorbance, wavelengths, epsilon_r, epsilon_p, path_length_cm=1):
    epsilon = np.column_stack((epsilon_r, epsilon_p)) * path_length_cm
    coefficients = np.array([nnls(epsilon, spectrum)[0] for spectrum in absorbance.T])
    totals = coefficients.sum(axis=1)
    fractions = np.divide(
        coefficients,
        totals[:, None],
        out=np.zeros_like(coefficients),
        where=totals[:, None] > 0,
    )
    total_concentration = np.trapezoid(absorbance[:, 0], wavelengths) / np.trapezoid(
        path_length_cm * (fractions[0, 0] * epsilon_r + fractions[0, 1] * epsilon_p),
        wavelengths,
    )
    concentrations = fractions * total_concentration
    reconstructed = concentrations @ epsilon.T
    return ConcentrationFit(concentrations, fractions, reconstructed, absorbance.T - reconstructed)


def fit_concentrations_regularized(absorbance, wavelengths, epsilon_r, epsilon_p,
                                   timestamps, path_length_cm=1,
                                   regularization_strength=1):
    """Recover a conserved concentration pair with a soft temporal trend.

    The exponential is only a regularizer: every timestamp retains its own fitted
    fraction. This prevents small spectral mismatches from forcing consecutive
    boundary values without making an exponential the kinetic model.
    """
    if regularization_strength <= 0:
        raise ValueError("regularization_strength must be positive")
    target = np.asarray(absorbance, float).T
    timestamps = np.asarray(timestamps, float)
    if target.shape[0] != len(timestamps):
        raise ValueError("One timestamp is required for each absorbance spectrum")

    independent = fit_concentrations(
        absorbance, wavelengths, epsilon_r, epsilon_p, path_length_cm
    )
    initial_fraction = independent.fractions[:, 0]
    concentration_scale = max(float(independent.concentrations[0].sum()),
                              np.finfo(float).eps)
    duration = np.ptp(timestamps)
    scaled_time = ((timestamps - timestamps[0]) / duration
                   if duration > 0 else np.zeros_like(timestamps))
    wavelength_count = target.shape[1]
    absorbance_scale = max(float(np.max(np.abs(target))), np.finfo(float).eps)

    epsilon_p_scaled = path_length_cm * epsilon_p
    epsilon_difference = path_length_cm * (epsilon_r - epsilon_p)
    difference_norm = np.dot(epsilon_difference, epsilon_difference)
    temporal_weight = regularization_strength * wavelength_count

    def evaluate(values, return_model=False):
        total = concentration_scale * np.exp(values[0])
        start, plateau, log_rate = values[1:]
        trend = plateau + (start - plateau) * np.exp(-np.exp(log_rate) * scaled_time)
        base = total * epsilon_p_scaled
        numerator = (total * ((target - base) @ epsilon_difference)
                     / absorbance_scale ** 2 + temporal_weight * trend)
        denominator = (total ** 2 * difference_norm / absorbance_scale ** 2
                       + temporal_weight)
        fraction_r = np.clip(numerator / denominator, 0, 1)
        model = base + total * fraction_r[:, None] * epsilon_difference
        if return_model:
            return fraction_r, total, model
        spectral = ((model - target) / absorbance_scale).ravel()
        temporal = (np.sqrt(regularization_strength * wavelength_count)
                    * (fraction_r - trend))
        return np.concatenate((spectral, temporal))

    initial = np.array([0, initial_fraction[0], initial_fraction[-1], 0])
    lower = np.array([np.log(0.1), 0, 0, -6])
    upper = np.array([np.log(10), 1, 1, 6])
    fit = least_squares(evaluate, initial, bounds=(lower, upper))
    fraction_r, total, reconstructed = evaluate(fit.x, return_model=True)
    fractions = np.column_stack((fraction_r, 1 - fraction_r))
    concentrations = total * fractions
    return ConcentrationFit(
        concentrations, fractions, reconstructed, target - reconstructed
    )
