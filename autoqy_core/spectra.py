"""Spectral preprocessing and concentration recovery."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import nnls
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
