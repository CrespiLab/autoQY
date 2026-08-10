"""Headless AutoQY analysis pipeline."""

from dataclasses import dataclass

import numpy as np

from .kinetics import YieldFit, fit_quantum_yields, fit_quantum_yields_absorbance
from .spectra import ConcentrationFit, fit_concentrations, interpolate_inputs, process_led


@dataclass(frozen=True)
class AnalysisInput:
    wavelengths: np.ndarray
    absorbance: np.ndarray
    timestamps: np.ndarray
    epsilon_r: tuple[np.ndarray, np.ndarray]
    epsilon_p: tuple[np.ndarray, np.ndarray]
    led: tuple[np.ndarray, np.ndarray]
    power_mw: float
    power_error_mw: float
    volume_ml: float
    thermal_rate: float = 0
    path_length_cm: float = 1
    wavelength_limits: tuple[float, float] = (250, 800)
    baseline_correct_led: bool = True
    led_smoothing_window: int = 12
    led_polynomial_order: int = 3
    baseline_exclusion_fwhm_multiplier: float = 10
    fit_method: str = "concentrations"
    emission_threshold_fraction: float = 0.01
    initial_yields: tuple[float, float] = (0.5, 0.5)
    yield_bounds: tuple[float, float] = (0, 1)


@dataclass(frozen=True)
class AnalysisResult:
    concentration_fit: ConcentrationFit
    yield_fit: YieldFit
    yield_errors: np.ndarray
    wavelengths: np.ndarray
    absorbance: np.ndarray
    epsilon_r: np.ndarray
    epsilon_p: np.ndarray
    fit_method: str


def run_analysis_pipeline(data):
    low = max(data.wavelength_limits[0], data.epsilon_r[0][0], data.epsilon_p[0][0])
    high = min(data.wavelength_limits[1], data.epsilon_r[0][-1], data.epsilon_p[0][-1])
    start = np.argmin(np.abs(data.wavelengths - low))
    stop = np.argmin(np.abs(data.wavelengths - high))
    wavelengths = data.wavelengths[start:stop]
    absorbance = data.absorbance[start:stop]
    led_processed = process_led(
        *data.led, data.baseline_correct_led, data.led_smoothing_window,
        data.led_polynomial_order, data.baseline_exclusion_fwhm_multiplier,
    )
    epsilon_r, epsilon_p, emission = interpolate_inputs(
        wavelengths, data.epsilon_r, data.epsilon_p, data.led[0], led_processed
    )
    concentration_fit = fit_concentrations(
        absorbance, wavelengths, epsilon_r, epsilon_p, data.path_length_cm
    )

    fits = []
    for power in (data.power_mw, data.power_mw + data.power_error_mw,
                  data.power_mw - data.power_error_mw):
        if data.fit_method == "concentrations":
            fits.append(fit_quantum_yields(
                wavelengths, emission, concentration_fit.concentrations,
                data.timestamps, epsilon_r, epsilon_p, power, data.volume_ml,
                data.thermal_rate, data.path_length_cm, data.initial_yields,
                data.yield_bounds,
            ))
        elif data.fit_method == "emission":
            threshold = emission.max() * data.emission_threshold_fraction
            active = np.flatnonzero(emission > threshold)
            if not len(active):
                raise ValueError("No LED-emission points exceed the configured threshold")
            fit_slice = slice(active[0], active[-1] + 1)
            fits.append(fit_quantum_yields_absorbance(
                wavelengths[fit_slice], emission[fit_slice], absorbance[fit_slice],
                data.timestamps, epsilon_r[fit_slice], epsilon_p[fit_slice], power,
                data.volume_ml, data.thermal_rate, data.path_length_cm,
                data.initial_yields, data.yield_bounds,
            ))
        else:
            raise ValueError(f"Unsupported fit method: {data.fit_method}")

    lower = fits[1].values - fits[1].standard_errors
    upper = fits[2].values + fits[2].standard_errors
    errors = np.maximum(fits[0].values - lower, upper - fits[0].values)
    return AnalysisResult(
        concentration_fit, fits[0], errors, wavelengths, absorbance,
        epsilon_r, epsilon_p, data.fit_method,
    )


def run_concentration_analysis(data):
    """Backward-compatible name for callers of the first extracted core."""
    return run_analysis_pipeline(data)
