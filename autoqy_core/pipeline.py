"""Headless AutoQY analysis pipeline."""

from dataclasses import dataclass, replace
import warnings

import numpy as np

from .kinetics import (YieldFit, analyze_nipe_time_windows,
                       extrapolate_photostationary_state,
                       fit_quantum_yields, fit_quantum_yields_absorbance,
                       fit_quantum_yields_nipe, fit_quantum_yields_ode_absorbance)
from .spectra import (ABModelDiagnostic, ConcentrationFit, assess_ab_model,
                      fit_concentrations, fit_concentrations_regularized,
                      fit_concentrations_variable_total, interpolate_inputs,
                      process_led)


@dataclass(frozen=True)
class AnalysisInput:
    wavelengths: np.ndarray
    absorbance: np.ndarray
    timestamps: np.ndarray
    epsilon_r: tuple[np.ndarray, np.ndarray]
    epsilon_p: tuple[np.ndarray, np.ndarray]
    led: tuple[np.ndarray, np.ndarray]
    power_mw: float | None
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
    regularization_strength: float = 1
    absorbance_baseline_order: int = 1
    robust_loss_scale: float = 0.02
    initial_yields: tuple[float, float] = (0.5, 0.5)
    yield_bounds: tuple[float, float] = (0, 1)
    thermal_forward_rate: float = 0
    photon_flux_mol_s: float | None = None
    photon_flux_error_mol_s: float = 0
    irradiation_wavelength_nm: float | None = None


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
    extrapolated_pss: np.ndarray
    epsilon_uncertainty: object | None = None
    ab_model_diagnostic: ABModelDiagnostic | None = None


def run_analysis_pipeline(data):
    low = max(data.wavelength_limits[0], data.epsilon_r[0][0], data.epsilon_p[0][0])
    high = min(data.wavelength_limits[1], data.epsilon_r[0][-1], data.epsilon_p[0][-1])
    start = np.argmin(np.abs(data.wavelengths - low))
    stop = np.argmin(np.abs(data.wavelengths - high))
    wavelengths = data.wavelengths[start:stop]
    absorbance = data.absorbance[start:stop]
    actinometer_mode = data.photon_flux_mol_s is not None
    if actinometer_mode:
        led_processed = np.zeros_like(data.led[1], dtype=float)
        nominal_index = np.argmin(np.abs(data.led[0] - data.irradiation_wavelength_nm))
        led_processed[nominal_index] = 100.0
    else:
        led_processed = process_led(
            *data.led, data.baseline_correct_led, data.led_smoothing_window,
            data.led_polynomial_order, data.baseline_exclusion_fwhm_multiplier,
        )
    epsilon_r, epsilon_p, emission = interpolate_inputs(
        wavelengths, data.epsilon_r, data.epsilon_p, data.led[0], led_processed
    )
    if actinometer_mode:
        emission = np.zeros_like(wavelengths, dtype=float)
        emission[np.argmin(np.abs(wavelengths - data.irradiation_wavelength_nm))] = 100.0
    variable_concentration_fit = fit_concentrations_variable_total(
        absorbance, epsilon_r, epsilon_p, data.path_length_cm
    )
    ab_model_diagnostic = assess_ab_model(variable_concentration_fit, absorbance)
    if data.fit_method == "nipe":
        concentration_fit = variable_concentration_fit
    elif data.fit_method in {"regularized_concentrations", "ode_absorbance"}:
        concentration_fit = fit_concentrations_regularized(
            absorbance, wavelengths, epsilon_r, epsilon_p, data.timestamps,
            data.path_length_cm, data.regularization_strength,
        )
    else:
        concentration_fit = fit_concentrations(
            absorbance, wavelengths, epsilon_r, epsilon_p, data.path_length_cm
        )

    initial_total = float(np.sum(concentration_fit.concentrations[0]))
    initial_product_fraction = (float(concentration_fit.concentrations[0, 1]) /
                                initial_total if initial_total > 0 else 0.0)
    if initial_product_fraction > 0.02:
        warnings.warn(
            f"Initial product is {initial_product_fraction:.1%} of total concentration; "
            "verify spectral references, baseline, and species assignment. This can be "
            "physical, but the legacy emission method assumes zero initial product.",
            RuntimeWarning,
            stacklevel=2,
        )

    kinetic_slice = slice(None)
    if data.fit_method == "emission":
        threshold = emission.max() * data.emission_threshold_fraction
        active = np.flatnonzero(emission > threshold)
        if not len(active):
            raise ValueError("No LED-emission points exceed the configured threshold")
        if actinometer_mode:
            kinetic_slice = slice(max(active[0] - 1, 0), min(active[-1] + 2, len(emission)))
        else:
            kinetic_slice = slice(active[0], active[-1] + 1)

    if actinometer_mode:
        levels = (
            data.photon_flux_mol_s,
            data.photon_flux_mol_s + data.photon_flux_error_mol_s,
            data.photon_flux_mol_s - data.photon_flux_error_mol_s,
        )
    else:
        levels = (
            data.power_mw,
            data.power_mw + data.power_error_mw,
            data.power_mw - data.power_error_mw,
        )
    fits = []
    for level in levels:
        power = None if actinometer_mode else level
        photon_flux = level if actinometer_mode else None
        if data.fit_method in {"concentrations", "regularized_concentrations"}:
            fits.append(fit_quantum_yields(
                wavelengths, emission, concentration_fit.concentrations,
                data.timestamps, epsilon_r, epsilon_p, power, data.volume_ml,
                data.thermal_rate, data.path_length_cm, data.initial_yields,
                data.yield_bounds, data.thermal_forward_rate, photon_flux,
                data.irradiation_wavelength_nm,
            ))
        elif data.fit_method == "emission":
            fits.append(fit_quantum_yields_absorbance(
                wavelengths[kinetic_slice], emission[kinetic_slice],
                absorbance[kinetic_slice], data.timestamps,
                epsilon_r[kinetic_slice], epsilon_p[kinetic_slice], power,
                data.volume_ml, data.thermal_rate, data.path_length_cm,
                data.initial_yields, data.yield_bounds, data.thermal_forward_rate,
                photon_flux, data.irradiation_wavelength_nm,
            ))
        elif data.fit_method == "ode_absorbance":
            fits.append(fit_quantum_yields_ode_absorbance(
                wavelengths, emission, absorbance, data.timestamps, epsilon_r,
                epsilon_p, power, data.volume_ml, data.thermal_rate,
                data.path_length_cm, data.initial_yields, data.yield_bounds,
                concentration_fit.concentrations[0], data.absorbance_baseline_order,
                data.robust_loss_scale, data.thermal_forward_rate, photon_flux,
                data.irradiation_wavelength_nm,
            ))
        elif data.fit_method == "nipe":
            fits.append(fit_quantum_yields_nipe(
                wavelengths, emission, absorbance,
                variable_concentration_fit.concentrations, data.timestamps,
                epsilon_r, epsilon_p, power, data.volume_ml, data.thermal_rate,
                data.path_length_cm, data.initial_yields, data.yield_bounds,
                data.thermal_forward_rate, photon_flux,
                data.irradiation_wavelength_nm,
            ))
        else:
            raise ValueError(f"Unsupported fit method: {data.fit_method}")

    if data.fit_method == "nipe":
        window_analysis = analyze_nipe_time_windows(
            wavelengths, emission, absorbance,
            variable_concentration_fit.concentrations, data.timestamps,
            epsilon_r, epsilon_p, data.power_mw, data.volume_ml,
            data.thermal_rate, data.path_length_cm, data.initial_yields,
            data.yield_bounds, data.thermal_forward_rate,
            data.photon_flux_mol_s, data.irradiation_wavelength_nm,
            full_trace_values=fits[0].values,
        )
        fits[0] = replace(
            fits[0], nipe=replace(fits[0].nipe, window_analysis=window_analysis)
        )

    lower = fits[1].values - fits[1].standard_errors
    upper = fits[2].values + fits[2].standard_errors
    errors = np.maximum(fits[0].values - lower, upper - fits[0].values)
    extrapolated_pss = extrapolate_photostationary_state(
        wavelengths[kinetic_slice], emission[kinetic_slice],
        fits[0].concentrations[0].sum(), fits[0].values,
        epsilon_r[kinetic_slice], epsilon_p[kinetic_slice], data.power_mw,
        data.volume_ml, data.thermal_rate, data.path_length_cm,
        data.thermal_forward_rate, data.photon_flux_mol_s,
        data.irradiation_wavelength_nm,
    )
    return AnalysisResult(
        concentration_fit, fits[0], errors, wavelengths, absorbance,
        epsilon_r, epsilon_p, data.fit_method, extrapolated_pss,
        ab_model_diagnostic=ab_model_diagnostic,
    )


def run_concentration_analysis(data):
    """Backward-compatible name for callers of the first extracted core."""
    return run_analysis_pipeline(data)
