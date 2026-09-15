"""Quantum-yield fitting independent of the user interface."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq, least_squares

H = 6.626e-34
C = 299792458
AVOGADRO = 6.022e23


@dataclass(frozen=True)
class NIPEWindowAnalysis:
    """Automatic short-window analysis before the concentration plateau."""

    plateau_detected: bool
    plateau_time_s: float
    analysis_end_time_s: float
    window_point_count: int
    window_duration_s: float
    window_start_s: np.ndarray
    window_end_s: np.ndarray
    window_midpoint_s: np.ndarray
    window_values: np.ndarray
    window_standard_errors: np.ndarray
    window_jacobian_conditions: np.ndarray
    extrapolated_values: np.ndarray
    extrapolated_standard_errors: np.ndarray
    full_trace_values: np.ndarray
    full_trace_change_fraction: np.ndarray


@dataclass(frozen=True)
class NIPEMetadata:
    """Audit information for a normalized integrated photokinetic fit."""

    interval_count: int
    minimum_interval_fraction: float
    normalized_residual_rmse: float
    window_analysis: NIPEWindowAnalysis | None = None


@dataclass(frozen=True)
class YieldFit:
    values: np.ndarray
    standard_errors: np.ndarray
    concentrations: np.ndarray
    absorbance_correction: np.ndarray | None = None
    optimizer_success: bool = True
    optimizer_message: str = ""
    jacobian_condition: float = np.nan
    active_bounds: tuple[bool, bool] = (False, False)
    nipe: NIPEMetadata | None = None


def fit_quantum_yields(wavelengths_nm, emission, concentrations, timestamps,
                       epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
                       path_length_cm=1, initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
                       thermal_forward_rate=0, photon_flux_mol_s=None,
                       irradiation_wavelength_nm=None):
    """Fit quantum yields to experimentally recovered concentration traces."""
    scale = max(float(np.max(np.abs(concentrations))), np.finfo(float).eps)

    def residual(model):
        return ((model - concentrations) / scale).ravel()

    return _fit(wavelengths_nm, emission, concentrations[0], timestamps,
                epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
                path_length_cm, initial_yields, yield_bounds, residual,
                thermal_forward_rate, photon_flux_mol_s,
                irradiation_wavelength_nm)


def fit_quantum_yields_absorbance(wavelengths_nm, emission, absorbance, timestamps,
                                  epsilon_r, epsilon_p, power_mw, volume_ml,
                                  thermal_rate, path_length_cm=1,
                                  initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
                                  thermal_forward_rate=0, photon_flux_mol_s=None,
                                  irradiation_wavelength_nm=None):
    """Fit quantum yields directly to measured absorbance spectra."""
    initial_r = np.trapezoid(absorbance[:, 0], wavelengths_nm) / np.trapezoid(
        epsilon_r * path_length_cm, wavelengths_nm
    )
    initial = np.array([initial_r, 0.0])
    target = absorbance.T
    scale = max(float(np.max(np.abs(target))), np.finfo(float).eps)

    def residual(model):
        model_absorbance = model @ np.vstack((epsilon_r, epsilon_p)) * path_length_cm
        return ((model_absorbance - target) / scale).ravel()

    return _fit(wavelengths_nm, emission, initial, timestamps, epsilon_r, epsilon_p,
                power_mw, volume_ml, thermal_rate, path_length_cm,
                initial_yields, yield_bounds, residual, thermal_forward_rate,
                photon_flux_mol_s, irradiation_wavelength_nm)


def fit_quantum_yields_ode_absorbance(
        wavelengths_nm, emission, absorbance, timestamps, epsilon_r, epsilon_p,
        power_mw, volume_ml, thermal_rate, path_length_cm=1,
        initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
        initial_concentrations=None, baseline_order=1, robust_loss_scale=0.02,
        thermal_forward_rate=0, photon_flux_mol_s=None,
        irradiation_wavelength_nm=None):
    """Jointly fit full-spectrum absorbance, yields, and initial composition."""
    if baseline_order not in {-1, 0, 1}:
        raise ValueError("baseline_order must be -1, 0, or 1")
    if robust_loss_scale <= 0:
        raise ValueError("robust_loss_scale must be positive")
    target = np.asarray(absorbance, float).T
    if initial_concentrations is None:
        initial_r = np.trapezoid(target[0], wavelengths_nm) / np.trapezoid(
            epsilon_r * path_length_cm, wavelengths_nm
        )
        initial_concentrations = np.array([initial_r, 0.0])
    initial_concentrations = np.asarray(initial_concentrations, float)
    concentration_scale = max(float(initial_concentrations.sum()), np.finfo(float).eps)
    initial_fraction = np.clip(initial_concentrations[0] / concentration_scale, 1e-6, 1 - 1e-6)
    irradiation = _irradiation_inputs(
        wavelengths_nm, emission, epsilon_r, epsilon_p, power_mw,
        photon_flux_mol_s, irradiation_wavelength_nm,
    )
    wavelengths_m, irradiation_epsilon_r, irradiation_epsilon_p, photons = irradiation
    absorbance_scale = max(float(np.max(np.abs(target))), np.finfo(float).eps)

    def evaluate(values, return_model=False):
        yields = values[:2]
        total = concentration_scale * np.exp(values[2])
        initial = total * np.array([values[3], 1 - values[3]])
        concentrations = _solve(
            yields, initial, timestamps, wavelengths_m,
            irradiation_epsilon_r, irradiation_epsilon_p,
            path_length_cm, photons, volume_ml, thermal_rate, thermal_forward_rate,
        )
        model = concentrations @ np.vstack((epsilon_r, epsilon_p)) * path_length_cm
        correction = _baseline_correction(target - model, wavelengths_nm, baseline_order)
        if return_model:
            return yields, concentrations, correction
        return ((model + correction - target) / absorbance_scale).ravel()

    initial = np.array([*initial_yields, 0.0, initial_fraction])
    lower = np.array([yield_bounds[0], yield_bounds[0], np.log(0.1), 0.0])
    upper = np.array([yield_bounds[1], yield_bounds[1], np.log(10.0), 1.0])
    fit = least_squares(
        evaluate, initial, bounds=(lower, upper), loss="soft_l1",
        f_scale=robust_loss_scale,
    )
    dof = max(fit.fun.size - fit.x.size, 1)
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * np.dot(fit.fun, fit.fun) / dof
    yields, concentrations, correction = evaluate(fit.x, return_model=True)
    return YieldFit(
        yields, np.sqrt(np.maximum(np.diag(covariance)[:2], 0)),
        concentrations, correction, bool(fit.success), str(fit.message),
        _jacobian_condition(fit.jac), _active_yield_bounds(yields, yield_bounds),
    )


def fit_quantum_yields_nipe(
        wavelengths_nm, emission, absorbance, concentrations, timestamps,
        epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
        path_length_cm=1, initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
        thermal_forward_rate=0, photon_flux_mol_s=None,
        irradiation_wavelength_nm=None, minimum_interval_fraction=0.05):
    """Fit apparent A <=> B yields with normalized integrated photon balances.

    Every retained ``t1``/``t2`` interval is divided by its tracked absorbed
    photon dose.  This is the ratio normalization at the heart of NIPE: no
    physical ``t=0`` or infinite-conversion boundary is required, and the
    measured time-dependent optical density supplies the inner-filter
    correction.  The reaction coordinate ``(B - A) / 2`` removes common-mode
    A+B drift so that an *apparent* isomerization yield can still be estimated
    when the separate NIPE balance diagnostic flags degradation or another
    side process.
    """
    concentrations = np.asarray(concentrations, float)
    absorbance = np.asarray(absorbance, float)
    timestamps = np.asarray(timestamps, float)
    if concentrations.shape != (len(timestamps), 2):
        raise ValueError("NIPE requires one A/B concentration pair per timestamp")
    if absorbance.shape != (len(wavelengths_nm), len(timestamps)):
        raise ValueError("NIPE absorbance dimensions do not match wavelengths and timestamps")
    if len(timestamps) < 3 or np.any(np.diff(timestamps) <= 0):
        raise ValueError("NIPE requires at least three strictly increasing timestamps")
    if not 0 <= minimum_interval_fraction < 1:
        raise ValueError("minimum_interval_fraction must be at least 0 and less than 1")

    irradiation = _irradiation_inputs(
        wavelengths_nm, emission, epsilon_r, epsilon_p, power_mw,
        photon_flux_mol_s, irradiation_wavelength_nm,
    )
    wavelengths_m, irradiation_epsilon_r, irradiation_epsilon_p, photons = irradiation
    if len(wavelengths_m) == 1:
        optical_density = np.array([
            np.interp(irradiation_wavelength_nm, wavelengths_nm, spectrum)
            for spectrum in absorbance.T
        ])[:, None]
    else:
        optical_density = absorbance.T
    optical_density = np.maximum(optical_density, 0)

    contribution_r = (path_length_cm * concentrations[:, 0, None]
                      * irradiation_epsilon_r[None, :])
    contribution_p = (path_length_cm * concentrations[:, 1, None]
                      * irradiation_epsilon_p[None, :])
    # Work wavelength by wavelength when the LED is broadband.
    fractions = np.stack((contribution_r, contribution_p), axis=-1)
    fractions = np.divide(
        fractions, optical_density[..., None], out=np.zeros_like(fractions),
        where=optical_density[..., None] > np.finfo(float).eps,
    )
    fractions = np.clip(fractions, 0, 1)
    assigned = fractions.sum(axis=-1)
    over_assigned = assigned > 1
    fractions[over_assigned] /= assigned[over_assigned, None]
    absorbed = ((1 - 10 ** -optical_density) * photons[None, :] / volume_ml)
    absorbed_r = fractions[..., 0] * absorbed
    absorbed_p = fractions[..., 1] * absorbed
    if len(wavelengths_m) == 1:
        rates = np.column_stack((absorbed_r[:, 0], absorbed_p[:, 0]))
    else:
        rates = np.column_stack((
            np.trapezoid(absorbed_r, wavelengths_m, axis=1),
            np.trapezoid(absorbed_p, wavelengths_m, axis=1),
        ))
    cumulative_photons = _cumulative_trapezoid(rates, timestamps)

    thermal_xi_rate = (thermal_forward_rate * concentrations[:, 0]
                       - thermal_rate * concentrations[:, 1])
    cumulative_thermal_xi = _cumulative_trapezoid(
        thermal_xi_rate[:, None], timestamps
    )[:, 0]
    xi = (concentrations[:, 1] - concentrations[:, 0]) / 2
    total_exposure = float(cumulative_photons[-1].sum())
    if not np.isfinite(total_exposure) or total_exposure <= 0:
        raise ValueError("NIPE found no positive absorbed photon exposure")
    minimum_exposure = total_exposure * minimum_interval_fraction

    design, observed = [], []
    for first in range(len(timestamps) - 1):
        for last in range(first + 1, len(timestamps)):
            exposure = cumulative_photons[last] - cumulative_photons[first]
            normalization = float(exposure.sum())
            if normalization <= max(minimum_exposure, np.finfo(float).eps):
                continue
            delta_xi = (xi[last] - xi[first]
                        - cumulative_thermal_xi[last] + cumulative_thermal_xi[first])
            design.append([exposure[0] / normalization,
                           -exposure[1] / normalization])
            observed.append(delta_xi / normalization)
    design = np.asarray(design, float)
    observed = np.asarray(observed, float)
    if len(observed) < 3 or np.linalg.matrix_rank(design) < 2:
        raise ValueError("NIPE intervals do not independently identify both quantum yields")

    initial = np.clip(np.asarray(initial_yields, float), yield_bounds[0], yield_bounds[1])
    scale = max(float(np.median(np.abs(observed))), 0.01)
    fit = least_squares(
        lambda values: design @ values - observed,
        initial, bounds=yield_bounds, loss="soft_l1", f_scale=scale,
    )
    residual = design @ fit.x - observed
    dof = max(len(observed) - len(fit.x), 1)
    covariance = (np.linalg.pinv(fit.jac.T @ fit.jac)
                  * np.dot(residual, residual) / dof)

    fitted_xi = (xi[0] + cumulative_thermal_xi
                 + cumulative_photons[:, 0] * fit.x[0]
                 - cumulative_photons[:, 1] * fit.x[1])
    totals = concentrations.sum(axis=1)
    fitted_concentrations = np.column_stack((
        totals / 2 - fitted_xi,
        totals / 2 + fitted_xi,
    ))
    metadata = NIPEMetadata(
        interval_count=len(observed),
        minimum_interval_fraction=minimum_interval_fraction,
        normalized_residual_rmse=float(np.sqrt(np.mean(residual ** 2))),
    )
    return YieldFit(
        fit.x, np.sqrt(np.maximum(np.diag(covariance), 0)),
        fitted_concentrations, optimizer_success=bool(fit.success),
        optimizer_message=str(fit.message),
        jacobian_condition=_jacobian_condition(fit.jac),
        active_bounds=_active_yield_bounds(fit.x, yield_bounds), nipe=metadata,
    )


def analyze_nipe_time_windows(
        wavelengths_nm, emission, absorbance, concentrations, timestamps,
        epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
        path_length_cm=1, initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
        thermal_forward_rate=0, photon_flux_mol_s=None,
        irradiation_wavelength_nm=None, full_trace_values=None,
        plateau_rate_fraction=0.05, plateau_run_length=3,
        maximum_window_condition=1e4):
    """Estimate zero-exposure yields from identifiable pre-plateau NIPE windows.

    Plateau onset is the first sustained run whose absolute composition change
    per acquisition is at most ``plateau_rate_fraction`` of the initial change.
    Windows remain fully before that onset, contain four to six spectra, span a
    measurable composition change, and are rejected when their two-yield
    Jacobian is poorly conditioned or a yield reaches a configured bound.
    """
    timestamps = np.asarray(timestamps, float)
    concentrations = np.asarray(concentrations, float)
    absorbance = np.asarray(absorbance, float)
    if len(timestamps) < 6:
        return None
    if not 0 < plateau_rate_fraction < 1:
        raise ValueError("plateau_rate_fraction must be greater than 0 and less than 1")
    if plateau_run_length < 2:
        raise ValueError("plateau_run_length must be at least 2")

    totals = concentrations.sum(axis=1)
    product_fraction = np.divide(
        concentrations[:, 1], totals, out=np.zeros_like(totals), where=totals > 0,
    )
    plateau_index = _nipe_plateau_index(
        product_fraction, plateau_rate_fraction, plateau_run_length
    )
    plateau_detected = plateau_index < len(timestamps)
    pre_plateau_count = plateau_index if plateau_detected else len(timestamps)
    if pre_plateau_count < 6:
        return None

    window_point_count = min(6, max(4, int(round(pre_plateau_count * 0.4))))
    starts = np.arange(pre_plateau_count - window_point_count + 1, dtype=int)
    if len(starts) > 8:
        starts = np.unique(np.linspace(0, starts[-1], 8).round().astype(int))
    pre_plateau_span = float(np.ptp(product_fraction[:pre_plateau_count]))
    minimum_window_span = max(0.01, 0.05 * pre_plateau_span)

    accepted = []
    for first in starts:
        last = first + window_point_count
        if np.ptp(product_fraction[first:last]) < minimum_window_span:
            continue
        try:
            fit = fit_quantum_yields_nipe(
                wavelengths_nm, emission, absorbance[:, first:last],
                concentrations[first:last], timestamps[first:last], epsilon_r,
                epsilon_p, power_mw, volume_ml, thermal_rate, path_length_cm,
                initial_yields, yield_bounds, thermal_forward_rate,
                photon_flux_mol_s, irradiation_wavelength_nm,
            )
        except ValueError:
            continue
        if (not fit.optimizer_success or any(fit.active_bounds)
                or not np.isfinite(fit.jacobian_condition)
                or fit.jacobian_condition > maximum_window_condition):
            continue
        accepted.append((first, last, fit))
    if len(accepted) < 3:
        return None

    starts = np.array([item[0] for item in accepted], int)
    stops = np.array([item[1] for item in accepted], int)
    values = np.array([item[2].values for item in accepted])
    errors = np.array([item[2].standard_errors for item in accepted])
    conditions = np.array([item[2].jacobian_condition for item in accepted])
    relative_time = timestamps - timestamps[0]
    midpoints = (relative_time[starts] + relative_time[stops - 1]) / 2
    design = np.column_stack((midpoints, np.ones_like(midpoints)))
    extrapolated, extrapolated_errors = [], []
    for index in range(2):
        coefficients = np.linalg.lstsq(design, values[:, index], rcond=None)[0]
        residual = values[:, index] - design @ coefficients
        dof = max(len(values) - 2, 1)
        covariance = (np.linalg.pinv(design.T @ design)
                      * np.dot(residual, residual) / dof)
        intercept = float(np.clip(coefficients[1], *yield_bounds))
        regression_error = float(np.sqrt(max(covariance[1, 1], 0)))
        local_error = float(np.sqrt(np.mean(errors[:, index] ** 2)))
        extrapolated.append(intercept)
        extrapolated_errors.append(np.hypot(regression_error, local_error))
    extrapolated = np.asarray(extrapolated)
    extrapolated_errors = np.asarray(extrapolated_errors)
    full_trace_values = np.asarray(
        full_trace_values if full_trace_values is not None else values[-1], float
    )
    change = np.divide(
        full_trace_values - extrapolated, extrapolated,
        out=np.full(2, np.nan), where=np.abs(extrapolated) > np.finfo(float).eps,
    )
    plateau_time = (float(relative_time[plateau_index]) if plateau_detected
                    else float(relative_time[-1]))
    return NIPEWindowAnalysis(
        plateau_detected=plateau_detected,
        plateau_time_s=plateau_time,
        analysis_end_time_s=float(relative_time[pre_plateau_count - 1]),
        window_point_count=window_point_count,
        window_duration_s=float(np.median(
            relative_time[stops - 1] - relative_time[starts]
        )),
        window_start_s=relative_time[starts],
        window_end_s=relative_time[stops - 1],
        window_midpoint_s=midpoints,
        window_values=values,
        window_standard_errors=errors,
        window_jacobian_conditions=conditions,
        extrapolated_values=extrapolated,
        extrapolated_standard_errors=extrapolated_errors,
        full_trace_values=full_trace_values,
        full_trace_change_fraction=change,
    )


def _nipe_plateau_index(product_fraction, rate_fraction=0.05, run_length=3):
    """Return the first plateau point, or the series length when none is found."""
    product_fraction = np.asarray(product_fraction, float)
    increments = np.abs(np.diff(product_fraction))
    if len(increments) < run_length + 2:
        return len(product_fraction)
    initial_change = float(np.max(increments[:min(2, len(increments))]))
    if not np.isfinite(initial_change) or initial_change <= np.finfo(float).eps:
        return len(product_fraction)
    tail_count = max(run_length, len(increments) // 5)
    noise = float(np.median(increments[-tail_count:]))
    threshold = max(
        rate_fraction * initial_change,
        min(3 * noise, 0.25 * initial_change),
    )
    below = increments <= threshold
    for start in range(2, len(below) - run_length + 1):
        if np.all(below[start:start + run_length]):
            return start + 1
    return len(product_fraction)


def _fit(wavelengths_nm, emission, initial, timestamps, epsilon_r, epsilon_p,
         power_mw, volume_ml, thermal_rate, path_length_cm, initial_yields,
         yield_bounds, residual_function, thermal_forward_rate=0,
         photon_flux_mol_s=None, irradiation_wavelength_nm=None):
    irradiation = _irradiation_inputs(
        wavelengths_nm, emission, epsilon_r, epsilon_p, power_mw,
        photon_flux_mol_s, irradiation_wavelength_nm,
    )
    wavelengths_m, irradiation_epsilon_r, irradiation_epsilon_p, photons = irradiation

    def residual(values):
        model = _solve(
            values, initial, timestamps, wavelengths_m,
            irradiation_epsilon_r, irradiation_epsilon_p, path_length_cm,
            photons, volume_ml, thermal_rate, thermal_forward_rate,
        )
        return residual_function(model)

    fit = least_squares(residual, initial_yields, bounds=yield_bounds)
    dof = max(fit.fun.size - fit.x.size, 1)
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * np.dot(fit.fun, fit.fun) / dof
    fitted = _solve(
        fit.x, initial, timestamps, wavelengths_m,
        irradiation_epsilon_r, irradiation_epsilon_p, path_length_cm,
        photons, volume_ml, thermal_rate, thermal_forward_rate,
    )
    return YieldFit(
        fit.x, np.sqrt(np.maximum(np.diag(covariance), 0)), fitted,
        optimizer_success=bool(fit.success), optimizer_message=str(fit.message),
        jacobian_condition=_jacobian_condition(fit.jac),
        active_bounds=_active_yield_bounds(fit.x, yield_bounds),
    )


def _cumulative_trapezoid(values, coordinates):
    values = np.asarray(values, float)
    coordinates = np.asarray(coordinates, float)
    increments = ((values[:-1] + values[1:]) / 2
                  * np.diff(coordinates)[:, None])
    return np.vstack((np.zeros((1, values.shape[1])), np.cumsum(increments, axis=0)))


def _jacobian_condition(jacobian):
    condition = float(np.linalg.cond(jacobian))
    return condition if np.isfinite(condition) else np.inf


def _active_yield_bounds(values, bounds):
    tolerance = max((bounds[1] - bounds[0]) * 1e-5, 1e-10)
    return tuple(bool(value <= bounds[0] + tolerance or value >= bounds[1] - tolerance)
                 for value in values[:2])


def _photon_flux(wavelengths_nm, emission, power_mw):
    wavelengths_m = wavelengths_nm * 1e-9
    area = np.trapezoid(emission, wavelengths_m)
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Processed LED emission must have a positive finite area")
    normalized_emission = emission / area
    photons = (normalized_emission * power_mw * 1e-3
               / (H * C / wavelengths_m) / AVOGADRO * 1000)
    return wavelengths_m, photons


def power_mw_to_photon_flux_mol_s(power_mw, wavelength_nm):
    """Convert monochromatic optical power to mol photons per second."""
    power_mw = float(power_mw)
    wavelength_m = float(wavelength_nm) * 1e-9
    if not np.isfinite(power_mw) or power_mw <= 0:
        raise ValueError("Optical power must be positive and finite")
    if not np.isfinite(wavelength_m) or wavelength_m <= 0:
        raise ValueError("Irradiation wavelength must be positive and finite")
    return power_mw * 1e-3 / (H * C / wavelength_m) / AVOGADRO


def photon_flux_mol_s_to_power_mw(photon_flux_mol_s, wavelength_nm):
    """Convert a monochromatic molar photon flux to equivalent optical power."""
    photon_flux_mol_s = float(photon_flux_mol_s)
    wavelength_m = float(wavelength_nm) * 1e-9
    if not np.isfinite(photon_flux_mol_s) or photon_flux_mol_s <= 0:
        raise ValueError("Photon flux must be positive and finite")
    if not np.isfinite(wavelength_m) or wavelength_m <= 0:
        raise ValueError("Irradiation wavelength must be positive and finite")
    return photon_flux_mol_s * AVOGADRO * H * C / wavelength_m * 1e3


def _irradiation_inputs(wavelengths_nm, emission, epsilon_r, epsilon_p,
                        power_mw, photon_flux_mol_s, irradiation_wavelength_nm):
    wavelengths_nm = np.asarray(wavelengths_nm, float)
    epsilon_r = np.asarray(epsilon_r, float)
    epsilon_p = np.asarray(epsilon_p, float)
    if photon_flux_mol_s is None:
        wavelengths_m, photons = _photon_flux(wavelengths_nm, emission, power_mw)
        return wavelengths_m, epsilon_r, epsilon_p, photons

    photon_flux_mol_s = float(photon_flux_mol_s)
    irradiation_wavelength_nm = float(irradiation_wavelength_nm)
    if not np.isfinite(photon_flux_mol_s) or photon_flux_mol_s <= 0:
        raise ValueError("Photon flux must be positive and finite")
    if (not np.isfinite(irradiation_wavelength_nm) or irradiation_wavelength_nm <= 0
            or irradiation_wavelength_nm < wavelengths_nm[0]
            or irradiation_wavelength_nm > wavelengths_nm[-1]):
        raise ValueError("Irradiation wavelength must lie within the fitted wavelength range")
    irradiation_epsilon_r = np.array([
        np.interp(irradiation_wavelength_nm, wavelengths_nm, epsilon_r)
    ])
    irradiation_epsilon_p = np.array([
        np.interp(irradiation_wavelength_nm, wavelengths_nm, epsilon_p)
    ])
    return (
        np.array([irradiation_wavelength_nm * 1e-9]),
        irradiation_epsilon_r,
        irradiation_epsilon_p,
        np.array([photon_flux_mol_s * 1000]),
    )


def extrapolate_photostationary_state(
        wavelengths_nm, emission, total_concentration, yields, epsilon_r,
        epsilon_p, power_mw, volume_ml, thermal_rate, path_length_cm=1,
        thermal_forward_rate=0, photon_flux_mol_s=None,
        irradiation_wavelength_nm=None):
    """Find the steady composition under continued constant irradiation."""
    total_concentration = float(total_concentration)
    if not np.isfinite(total_concentration) or total_concentration <= 0:
        raise ValueError("Total concentration must be positive and finite")
    irradiation = _irradiation_inputs(
        wavelengths_nm, emission, epsilon_r, epsilon_p, power_mw,
        photon_flux_mol_s, irradiation_wavelength_nm,
    )
    wavelengths_m, irradiation_epsilon_r, irradiation_epsilon_p, photons = irradiation

    def reactant_rate(reactant_fraction):
        concentrations = total_concentration * np.array(
            [reactant_fraction, 1 - reactant_fraction]
        )
        return _rates(
            concentrations, 0, wavelengths_m, *yields,
            irradiation_epsilon_r, irradiation_epsilon_p,
            path_length_cm, photons, volume_ml, thermal_rate, thermal_forward_rate,
        )[0]

    rate_at_product = reactant_rate(0.0)
    rate_at_reactant = reactant_rate(1.0)
    if rate_at_product == 0:
        reactant_fraction = 0.0
    elif rate_at_reactant == 0:
        reactant_fraction = 1.0
    elif rate_at_product * rate_at_reactant < 0:
        reactant_fraction = brentq(reactant_rate, 0.0, 1.0)
    else:
        raise ValueError("The fitted model has no physical steady composition")
    return total_concentration * np.array(
        [reactant_fraction, 1 - reactant_fraction]
    )


def _baseline_correction(delta, wavelengths_nm, order):
    if order < 0:
        return np.zeros_like(delta)
    coordinate = np.asarray(wavelengths_nm, float)
    span = np.ptp(coordinate)
    coordinate = ((coordinate - coordinate.mean()) / span
                  if span > 0 else np.zeros_like(coordinate))
    design = np.column_stack([coordinate ** degree for degree in range(order + 1)])
    coefficients = np.linalg.lstsq(design, delta.T, rcond=None)[0]
    return (design @ coefficients).T


def _solve(yields, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
           path_length_cm, photons, volume_ml, thermal_rate, thermal_forward_rate=0):
    args = (wavelengths_m, *yields, epsilon_r, epsilon_p, path_length_cm,
            photons, volume_ml, thermal_rate, thermal_forward_rate)
    return odeint(_rates, initial, timestamps, args=args)


def _rates(concentrations, _, wavelengths_m, yield_rp, yield_pr,
           epsilon_r, epsilon_p, path_length_cm, photons, volume_ml, thermal_rate,
           thermal_forward_rate=0):
    reactant, product = concentrations
    absorbance_r = reactant * epsilon_r
    absorbance_p = product * epsilon_p
    total_species = absorbance_r + absorbance_p
    total = path_length_cm * total_species
    fraction_r = np.divide(absorbance_r, total_species, out=np.zeros_like(total),
                           where=total_species != 0)
    fraction_p = np.divide(absorbance_p, total_species, out=np.zeros_like(total),
                           where=total_species != 0)
    absorbed = (1 - 10 ** -total) * photons / volume_ml
    if len(wavelengths_m) == 1:
        rate_r = -float(fraction_r[0] * absorbed[0])
        rate_p = float(fraction_p[0] * absorbed[0])
    else:
        rate_r = -np.trapezoid(fraction_r * absorbed, wavelengths_m)
        rate_p = np.trapezoid(fraction_p * absorbed, wavelengths_m)
    change = (yield_rp * rate_r + yield_pr * rate_p + thermal_rate * product
              - thermal_forward_rate * reactant)
    return change, -change
