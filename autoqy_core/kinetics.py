"""Quantum-yield fitting independent of the user interface."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq, least_squares

H = 6.626e-34
C = 299792458
AVOGADRO = 6.022e23


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


def fit_quantum_yields(wavelengths_nm, emission, concentrations, timestamps,
                       epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
                       path_length_cm=1, initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
                       thermal_forward_rate=0):
    """Fit quantum yields to experimentally recovered concentration traces."""
    scale = max(float(np.max(np.abs(concentrations))), np.finfo(float).eps)

    def residual(model):
        return ((model - concentrations) / scale).ravel()

    return _fit(wavelengths_nm, emission, concentrations[0], timestamps,
                epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
                path_length_cm, initial_yields, yield_bounds, residual,
                thermal_forward_rate)


def fit_quantum_yields_absorbance(wavelengths_nm, emission, absorbance, timestamps,
                                  epsilon_r, epsilon_p, power_mw, volume_ml,
                                  thermal_rate, path_length_cm=1,
                                  initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
                                  thermal_forward_rate=0):
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
                initial_yields, yield_bounds, residual, thermal_forward_rate)


def fit_quantum_yields_ode_absorbance(
        wavelengths_nm, emission, absorbance, timestamps, epsilon_r, epsilon_p,
        power_mw, volume_ml, thermal_rate, path_length_cm=1,
        initial_yields=(0.5, 0.5), yield_bounds=(0, 1),
        initial_concentrations=None, baseline_order=1, robust_loss_scale=0.02,
        thermal_forward_rate=0):
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
    wavelengths_m, photons = _photon_flux(wavelengths_nm, emission, power_mw)
    absorbance_scale = max(float(np.max(np.abs(target))), np.finfo(float).eps)

    def evaluate(values, return_model=False):
        yields = values[:2]
        total = concentration_scale * np.exp(values[2])
        initial = total * np.array([values[3], 1 - values[3]])
        concentrations = _solve(
            yields, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
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


def _fit(wavelengths_nm, emission, initial, timestamps, epsilon_r, epsilon_p,
         power_mw, volume_ml, thermal_rate, path_length_cm, initial_yields,
         yield_bounds, residual_function, thermal_forward_rate=0):
    wavelengths_m, photons = _photon_flux(wavelengths_nm, emission, power_mw)

    def residual(values):
        model = _solve(values, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
                       path_length_cm, photons, volume_ml, thermal_rate,
                       thermal_forward_rate)
        return residual_function(model)

    fit = least_squares(residual, initial_yields, bounds=yield_bounds)
    dof = max(fit.fun.size - fit.x.size, 1)
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * np.dot(fit.fun, fit.fun) / dof
    fitted = _solve(fit.x, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
                    path_length_cm, photons, volume_ml, thermal_rate,
                    thermal_forward_rate)
    return YieldFit(
        fit.x, np.sqrt(np.maximum(np.diag(covariance), 0)), fitted,
        optimizer_success=bool(fit.success), optimizer_message=str(fit.message),
        jacobian_condition=_jacobian_condition(fit.jac),
        active_bounds=_active_yield_bounds(fit.x, yield_bounds),
    )


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


def extrapolate_photostationary_state(
        wavelengths_nm, emission, total_concentration, yields, epsilon_r,
        epsilon_p, power_mw, volume_ml, thermal_rate, path_length_cm=1,
        thermal_forward_rate=0):
    """Find the steady composition under continued constant irradiation."""
    total_concentration = float(total_concentration)
    if not np.isfinite(total_concentration) or total_concentration <= 0:
        raise ValueError("Total concentration must be positive and finite")
    wavelengths_m, photons = _photon_flux(wavelengths_nm, emission, power_mw)

    def reactant_rate(reactant_fraction):
        concentrations = total_concentration * np.array(
            [reactant_fraction, 1 - reactant_fraction]
        )
        return _rates(
            concentrations, 0, wavelengths_m, *yields, epsilon_r, epsilon_p,
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
    rate_r = -np.trapezoid(fraction_r * absorbed, wavelengths_m)
    rate_p = np.trapezoid(fraction_p * absorbed, wavelengths_m)
    change = (yield_rp * rate_r + yield_pr * rate_p + thermal_rate * product
              - thermal_forward_rate * reactant)
    return change, -change
