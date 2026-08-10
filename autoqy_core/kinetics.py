"""Quantum-yield fitting independent of the user interface."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import least_squares

H = 6.626e-34
C = 299792458
AVOGADRO = 6.022e23


@dataclass(frozen=True)
class YieldFit:
    values: np.ndarray
    standard_errors: np.ndarray
    concentrations: np.ndarray


def fit_quantum_yields(wavelengths_nm, emission, concentrations, timestamps,
                       epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
                       path_length_cm=1, initial_yields=(0.5, 0.5), yield_bounds=(0, 1)):
    """Fit quantum yields to experimentally recovered concentration traces."""
    scale = max(float(np.max(np.abs(concentrations))), np.finfo(float).eps)

    def residual(model):
        return ((model - concentrations) / scale).ravel()

    return _fit(wavelengths_nm, emission, concentrations[0], timestamps,
                epsilon_r, epsilon_p, power_mw, volume_ml, thermal_rate,
                path_length_cm, initial_yields, yield_bounds, residual)


def fit_quantum_yields_absorbance(wavelengths_nm, emission, absorbance, timestamps,
                                  epsilon_r, epsilon_p, power_mw, volume_ml,
                                  thermal_rate, path_length_cm=1,
                                  initial_yields=(0.5, 0.5), yield_bounds=(0, 1)):
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
                initial_yields, yield_bounds, residual)


def _fit(wavelengths_nm, emission, initial, timestamps, epsilon_r, epsilon_p,
         power_mw, volume_ml, thermal_rate, path_length_cm, initial_yields,
         yield_bounds, residual_function):
    wavelengths_m = wavelengths_nm * 1e-9
    area = np.trapezoid(emission, wavelengths_m)
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Processed LED emission must have a positive finite area")
    normalized_emission = emission / area
    photons = (normalized_emission * power_mw * 1e-3
               / (H * C / wavelengths_m) / AVOGADRO * 1000)

    def residual(values):
        model = _solve(values, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
                       path_length_cm, photons, volume_ml, thermal_rate)
        return residual_function(model)

    fit = least_squares(residual, initial_yields, bounds=yield_bounds)
    dof = max(fit.fun.size - fit.x.size, 1)
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * np.dot(fit.fun, fit.fun) / dof
    fitted = _solve(fit.x, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
                    path_length_cm, photons, volume_ml, thermal_rate)
    return YieldFit(fit.x, np.sqrt(np.maximum(np.diag(covariance), 0)), fitted)


def _solve(yields, initial, timestamps, wavelengths_m, epsilon_r, epsilon_p,
           path_length_cm, photons, volume_ml, thermal_rate):
    args = (wavelengths_m, *yields, epsilon_r, epsilon_p, path_length_cm,
            photons, volume_ml, thermal_rate)
    return odeint(_rates, initial, timestamps, args=args)


def _rates(concentrations, _, wavelengths_m, yield_rp, yield_pr,
           epsilon_r, epsilon_p, path_length_cm, photons, volume_ml, thermal_rate):
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
    change = yield_rp * rate_r + yield_pr * rate_p + thermal_rate * product
    return change, -change
