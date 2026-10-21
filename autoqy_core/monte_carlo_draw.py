"""One Monte Carlo input draw and its deterministic kinetic fit."""

from dataclasses import replace
import warnings

import numpy as np

from .epsilon_uncertainty import _absorbance_residual_rmse


def prepare_trial(data, reactant, product, model, residuals, rows, rng, sample_pair):
    from .epsilon_monte_carlo import _positive_normal

    er, ep = sample_pair(rng)
    absorbance = data.absorbance.copy()
    bootstrap = residuals[rng.integers(0, len(residuals), size=len(residuals))]
    absorbance[rows] = (model + bootstrap).T
    power, flux = data.power_mw, data.photon_flux_mol_s
    if flux is None:
        power = float(_positive_normal(power, data.power_error_mw, rng.random()))
    else:
        flux = float(_positive_normal(flux, data.photon_flux_error_mol_s, rng.random()))
    return replace(
        data, epsilon_r=(reactant.wavelengths, er), epsilon_p=(product.wavelengths, ep),
        absorbance=absorbance, power_mw=power, photon_flux_mol_s=flux,
        power_error_mw=0, photon_flux_error_mol_s=0,
    )


def fit_trial(trial_data):
    from .pipeline import run_analysis_pipeline

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Initial product is .*", category=RuntimeWarning)
        trial = run_analysis_pipeline(trial_data, include_power_variations=False)
    fit = trial.yield_fit
    if not fit.optimizer_success or not np.isfinite(fit.values).all():
        raise ValueError("Kinetic optimizer failed or returned non-finite yields")
    total = fit.concentrations.sum(axis=1)
    fraction = np.divide(fit.concentrations[:, 0], total,
                         out=np.zeros_like(total), where=total != 0)
    residual = trial.concentration_fit.fractions[:, 0] - fraction
    return (fit.values, fit.standard_errors, fit.active_bounds,
            fit.jacobian_condition, trial.concentration_fit.concentrations,
            fit.concentrations, residual,
            _absorbance_residual_rmse(trial, trial_data.path_length_cm))
