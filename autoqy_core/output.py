"""Write AutoQY results without GUI dependencies."""

import math
from pathlib import Path

import numpy as np
import pandas as pd


def format_value_uncertainty(value, uncertainty, two_digit_threshold=3):
    """Round uncertainty and value to a scientifically matching decimal place."""
    value, uncertainty = float(value), abs(float(uncertainty))
    if not math.isfinite(uncertainty) or uncertainty == 0:
        return f"{value:g}", f"{uncertainty:g}"
    exponent = math.floor(math.log10(uncertainty))
    leading = uncertainty / 10 ** exponent
    significant_digits = 2 if leading < two_digit_threshold else 1
    place = exponent - significant_digits + 1
    rounded_value = round(value, -place)
    rounded_uncertainty = round(uncertainty, -place)
    decimals = max(0, -place)
    return f"{rounded_value:.{decimals}f}", f"{rounded_uncertainty:.{decimals}f}"


def result_summary(result, data, irradiation_wavelength_nm):
    values = result.yield_fit.values * 100
    errors = result.yield_errors * 100
    last_composition = result.yield_fit.concentrations[-1]
    last_composition = last_composition / last_composition.sum() * 100
    extrapolated_pss = result.extrapolated_pss
    extrapolated_pss = extrapolated_pss / extrapolated_pss.sum() * 100
    formatted = [format_value_uncertainty(value, error, two_digit_threshold=2)
                 for value, error in zip(values, errors)]
    summary = {
        "schema_version": 2,
        "fit_method": result.fit_method,
        "quantum_yield_percent": {"R_to_P": float(values[0]), "P_to_R": float(values[1])},
        "quantum_yield_error_percent": {"R_to_P": float(errors[0]), "P_to_R": float(errors[1])},
        "quantum_yield_formatted_percent": {
            "R_to_P": {"value": formatted[0][0], "error": formatted[0][1]},
            "P_to_R": {"value": formatted[1][0], "error": formatted[1][1]},
        },
        "composition_at_last_timestamp_percent": {
            "time_s": float(data.timestamps[-1]),
            "reactant": float(last_composition[0]),
            "product": float(last_composition[1]),
        },
        "extrapolated_pss_percent": {
            "reactant": float(extrapolated_pss[0]),
            "product": float(extrapolated_pss[1]),
        },
        "experiment": {
            "volume_ml": data.volume_ml,
            "power_mw": data.power_mw,
            "power_error_mw": data.power_error_mw,
            "thermal_back_reaction_s_1": data.thermal_rate,
            "irradiation_wavelength_nm": irradiation_wavelength_nm,
            "path_length_cm": data.path_length_cm,
        },
    }
    uncertainty = result.epsilon_uncertainty
    if uncertainty is not None:
        optimizer_power = uncertainty.optimizer_power_errors * 100
        epsilon = uncertainty.epsilon_errors * 100
        combined = uncertainty.combined_errors * 100
        summary["quantum_yield_error_components_percent"] = {
            "optimizer_and_power": _yield_pair(optimizer_power),
            "epsilon": _yield_pair(epsilon),
            "combined": _yield_pair(combined),
        }
        summary["epsilon_uncertainty"] = {
            "method": uncertainty.method,
            "error_metric": uncertainty.error_metric,
            "bound_combination_count": uncertainty.bound_combination_count,
            "reactant_source_schema": uncertainty.reactant_source_schema,
            "product_source_schema": uncertainty.product_source_schema,
            "reactant_source_path": uncertainty.reactant_source_path,
            "product_source_path": uncertainty.product_source_path,
            "reactant_error_metric": uncertainty.reactant_error_metric,
            "product_error_metric": uncertainty.product_error_metric,
            "constrained_negative_points": {
                "reactant": uncertainty.constrained_negative_points[0],
                "product": uncertainty.constrained_negative_points[1],
            },
            "quantum_yield_minimum_percent": _yield_pair(
                uncertainty.epsilon_yield_minimum * 100
            ),
            "quantum_yield_maximum_percent": _yield_pair(
                uncertainty.epsilon_yield_maximum * 100
            ),
            "absorbance_residual_rmse_range": {
                "minimum": uncertainty.absorbance_residual_rmse_minimum,
                "maximum": uncertainty.absorbance_residual_rmse_maximum,
            },
        }
    return summary


def write_results(path, result, data, irradiation_wavelength_nm):
    summary = result_summary(result, data, irradiation_wavelength_nm)
    formatted = summary["quantum_yield_formatted_percent"]
    last_composition = summary["composition_at_last_timestamp_percent"]
    extrapolated_pss = summary["extrapolated_pss_percent"]
    low, high = data.wavelength_limits
    method = {
        "concentrations": "Concentrations",
        "emission": "Emission (legacy)",
        "regularized_concentrations": "Regularized concentrations",
        "ode_absorbance": "Full-spectrum ODE absorbance",
    }[result.fit_method]
    epsilon_text = "\n"
    if result.epsilon_uncertainty is not None:
        components = summary["quantum_yield_error_components_percent"]
        metadata = summary["epsilon_uncertainty"]
        epsilon_text = f"""Error component optimizer + power R_to_P (%): {components['optimizer_and_power']['R_to_P']:g}
Error component optimizer + power P_to_R (%): {components['optimizer_and_power']['P_to_R']:g}
Error component epsilon R_to_P (%): {components['epsilon']['R_to_P']:g}
Error component epsilon P_to_R (%): {components['epsilon']['P_to_R']:g}
Epsilon uncertainty method: {metadata['method']}
Reactant epsilon error metric: {metadata['reactant_error_metric']}
Product epsilon error metric: {metadata['product_error_metric']}
Epsilon bound combinations: {metadata['bound_combination_count']}
Reactant epsilon values constrained to zero: {metadata['constrained_negative_points']['reactant']}
NMR epsilon values constrained to zero: {metadata['constrained_negative_points']['product']}

"""
    text = f"""Composition at the last timestamp (s): {last_composition['time_s']:g}
Composition at the last timestamp - Reactant (%): {last_composition['reactant']:.1f}
Composition at the last timestamp - Product (%): {last_composition['product']:.1f}
Extrapolated PSS - Reactant (%): {extrapolated_pss['reactant']:.1f}
Extrapolated PSS - Product (%): {extrapolated_pss['product']:.1f}
QY_AB_opt (%): {formatted['R_to_P']['value']}
QY_BA_opt (%): {formatted['P_to_R']['value']}
error_QY_AB (%): {formatted['R_to_P']['error']}
error_QY_BA (%): {formatted['P_to_R']['error']}
{epsilon_text}Volume (ml): {data.volume_ml:g}
k thermal back-reaction (s-1): {data.thermal_rate:g}
Power average (mW): {data.power_mw:g}
Power error (mW): {data.power_error_mw:g}
Wavelength of irradiation: {irradiation_wavelength_nm:g}

Calculation Method: Integration
ODE Solving Method: {method}
Baseline Correction LED Spectrum: {'ON' if data.baseline_correct_led else 'OFF'}
Wavelength Range: {low:g}-{high:g}
"""
    Path(path).write_text(text, encoding="utf-8")


def _yield_pair(values):
    return {"R_to_P": float(values[0]), "P_to_R": float(values[1])}


def write_detailed_data(stem, result, data):
    """Write time traces and long-form spectral residuals."""
    stem = Path(stem)
    measured = result.concentration_fit.concentrations
    fitted = result.yield_fit.concentrations
    measured_fractions = result.concentration_fit.fractions
    fitted_totals = fitted.sum(axis=1)
    fitted_fractions = np.divide(fitted, fitted_totals[:, None],
                                 out=np.zeros_like(fitted),
                                 where=fitted_totals[:, None] != 0)
    concentration_residual = measured - fitted
    fraction_residual = measured_fractions - fitted_fractions
    columns = {
        "time_s": data.timestamps,
        "reactant_concentration_data_M": measured[:, 0],
        "product_concentration_data_M": measured[:, 1],
        "reactant_concentration_fit_M": fitted[:, 0],
        "product_concentration_fit_M": fitted[:, 1],
        "reactant_concentration_residual_M": concentration_residual[:, 0],
        "product_concentration_residual_M": concentration_residual[:, 1],
        "reactant_fraction_data": measured_fractions[:, 0],
        "product_fraction_data": measured_fractions[:, 1],
        "reactant_fraction_fit": fitted_fractions[:, 0],
        "product_fraction_fit": fitted_fractions[:, 1],
        "reactant_fraction_residual": fraction_residual[:, 0],
        "product_fraction_residual": fraction_residual[:, 1],
    }
    uncertainty = result.epsilon_uncertainty
    if uncertainty is not None:
        for species, index in (("reactant", 0), ("product", 1)):
            columns[f"{species}_concentration_data_epsilon_min_M"] = (
                uncertainty.concentration_data_minimum[:, index]
            )
            columns[f"{species}_concentration_data_epsilon_max_M"] = (
                uncertainty.concentration_data_maximum[:, index]
            )
            columns[f"{species}_concentration_fit_epsilon_min_M"] = (
                uncertainty.concentration_fit_minimum[:, index]
            )
            columns[f"{species}_concentration_fit_epsilon_max_M"] = (
                uncertainty.concentration_fit_maximum[:, index]
            )
        columns["reactant_fraction_residual_epsilon_min"] = (
            uncertainty.fraction_residual_minimum
        )
        columns["reactant_fraction_residual_epsilon_max"] = (
            uncertainty.fraction_residual_maximum
        )
    traces = pd.DataFrame(columns)
    traces_path = stem.parent / f"{stem.name}_traces.tsv"
    traces.to_csv(traces_path, sep="\t", index=False)

    measured_absorbance = result.absorbance.T
    concentration_reconstruction = result.concentration_fit.reconstructed_absorbance
    kinetic_model = (fitted @ np.vstack((result.epsilon_r, result.epsilon_p))
                     * data.path_length_cm)
    correction = result.yield_fit.absorbance_correction
    if correction is None:
        correction = np.zeros_like(kinetic_model)
    kinetic_fit = kinetic_model + correction
    count_time, count_wavelength = measured_absorbance.shape
    spectra = pd.DataFrame({
        "time_s": np.repeat(data.timestamps, count_wavelength),
        "wavelength_nm": np.tile(result.wavelengths, count_time),
        "absorbance_measured": measured_absorbance.ravel(),
        "absorbance_concentration_reconstruction": concentration_reconstruction.ravel(),
        "absorbance_kinetic_model": kinetic_model.ravel(),
        "absorbance_baseline_correction": correction.ravel(),
        "absorbance_kinetic_fit": kinetic_fit.ravel(),
        "concentration_reconstruction_residual": (
            measured_absorbance - concentration_reconstruction).ravel(),
        "kinetic_fit_residual": (measured_absorbance - kinetic_fit).ravel(),
    })
    spectra_path = stem.parent / f"{stem.name}_spectra.tsv"
    spectra.to_csv(spectra_path, sep="\t", index=False)
    return traces_path, spectra_path
