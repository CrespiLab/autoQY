"""Write AutoQY results without GUI dependencies."""

import math
from pathlib import Path

import numpy as np
import pandas as pd


def format_value_uncertainty(value, uncertainty):
    """Round an uncertainty to 1-2 significant digits and its value to the same place."""
    value, uncertainty = float(value), abs(float(uncertainty))
    if not math.isfinite(uncertainty) or uncertainty == 0:
        return f"{value:g}", f"{uncertainty:g}"
    exponent = math.floor(math.log10(uncertainty))
    leading = uncertainty / 10 ** exponent
    significant_digits = 2 if leading < 3 else 1
    place = exponent - significant_digits + 1
    rounded_value = round(value, -place)
    rounded_uncertainty = round(uncertainty, -place)
    decimals = max(0, -place)
    return f"{rounded_value:.{decimals}f}", f"{rounded_uncertainty:.{decimals}f}"


def result_summary(result, data, irradiation_wavelength_nm):
    values = result.yield_fit.values * 100
    errors = result.yield_errors * 100
    pss = result.yield_fit.concentrations[-1]
    pss = pss / pss.sum() * 100
    formatted = [format_value_uncertainty(value, error)
                 for value, error in zip(values, errors)]
    return {
        "schema_version": 1,
        "fit_method": result.fit_method,
        "quantum_yield_percent": {"R_to_P": float(values[0]), "P_to_R": float(values[1])},
        "quantum_yield_error_percent": {"R_to_P": float(errors[0]), "P_to_R": float(errors[1])},
        "quantum_yield_formatted_percent": {
            "R_to_P": {"value": formatted[0][0], "error": formatted[0][1]},
            "P_to_R": {"value": formatted[1][0], "error": formatted[1][1]},
        },
        "photostationary_state_percent": {
            "reactant": float(pss[0]), "product": float(pss[1])
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


def write_results(path, result, data, irradiation_wavelength_nm):
    summary = result_summary(result, data, irradiation_wavelength_nm)
    formatted = summary["quantum_yield_formatted_percent"]
    pss = summary["photostationary_state_percent"]
    low, high = data.wavelength_limits
    method = "Concentrations" if result.fit_method == "concentrations" else "Emission"
    text = f"""PSS_Reactant (%): {pss['reactant']:.1f}
PSS_Product (%): {pss['product']:.1f}
QY_AB_opt (%): {formatted['R_to_P']['value']}
QY_BA_opt (%): {formatted['P_to_R']['value']}
error_QY_AB (%): {formatted['R_to_P']['error']}
error_QY_BA (%): {formatted['P_to_R']['error']}

Volume (ml): {data.volume_ml:g}
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
    traces = pd.DataFrame({
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
    })
    traces_path = stem.parent / f"{stem.name}_traces.tsv"
    traces.to_csv(traces_path, sep="\t", index=False)

    measured_absorbance = result.absorbance.T
    concentration_reconstruction = result.concentration_fit.reconstructed_absorbance
    kinetic_fit = (fitted @ np.vstack((result.epsilon_r, result.epsilon_p))
                   * data.path_length_cm)
    count_time, count_wavelength = measured_absorbance.shape
    spectra = pd.DataFrame({
        "time_s": np.repeat(data.timestamps, count_wavelength),
        "wavelength_nm": np.tile(result.wavelengths, count_time),
        "absorbance_measured": measured_absorbance.ravel(),
        "absorbance_concentration_reconstruction": concentration_reconstruction.ravel(),
        "absorbance_kinetic_fit": kinetic_fit.ravel(),
        "concentration_reconstruction_residual": (
            measured_absorbance - concentration_reconstruction).ravel(),
        "kinetic_fit_residual": (measured_absorbance - kinetic_fit).ravel(),
    })
    spectra_path = stem.parent / f"{stem.name}_spectra.tsv"
    spectra.to_csv(spectra_path, sep="\t", index=False)
    return traces_path, spectra_path
