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
    rounded_uncertainty = round(uncertainty, -place)
    # Rounding can carry the uncertainty into the next decade (0.099 -> 0.1).
    # Recalculate the reporting place so the value still matches the displayed
    # uncertainty rather than retaining spurious decimal places.
    rounded_exponent = math.floor(math.log10(rounded_uncertainty))
    place = rounded_exponent - significant_digits + 1
    rounded_uncertainty = round(rounded_uncertainty, -place)
    rounded_value = round(value, -place)
    decimals = max(0, -place)
    return f"{rounded_value:.{decimals}f}", f"{rounded_uncertainty:.{decimals}f}"


def format_scientific_value_uncertainty(value, uncertainty, two_digit_threshold=3):
    """Format a value/error pair with one shared base-ten exponent."""
    value, uncertainty = float(value), abs(float(uncertainty))
    reference = abs(value) if value else uncertainty
    if not math.isfinite(reference) or reference == 0:
        formatted = format_value_uncertainty(value, uncertainty, two_digit_threshold)
        return formatted[0], formatted[1], 0
    exponent = math.floor(math.log10(reference))
    scale = 10 ** exponent
    formatted = format_value_uncertainty(
        value / scale, uncertainty / scale, two_digit_threshold
    )
    return formatted[0], formatted[1], exponent


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
            "irradiation_source": ("chemical_actinometer"
                                   if data.photon_flux_mol_s is not None
                                   else "optical_power"),
            "power_mw": data.power_mw,
            "power_error_mw": data.power_error_mw,
            "photon_flux_mol_s": data.photon_flux_mol_s,
            "photon_flux_error_mol_s": data.photon_flux_error_mol_s,
            "thermal_back_reaction_s_1": data.thermal_rate,
            "thermal_forward_reaction_s_1": getattr(data, "thermal_forward_rate", 0),
            "irradiation_wavelength_nm": irradiation_wavelength_nm,
            "path_length_cm": data.path_length_cm,
        },
    }
    diagnostic = getattr(result, "ab_model_diagnostic", None)
    if diagnostic is not None:
        apparent = not diagnostic.model_supported
        summary["ab_model_assessment"] = {
            "status": diagnostic.level,
            "model_supported": diagnostic.model_supported,
            "tracked_balance_change_percent": diagnostic.balance_change_fraction * 100,
            "tracked_balance_span_percent": diagnostic.balance_span_fraction * 100,
            "spectral_relative_rmse_percent": diagnostic.spectral_relative_rmse * 100,
            "quantum_yield_interpretation": (
                "apparent estimate; degradation or another side process breaks the closed "
                "A <=> B model" if apparent else "closed A <=> B model supported"
            ),
        }
    nipe = getattr(result.yield_fit, "nipe", None)
    if nipe is not None:
        summary["nipe"] = {
            "reference": "Vorobyev, Lim and Lee, J. Photochem. Photobiol. A 478 (2026) 117229",
            "interval_count": nipe.interval_count,
            "minimum_interval_fraction": nipe.minimum_interval_fraction,
            "normalized_residual_rmse": nipe.normalized_residual_rmse,
        }
        windows = nipe.window_analysis
        if windows is not None:
            summary["nipe"]["pre_plateau_window_analysis"] = {
                "plateau_detected": windows.plateau_detected,
                "plateau_time_s": windows.plateau_time_s,
                "analysis_end_time_s": windows.analysis_end_time_s,
                "window_point_count": windows.window_point_count,
                "window_duration_s": windows.window_duration_s,
                "window_count": len(windows.window_start_s),
                "extrapolated_zero_exposure_yield_percent": {
                    "R_to_P": float(windows.extrapolated_values[0] * 100),
                    "P_to_R": float(windows.extrapolated_values[1] * 100),
                },
                "extrapolated_standard_error_percent": {
                    "R_to_P": float(windows.extrapolated_standard_errors[0] * 100),
                    "P_to_R": float(windows.extrapolated_standard_errors[1] * 100),
                },
                "full_trace_change_percent": {
                    "R_to_P": float(windows.full_trace_change_fraction[0] * 100),
                    "P_to_R": float(windows.full_trace_change_fraction[1] * 100),
                },
                "windows": [
                    {
                        "start_s": float(start),
                        "end_s": float(end),
                        "midpoint_s": float(midpoint),
                        "R_to_P_percent": float(values[0] * 100),
                        "P_to_R_percent": float(values[1] * 100),
                        "R_to_P_standard_error_percent": float(errors[0] * 100),
                        "P_to_R_standard_error_percent": float(errors[1] * 100),
                        "jacobian_condition": float(condition),
                    }
                    for start, end, midpoint, values, errors, condition in zip(
                        windows.window_start_s, windows.window_end_s,
                        windows.window_midpoint_s, windows.window_values,
                        windows.window_standard_errors,
                        windows.window_jacobian_conditions,
                    )
                ],
            }
    if data.photon_flux_mol_s is not None:
        flux = format_scientific_value_uncertainty(
            data.photon_flux_mol_s, data.photon_flux_error_mol_s
        )
        summary["experiment"]["photon_flux_formatted_mol_s"] = {
            "value": flux[0], "error": flux[1], "exponent": flux[2]
        }
    uncertainty = result.epsilon_uncertainty
    if uncertainty is not None:
        optimizer_power = uncertainty.optimizer_power_errors * 100
        epsilon = uncertainty.epsilon_errors * 100
        combined = uncertainty.combined_errors * 100
        irradiation_error_key = ("optimizer_and_photon_flux"
                                 if data.photon_flux_mol_s is not None
                                 else "optimizer_and_power")
        summary["quantum_yield_error_components_percent"] = {
            irradiation_error_key: _yield_pair(optimizer_power),
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
        "concentrations": "Concentrations (legacy pure NNLS)",
        "emission": "Emission (legacy)",
        "regularized_concentrations": "Regularized concentrations",
        "ode_absorbance": "Full-spectrum ODE absorbance",
        "nipe": "NIPE normalized integrated photokinetic equation",
    }[result.fit_method]
    assessment_text = ""
    assessment = summary.get("ab_model_assessment")
    if assessment is not None:
        assessment_text = f"""A<=>B model status: {assessment['status'].upper()}
A<=>B model supported: {'YES' if assessment['model_supported'] else 'NO'}
Tracked A+B balance change (%): {assessment['tracked_balance_change_percent']:.3g}
Tracked A+B balance span (%): {assessment['tracked_balance_span_percent']:.3g}
Two-reference spectral relative RMSE (%): {assessment['spectral_relative_rmse_percent']:.3g}
Quantum-yield interpretation: {assessment['quantum_yield_interpretation']}

"""
    nipe_text = ""
    if "nipe" in summary:
        nipe = summary["nipe"]
        nipe_text = f"""NIPE normalized intervals: {nipe['interval_count']}
NIPE normalized residual RMSE: {nipe['normalized_residual_rmse']:.6g}

"""
        windows = nipe.get("pre_plateau_window_analysis")
        if windows is not None:
            early = windows["extrapolated_zero_exposure_yield_percent"]
            early_error = windows["extrapolated_standard_error_percent"]
            early_rp = format_value_uncertainty(
                early["R_to_P"], early_error["R_to_P"], two_digit_threshold=2,
            )
            early_pr = format_value_uncertainty(
                early["P_to_R"], early_error["P_to_R"], two_digit_threshold=2,
            )
            change = windows["full_trace_change_percent"]
            nipe_text += f"""NIPE pre-plateau window count: {windows['window_count']}
NIPE automatic plateau time (s): {windows['plateau_time_s']:.6g}
NIPE last pre-plateau time (s): {windows['analysis_end_time_s']:.6g}
NIPE window duration (s): {windows['window_duration_s']:.6g}
NIPE zero-exposure apparent QY R_to_P (%): {early_rp[0]} +/- {early_rp[1]}
NIPE zero-exposure apparent QY P_to_R (%): {early_pr[0]} +/- {early_pr[1]}
NIPE full-trace change from early R_to_P (%): {change['R_to_P']:.6g}
NIPE full-trace change from early P_to_R (%): {change['P_to_R']:.6g}

"""
    epsilon_text = "\n"
    if result.epsilon_uncertainty is not None:
        components = summary["quantum_yield_error_components_percent"]
        metadata = summary["epsilon_uncertainty"]
        irradiation_error_key = ("optimizer_and_photon_flux"
                                 if data.photon_flux_mol_s is not None
                                 else "optimizer_and_power")
        irradiation_error_label = ("photon flux" if data.photon_flux_mol_s is not None
                                   else "power")
        epsilon_text = f"""Error component optimizer + {irradiation_error_label} R_to_P (%): {components[irradiation_error_key]['R_to_P']:g}
Error component optimizer + {irradiation_error_label} P_to_R (%): {components[irradiation_error_key]['P_to_R']:g}
Error component epsilon R_to_P (%): {components['epsilon']['R_to_P']:g}
Error component epsilon P_to_R (%): {components['epsilon']['P_to_R']:g}
Epsilon uncertainty method: {metadata['method']}
Reactant epsilon error metric: {metadata['reactant_error_metric']}
Product epsilon error metric: {metadata['product_error_metric']}
Epsilon bound combinations: {metadata['bound_combination_count']}
Reactant epsilon values constrained to zero: {metadata['constrained_negative_points']['reactant']}
NMR epsilon values constrained to zero: {metadata['constrained_negative_points']['product']}

"""
    if data.photon_flux_mol_s is None:
        irradiation_input = (
            f"Power average (mW): {data.power_mw:g}\n"
            f"Power error (mW): {data.power_error_mw:g}"
        )
    else:
        flux = summary["experiment"]["photon_flux_formatted_mol_s"]
        irradiation_input = (
            "Irradiation source: Chemical actinometer\n"
            f"Photon flux (mol photons/s): ({flux['value']} +/- {flux['error']}) "
            f"x 10^{flux['exponent']}"
        )
    text = f"""Composition at the last timestamp (s): {last_composition['time_s']:g}
Composition at the last timestamp - Reactant (%): {last_composition['reactant']:.1f}
Composition at the last timestamp - Product (%): {last_composition['product']:.1f}
Extrapolated PSS - Reactant (%): {extrapolated_pss['reactant']:.1f}
Extrapolated PSS - Product (%): {extrapolated_pss['product']:.1f}
QY_AB_opt (%): {formatted['R_to_P']['value']}
QY_BA_opt (%): {formatted['P_to_R']['value']}
error_QY_AB (%): {formatted['R_to_P']['error']}
error_QY_BA (%): {formatted['P_to_R']['error']}
{assessment_text}{nipe_text}{epsilon_text}Volume (ml): {data.volume_ml:g}
k thermal back-reaction (s-1): {data.thermal_rate:g}
k thermal forward-reaction (s-1): {getattr(data, 'thermal_forward_rate', 0):g}
{irradiation_input}
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
    traces_path = stem.parent / f"{stem.name}_traces.csv"
    traces.to_csv(traces_path, index=False)

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
    spectra_path = stem.parent / f"{stem.name}_spectra.csv"
    spectra.to_csv(spectra_path, index=False)
    return traces_path, spectra_path
