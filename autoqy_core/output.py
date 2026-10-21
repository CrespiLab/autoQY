"""Write AutoQY results without GUI dependencies."""

import math
from pathlib import Path

import numpy as np
import pandas as pd


def uncertainty_decimal_place(uncertainty, two_digit_threshold=3):
    """Decimal place of the last reported uncertainty digit (for a positive SD)."""
    uncertainty = abs(float(uncertainty))
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
    return place


def format_value_uncertainty(value, uncertainty, two_digit_threshold=3):
    """Round uncertainty and value to a scientifically matching decimal place."""
    value, uncertainty = float(value), abs(float(uncertainty))
    if not math.isfinite(uncertainty) or uncertainty == 0:
        return f"{value:g}", f"{uncertainty:g}"
    place = uncertainty_decimal_place(uncertainty, two_digit_threshold)
    rounded_uncertainty = round(uncertainty, -place)
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


def quantum_yield_values(result):
    """Use the MC mean for reporting, keeping the original fitted values on the result."""
    uncertainty = result.epsilon_uncertainty
    if getattr(uncertainty, "method", "") == "monte_carlo":
        return np.asarray(uncertainty.sample_mean, float)
    return np.asarray(result.yield_fit.values, float)


def formatted_quantum_yields(result):
    return [format_value_uncertainty(value * 100, error * 100)
            for value, error in zip(quantum_yield_values(result), result.yield_errors)]


def monte_carlo_messages(uncertainty):
    """Explain completion, reliability, and the next action in GUI and text exports."""
    u = uncertainty
    if u.failed_samples:
        messages = [f"{u.bound_combination_count} of {u.attempted_samples} simulations finished successfully.",
                    f"{u.failed_samples} fits failed, so the uncertainty estimate is incomplete. "
                    "Review the failed fits in Technical details before using these results."]
    else:
        messages = [f"All {u.attempted_samples} simulations finished successfully."]
    messages.extend(getattr(u, "assumed_input_warnings", ()))
    if getattr(u, "requested_workers", 1) > getattr(u, "workers", 1):
        messages.append(f"Using {u.workers} CPU workers instead of the requested {u.requested_workers} "
                        "to leave capacity for other applications.")
    if u.converged:
        messages.append("The mean and error passed the precision check.")
    elif not u.failed_samples:
        messages.append("The simulation limit was reached before the mean and error passed the precision check. "
                        "Their remaining simulation error is shown in Technical details.")
    if np.any(u.boundary_fractions):
        limits = "; ".join(f"{name}: {fraction * 100:.2g}%" for name, fraction in
                           zip(("R→P", "P→R"), u.boundary_fractions) if fraction > 0)
        messages.append(f"Some simulated yields reached an allowed minimum or maximum ({limits} of fits). "
                        "These limits affect the uncertainty range; review the yield limits and input data.")
    if getattr(u, "reused_samples", 0):
        messages.append(f"Reused {u.reused_samples} saved simulations; "
                        f"completed {u.attempted_samples - u.reused_samples} new simulations.")
    messages.append("The ± value is the estimated error (1 SD).")
    return messages


def result_summary(result, data, irradiation_wavelength_nm):
    values = quantum_yield_values(result) * 100
    errors = result.yield_errors * 100
    last_composition = result.yield_fit.concentrations[-1]
    last_composition = last_composition / last_composition.sum() * 100
    extrapolated_pss = result.extrapolated_pss
    extrapolated_pss = extrapolated_pss / extrapolated_pss.sum() * 100
    formatted = formatted_quantum_yields(result)
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
    if uncertainty is not None and uncertainty.method == "monte_carlo":
        summary.pop("quantum_yield_error_components_percent", None)
        summary["quantum_yield_estimate_kind"] = "monte_carlo_mean"
        summary["quantum_yield_error_kind"] = "Monte Carlo standard uncertainty (1 SD)"
        summary["quantum_yield_interval_percent"] = {
            "coverage_probability": uncertainty.coverage_probability,
            "endpoints_stable": uncertainty.interval_converged,
            "lower": _yield_pair(uncertainty.interval_lower * 100),
            "upper": _yield_pair(uncertainty.interval_upper * 100),
        }
        summary["monte_carlo"] = {
            "mean_quantum_yield_percent": _yield_pair(uncertainty.sample_mean * 100),
            "standard_uncertainty_percentage_points": _yield_pair(uncertainty.sample_standard_deviation * 100),
            "nominal_quantum_yield_percent": _yield_pair(result.yield_fit.values * 100),
            "requested_maximum_simulations": uncertainty.requested_samples,
            "workers": uncertainty.workers,
            "requested_workers": uncertainty.requested_workers,
            "available_cpus": uncertainty.available_cpus,
            "attempted_simulations": uncertainty.attempted_samples,
            "successful_simulations": uncertainty.bound_combination_count,
            "failed_simulations": uncertainty.failed_samples,
            "reused_simulations": uncertainty.reused_samples,
            "checkpoint_file": uncertainty.checkpoint_path,
            "boundary_hit_fraction": _yield_pair(uncertainty.boundary_fractions),
            "seed": uncertainty.random_seed,
            "mean_and_error_stable": uncertainty.converged,
            "endpoints_stable": uncertainty.interval_converged,
            "summary_mc_standard_error_percentage_points": {
                "mean": _yield_pair(uncertainty.summary_mc_error[0] * 100),
                "standard_uncertainty": _yield_pair(uncertainty.summary_mc_error[1] * 100),
            },
            "summary_mc_tolerance_percentage_points": _yield_pair(uncertainty.summary_mc_tolerance * 100),
            "endpoint_mc_standard_error_percentage_points": {
                "lower": _yield_pair(uncertainty.endpoint_mc_error[0] * 100),
                "upper": _yield_pair(uncertainty.endpoint_mc_error[1] * 100),
            },
            "scope": uncertainty.scope,
            "sampling_assumptions": list(uncertainty.sampling_assumptions),
            "assumed_input_warnings": list(uncertainty.assumed_input_warnings),
            "failures": list(uncertainty.failure_messages),
        }
        summary["epsilon_uncertainty"]["interval_kind"] = "Monte Carlo coverage interval"
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
        "ode_absorbance": "Full spectrum ODE",
    }[result.fit_method]
    uncertainty = result.epsilon_uncertainty
    uncertainty_text = "Error: fit and irradiation-input uncertainty.\n"
    if getattr(uncertainty, "method", "") == "monte_carlo":
        uncertainty_text = "\n".join(monte_carlo_messages(uncertainty)) + "\n"
        uncertainty_text += f"CPU workers: {uncertainty.workers}.\n"
        uncertainty_text += ("Included errors: reference spectra, NMR composition, irradiation input and kinetic fitting.\n"
                             "The PSS reconstruction absorbance is held fixed.\n")
        mc_error = uncertainty.summary_mc_error * 100
        uncertainty_text += (f"Interval endpoint precision check: {'passed' if uncertainty.interval_converged else 'not passed; interval limits are approximate'}.\n"
                             "Remaining simulation error (1 SD, percentage points): "
                             f"mean R->P {mc_error[0, 0]:.2g}, P->R {mc_error[0, 1]:.2g}; "
                             f"SD R->P {mc_error[1, 0]:.2g}, P->R {mc_error[1, 1]:.2g}.\n")
    elif uncertainty is not None:
        uncertainty_text = "Error: combined ε-range, fit, and irradiation-input sensitivity.\n"
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
Φ R→P: {formatted['R_to_P']['value']} ± {formatted['R_to_P']['error']}%
Φ P→R: {formatted['P_to_R']['value']} ± {formatted['P_to_R']['error']}%
{uncertainty_text}
Volume (ml): {data.volume_ml:g}
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
    if getattr(uncertainty, "method", "") == "monte_carlo":
        columns = {name.replace("_epsilon_min", "_mc_interval_lower").replace(
            "_epsilon_max", "_mc_interval_upper"): value for name, value in columns.items()}
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
    paths = [traces_path, spectra_path]
    if getattr(uncertainty, "method", "") == "monte_carlo":
        samples = pd.DataFrame({
            "successful_simulation": np.arange(1, len(uncertainty.bound_yields) + 1),
            "R_to_P_percent": uncertainty.bound_yields[:, 0] * 100,
            "P_to_R_percent": uncertainty.bound_yields[:, 1] * 100,
            "R_to_P_at_yield_limit": uncertainty.bound_active_bounds[:, 0],
            "P_to_R_at_yield_limit": uncertainty.bound_active_bounds[:, 1],
            "jacobian_condition": uncertainty.bound_jacobian_conditions,
            "absorbance_rmse": uncertainty.bound_absorbance_rmse,
        })
        samples_path = stem.parent / f"{stem.name}_uncertainty.csv"
        samples.to_csv(samples_path, index=False)
        paths.append(samples_path)
    return tuple(paths)
