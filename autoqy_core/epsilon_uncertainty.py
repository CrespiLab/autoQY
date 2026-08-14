"""Deterministic propagation of wavelength-resolved absorptivity errors."""

from dataclasses import dataclass, replace
from pathlib import Path
import warnings

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class EpsilonEnvelope:
    wavelengths: np.ndarray
    nominal: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    source_schema: str
    error_metric: str
    source_path: str
    constrained_negative_points: int = 0

    @property
    def bounds(self):
        return self.lower, self.nominal, self.upper


@dataclass(frozen=True)
class EpsilonUncertaintySummary:
    method: str
    error_metric: str
    bound_combination_count: int
    reactant_source_schema: str
    product_source_schema: str
    reactant_source_path: str
    product_source_path: str
    reactant_error_metric: str
    product_error_metric: str
    constrained_negative_points: tuple[int, int]
    optimizer_power_errors: np.ndarray
    epsilon_errors: np.ndarray
    combined_errors: np.ndarray
    epsilon_yield_minimum: np.ndarray
    epsilon_yield_maximum: np.ndarray
    concentration_data_minimum: np.ndarray
    concentration_data_maximum: np.ndarray
    concentration_fit_minimum: np.ndarray
    concentration_fit_maximum: np.ndarray
    fraction_residual_minimum: np.ndarray
    fraction_residual_maximum: np.ndarray
    absorbance_residual_rmse_minimum: float
    absorbance_residual_rmse_maximum: float
    bound_labels: tuple[str, ...]
    bound_yields: np.ndarray
    bound_optimizer_errors: np.ndarray
    bound_fraction_rmse: np.ndarray
    bound_absorbance_rmse: np.ndarray
    bound_active_bounds: np.ndarray
    bound_jacobian_conditions: np.ndarray


def load_epsilon_envelope(path, error_metric="sd"):
    """Load nominal and wavelength-resolved bounds from Spectral Treatment TSV."""
    path = Path(path)
    frame = pd.read_csv(path, sep="\t")
    wavelengths = _column(frame, "Wavelength_nm", path)
    metric = str(error_metric).lower()
    if metric not in {"sd", "sem"}:
        raise ValueError("epsilon error_metric must be sd or sem")

    if "Epsilon_mean_M-1_cm-1" in frame:
        nominal = _column(frame, "Epsilon_mean_M-1_cm-1", path)
        error_column = {
            "sd": "Epsilon_SD_M-1_cm-1",
            "sem": "Epsilon_SEM_M-1_cm-1",
        }[metric]
        error = np.abs(_column(frame, error_column, path))
        if metric == "sd" and {
            "Epsilon_lower_nonnegative_M-1_cm-1",
            "Epsilon_upper_nonnegative_M-1_cm-1",
        }.issubset(frame.columns):
            lower = _column(frame, "Epsilon_lower_nonnegative_M-1_cm-1", path)
            upper = _column(frame, "Epsilon_upper_nonnegative_M-1_cm-1", path)
        else:
            lower = np.maximum(nominal - error, 0.0)
            upper = np.maximum(nominal + error, 0.0)
        source_schema = "epsilon_spectra"
        negative_points = int(np.count_nonzero(nominal < 0))
        if negative_points:
            warnings.warn(
                f"Negative mean epsilon values in {path} were constrained to zero",
                RuntimeWarning,
                stacklevel=2,
            )
            nominal = np.maximum(nominal, 0.0)
    elif "Product_epsilon_M-1_cm-1" in frame:
        nominal = _column(frame, "Product_epsilon_M-1_cm-1", path)
        if np.min(nominal) < -500:
            raise ValueError(
                f"NMR-derived product epsilon in {path} falls below -500 M-1 cm-1"
            )
        negative_points = int(np.count_nonzero(nominal < 0))
        if negative_points:
            warnings.warn(
                f"Negative NMR-derived epsilon values in {path} were constrained to zero",
                RuntimeWarning,
                stacklevel=2,
            )
            nominal = np.maximum(nominal, 0.0)
        lower = _column(frame, "Product_lower_nonnegative_M-1_cm-1", path)
        upper = _column(frame, "Product_upper_nonnegative_M-1_cm-1", path)
        source_schema = "nmr_derived_product"
        metric = "asymmetric_bounds"
    else:
        raise ValueError(
            f"{path} is not an AutoQY epsilon or NMR-derived product TSV"
        )

    _validate_envelope(path, wavelengths, nominal, lower, upper)
    return EpsilonEnvelope(
        wavelengths, nominal, lower, upper, source_schema, metric,
        str(path.resolve()), negative_points,
    )


def load_epsilon_nominal(path):
    """Return the nominal curve from an AutoQY epsilon TSV, or None otherwise."""
    path = Path(path)
    try:
        columns = set(pd.read_csv(path, sep="\t", nrows=0).columns)
    except (OSError, UnicodeError, pd.errors.ParserError):
        return None
    if not ({"Wavelength_nm", "Epsilon_mean_M-1_cm-1"}.issubset(columns) or
            {"Wavelength_nm", "Product_epsilon_M-1_cm-1"}.issubset(columns)):
        return None
    envelope = load_epsilon_envelope(path, "sd")
    return envelope.wavelengths, envelope.nominal


def run_with_epsilon_uncertainty(data, reactant, product, error_metric="sd"):
    """Run distinct low/mean/high epsilon-bound combinations through the pipeline."""
    from .pipeline import run_analysis_pipeline

    nominal_data = replace(
        data,
        epsilon_r=(reactant.wavelengths, reactant.nominal),
        epsilon_p=(product.wavelengths, product.nominal),
    )
    nominal_result = run_analysis_pipeline(nominal_data)
    bound_pairs = []
    level_names = ("lower", "nominal", "upper")
    for reactant_level, reactant_curve in zip(level_names, reactant.bounds):
        for product_level, product_curve in zip(level_names, product.bounds):
            if not any(np.array_equal(reactant_curve, existing_reactant) and
                       np.array_equal(product_curve, existing_product)
                       for existing_reactant, existing_product, _ in bound_pairs):
                label = f"reactant {reactant_level} / product {product_level}"
                bound_pairs.append((reactant_curve, product_curve, label))
    results = []
    for reactant_curve, product_curve, _ in bound_pairs:
        if (np.array_equal(reactant_curve, reactant.nominal) and
                np.array_equal(product_curve, product.nominal)):
            result = nominal_result
        else:
            bound_data = replace(
                data,
                epsilon_r=(reactant.wavelengths, reactant_curve),
                epsilon_p=(product.wavelengths, product_curve),
            )
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", message="Initial product is .*", category=RuntimeWarning
                )
                result = run_analysis_pipeline(bound_data)
        results.append(result)

    values = np.asarray([result.yield_fit.values for result in results])
    bound_errors = np.asarray([result.yield_errors for result in results])
    nominal_values = nominal_result.yield_fit.values
    epsilon_minimum = np.min(values, axis=0)
    epsilon_maximum = np.max(values, axis=0)
    epsilon_errors = np.maximum(
        nominal_values - epsilon_minimum, epsilon_maximum - nominal_values
    )
    combined_lower = np.min(values - bound_errors, axis=0)
    combined_upper = np.max(values + bound_errors, axis=0)
    combined_errors = np.maximum(
        nominal_values - combined_lower, combined_upper - nominal_values
    )
    concentration_data = np.asarray([
        result.concentration_fit.concentrations for result in results
    ])
    concentration_fits = np.asarray([
        result.yield_fit.concentrations for result in results
    ])
    fitted_totals = concentration_fits.sum(axis=2)
    fitted_fractions = np.divide(
        concentration_fits[:, :, 0], fitted_totals,
        out=np.zeros_like(fitted_totals), where=fitted_totals != 0,
    )
    fraction_residuals = np.asarray([
        result.concentration_fit.fractions[:, 0] for result in results
    ]) - fitted_fractions
    residual_rmse = np.asarray([_absorbance_residual_rmse(result, data.path_length_cm)
                                for result in results])
    bound_fraction_rmse = np.sqrt(np.mean(fraction_residuals ** 2, axis=1))
    summary = EpsilonUncertaintySummary(
        method="deterministic_extremes",
        error_metric=str(error_metric).lower(),
        bound_combination_count=len(results),
        reactant_source_schema=reactant.source_schema,
        product_source_schema=product.source_schema,
        reactant_source_path=reactant.source_path,
        product_source_path=product.source_path,
        reactant_error_metric=reactant.error_metric,
        product_error_metric=product.error_metric,
        constrained_negative_points=(
            reactant.constrained_negative_points,
            product.constrained_negative_points,
        ),
        optimizer_power_errors=nominal_result.yield_errors,
        epsilon_errors=epsilon_errors,
        combined_errors=combined_errors,
        epsilon_yield_minimum=epsilon_minimum,
        epsilon_yield_maximum=epsilon_maximum,
        concentration_data_minimum=np.min(concentration_data, axis=0),
        concentration_data_maximum=np.max(concentration_data, axis=0),
        concentration_fit_minimum=np.min(concentration_fits, axis=0),
        concentration_fit_maximum=np.max(concentration_fits, axis=0),
        fraction_residual_minimum=np.min(fraction_residuals, axis=0),
        fraction_residual_maximum=np.max(fraction_residuals, axis=0),
        absorbance_residual_rmse_minimum=float(np.min(residual_rmse)),
        absorbance_residual_rmse_maximum=float(np.max(residual_rmse)),
        bound_labels=tuple(item[2] for item in bound_pairs),
        bound_yields=values,
        bound_optimizer_errors=bound_errors,
        bound_fraction_rmse=bound_fraction_rmse,
        bound_absorbance_rmse=residual_rmse,
        bound_active_bounds=np.asarray([
            result.yield_fit.active_bounds for result in results
        ], dtype=bool),
        bound_jacobian_conditions=np.asarray([
            result.yield_fit.jacobian_condition for result in results
        ], dtype=float),
    )
    return replace(
        nominal_result,
        yield_errors=combined_errors,
        epsilon_uncertainty=summary,
    )


def _absorbance_residual_rmse(result, path_length_cm):
    epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
    fitted = result.yield_fit.concentrations @ epsilon * path_length_cm
    if result.yield_fit.absorbance_correction is not None:
        fitted = fitted + result.yield_fit.absorbance_correction
    residual = result.absorbance.T - fitted
    return float(np.sqrt(np.mean(residual ** 2)))


def _column(frame, name, path):
    if name not in frame:
        raise ValueError(f"Missing {name} in {path}")
    values = pd.to_numeric(frame[name], errors="raise").to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError(f"{name} in {path} must contain a finite value at every wavelength")
    return values


def _validate_envelope(path, wavelengths, nominal, lower, upper):
    if len(wavelengths) < 2 or any(
            values.shape != wavelengths.shape for values in (nominal, lower, upper)):
        raise ValueError(f"Epsilon columns in {path} must have the same non-trivial length")
    if np.any(np.diff(wavelengths) <= 0):
        raise ValueError(f"Wavelengths in {path} must increase strictly")
    if np.any(lower < 0) or np.any(upper < 0):
        raise ValueError(f"Epsilon uncertainty bounds in {path} must be non-negative")
    tolerance = np.finfo(float).eps * np.maximum(1.0, np.abs(nominal)) * 16
    if np.any(lower > nominal + tolerance) or np.any(upper < nominal - tolerance):
        raise ValueError(f"Epsilon bounds in {path} must enclose the nominal spectrum")
