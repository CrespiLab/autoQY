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
    def scenarios(self):
        return self.lower, self.nominal, self.upper


@dataclass(frozen=True)
class EpsilonUncertaintySummary:
    method: str
    error_metric: str
    scenario_count: int
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
        negative_points = 0
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


def run_with_epsilon_uncertainty(data, reactant, product, error_metric="sd"):
    """Run all low/mean/high reactant-product combinations through the pipeline."""
    from .pipeline import run_analysis_pipeline

    nominal_data = replace(
        data,
        epsilon_r=(reactant.wavelengths, reactant.nominal),
        epsilon_p=(product.wavelengths, product.nominal),
    )
    nominal_result = run_analysis_pipeline(nominal_data)
    results = []
    for reactant_index, reactant_curve in enumerate(reactant.scenarios):
        for product_index, product_curve in enumerate(product.scenarios):
            if reactant_index == 1 and product_index == 1:
                result = nominal_result
            else:
                scenario_data = replace(
                    data,
                    epsilon_r=(reactant.wavelengths, reactant_curve),
                    epsilon_p=(product.wavelengths, product_curve),
                )
                result = run_analysis_pipeline(scenario_data)
            results.append(result)

    values = np.asarray([result.yield_fit.values for result in results])
    scenario_errors = np.asarray([result.yield_errors for result in results])
    nominal_values = nominal_result.yield_fit.values
    epsilon_minimum = np.min(values, axis=0)
    epsilon_maximum = np.max(values, axis=0)
    epsilon_errors = np.maximum(
        nominal_values - epsilon_minimum, epsilon_maximum - nominal_values
    )
    combined_lower = np.min(values - scenario_errors, axis=0)
    combined_upper = np.max(values + scenario_errors, axis=0)
    combined_errors = np.maximum(
        nominal_values - combined_lower, combined_upper - nominal_values
    )
    summary = EpsilonUncertaintySummary(
        method="deterministic_extremes",
        error_metric=str(error_metric).lower(),
        scenario_count=len(results),
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
    )
    return replace(
        nominal_result,
        yield_errors=combined_errors,
        epsilon_uncertainty=summary,
    )


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
