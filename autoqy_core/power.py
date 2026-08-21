"""Headless processing of optical power-monitor traces."""

from dataclasses import asdict, dataclass
from io import StringIO
import json
from pathlib import Path

import numpy as np
import pandas as pd

REGION_NAMES = (
    "open_off_before", "open_on", "open_off_after",
    "cuvette_off_before", "cuvette_on", "cuvette_off_after",
)


@dataclass(frozen=True)
class PowerTraceResult:
    path: str
    power_mw: float
    repeatability_sd_mw: float
    standard_error_mw: float
    open_beam_power_mw: float
    cuvette_power_mw: float


@dataclass(frozen=True)
class PowerResult:
    power_mw: float
    power_error_mw: float
    uncertainty_output: str
    repeatability_sd_mw: float
    standard_error_mw: float
    between_replicate_sd_mw: float
    traces: tuple[PowerTraceResult, ...]

    def as_dict(self):
        result = asdict(self)
        result["traces"] = [asdict(trace) for trace in self.traces]
        return result


@dataclass(frozen=True)
class PowerBaseline:
    """Full baseline traces used by graphical and scripted clients."""

    values_mw: np.ndarray
    open_baseline_mw: np.ndarray
    open_corrected_mw: np.ndarray
    cuvette_baseline_mw: np.ndarray
    cuvette_corrected_mw: np.ndarray
    regions: dict


def run_power_analysis(config, output_path=None):
    """Run a standalone power analysis from a JSON path or dictionary."""
    if isinstance(config, (str, Path)):
        source = Path(config).resolve()
        values = json.loads(source.read_text(encoding="utf-8"))
        base_directory = source.parent
    else:
        values = dict(config)
        base_directory = Path.cwd()
    _validate_configuration(values, base_directory)
    result = process_power_configuration(values, base_directory)
    configured_output = values.get("output")
    if output_path is not None or configured_output:
        target = Path(output_path) if output_path is not None else Path(configured_output)
        if not target.is_absolute():
            target = base_directory / target
        if target.exists() and not values.get("overwrite", False):
            raise FileExistsError(f"Output exists and overwrite is false: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(result.as_dict(), indent=2) + "\n", encoding="utf-8")
    return result


def process_power_configuration(configuration, base_directory):
    """Process every configured OPM trace and combine replicate measurements."""
    degree = configuration.get("baseline_polynomial_degree", 3)
    traces = []
    for measurement in configuration["measurements"]:
        path = Path(measurement["path"])
        if not path.is_absolute():
            path = Path(base_directory) / path
        samples = (load_generic_power_csv(path)
                   if configuration["format"] == "generic_csv"
                   else load_thorlabs_opm(path))
        traces.append(process_power_trace(samples, measurement["regions"], degree, path))

    return combine_power_traces(
        traces, configuration.get("uncertainty_output", "repeatability_sd")
    )


def combine_power_traces(traces, uncertainty_output="repeatability_sd"):
    """Combine one to three independently processed power traces."""
    traces = tuple(traces)
    if not traces:
        raise ValueError("At least one processed power trace is required")
    if uncertainty_output not in {"repeatability_sd", "standard_error"}:
        raise ValueError("uncertainty_output must be repeatability_sd or standard_error")
    powers = np.array([trace.power_mw for trace in traces])
    repeatability = float(np.sqrt(np.mean(
        np.square([trace.repeatability_sd_mw for trace in traces])
    )))
    within_sem = np.array([trace.standard_error_mw for trace in traces])
    between_sd = float(np.std(powers, ddof=1)) if len(powers) > 1 else 0.0
    standard_error = float(np.sqrt(np.sum(within_sem ** 2) / len(traces) ** 2
                                   + between_sd ** 2 / len(traces)))
    selected = repeatability if uncertainty_output == "repeatability_sd" else standard_error
    return PowerResult(float(np.mean(powers)), selected, uncertainty_output,
                       repeatability, standard_error, between_sd, tuple(traces))


def load_thorlabs_opm(path):
    data = pd.read_csv(path, skiprows=14, float_precision="round_trip")
    return _power_values(data, path)


def load_thorlabs_opm_text(text, source="uploaded CSV"):
    """Read an uploaded Thorlabs OPM CSV without writing it to disk."""
    data = pd.read_csv(StringIO(text), skiprows=14, float_precision="round_trip")
    return _power_values(data, source)


def load_generic_power_csv(path):
    """Load a plain CSV containing a ``power_mw`` column."""
    data = pd.read_csv(path, float_precision="round_trip")
    if "power_mw" not in data:
        raise ValueError(f"Generic power CSV has no 'power_mw' column: {path}")
    values = pd.to_numeric(data["power_mw"], errors="raise").to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite power values found in {path}")
    return values


def _power_values(data, source):
    if "Power (W)" not in data:
        raise ValueError(f"Thorlabs OPM file has no 'Power (W)' column: {source}")
    values = pd.to_numeric(data["Power (W)"], errors="raise").to_numpy(float) * 1000
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite power values found in {source}")
    return values


def process_power_trace(power_mw, regions, baseline_degree=3, path=""):
    """Baseline one trace using two off regions around each on region."""
    baseline = baseline_power_trace(power_mw, regions, baseline_degree)
    parsed = baseline.regions
    open_values = baseline.open_corrected_mw[slice(*parsed["open_on"])]
    cuvette_values = baseline.cuvette_corrected_mw[slice(*parsed["cuvette_on"])]
    means = np.array([np.mean(open_values), np.mean(cuvette_values)])
    standard_deviations = np.array([
        np.std(open_values, ddof=1), np.std(cuvette_values, ddof=1)
    ])
    standard_errors = standard_deviations / np.sqrt(
        [len(open_values), len(cuvette_values)]
    )
    return PowerTraceResult(
        str(path), float(np.mean(means)),
        float(np.sqrt(np.mean(standard_deviations ** 2))),
        float(np.sqrt(np.sum(standard_errors ** 2)) / 2),
        float(means[0]), float(means[1]),
    )


def baseline_power_trace(power_mw, regions, baseline_degree=3):
    """Return the two fitted baselines and full corrected traces."""
    values = np.asarray(power_mw, float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Power data must be a finite one-dimensional sequence")
    if not isinstance(baseline_degree, int) or baseline_degree < 0:
        raise ValueError("baseline_degree must be a nonnegative integer")
    if not isinstance(regions, dict) or set(regions) != set(REGION_NAMES):
        raise ValueError("Exactly six named power regions are required")
    parsed = {name: _region(regions[name], len(values), name) for name in REGION_NAMES}
    x = np.arange(len(values), dtype=float)
    open_baseline = _fit_baseline(
        x, values, parsed["open_off_before"], parsed["open_off_after"],
        baseline_degree,
    )
    cuvette_baseline = _fit_baseline(
        x, values, parsed["cuvette_off_before"], parsed["cuvette_off_after"],
        baseline_degree,
    )
    return PowerBaseline(
        values, open_baseline, values - open_baseline,
        cuvette_baseline, values - cuvette_baseline, parsed,
    )


def _fit_baseline(x, values, off_before, off_after, degree):
    baseline_indices = np.r_[np.arange(*off_before), np.arange(*off_after)]
    if len(baseline_indices) <= degree:
        raise ValueError("Power baseline regions contain too few points")
    coefficients = np.polyfit(x[baseline_indices], values[baseline_indices], degree)
    return np.polyval(coefficients, x)


def _region(value, length, name):
    if not (isinstance(value, list) and len(value) == 2
            and all(isinstance(index, int) for index in value)
            and 0 <= value[0] < value[1] <= length):
        raise ValueError(f"Invalid power region {name}: {value}")
    return tuple(value)


def _validate_configuration(values, base_directory):
    if values.get("schema_version") != 1:
        raise ValueError("Power schema_version must be 1")
    if values.get("format") not in {"generic_csv", "thorlabs_opm_csv"}:
        raise ValueError("Power format must be generic_csv or thorlabs_opm_csv")
    degree = values.get("baseline_polynomial_degree", 3)
    if not isinstance(degree, int) or degree < 0:
        raise ValueError("baseline_polynomial_degree must be a nonnegative integer")
    if values.get("uncertainty_output", "repeatability_sd") not in {
            "repeatability_sd", "standard_error"}:
        raise ValueError("uncertainty_output must be repeatability_sd or standard_error")
    measurements = values.get("measurements")
    if not isinstance(measurements, list) or not measurements:
        raise ValueError("measurements must be a nonempty list")
    for index, measurement in enumerate(measurements):
        if not isinstance(measurement, dict) or "path" not in measurement:
            raise ValueError(f"Power measurement {index} requires a path")
        path = Path(measurement["path"])
        path = path if path.is_absolute() else Path(base_directory) / path
        if not path.is_file():
            raise ValueError(f"Power measurement does not exist: {path}")
        regions = measurement.get("regions")
        if not isinstance(regions, dict) or set(regions) != set(REGION_NAMES):
            raise ValueError(f"Power measurement {index} must define exactly six named regions")
