"""Run AutoQY analyses from validated configurations."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .config import AnalysisConfig, input_format, load_config, validate_config
from .io import load_spectra, load_spectrum, load_timestamps
from .output import result_summary, write_detailed_data, write_results
from .pipeline import AnalysisInput, run_analysis_pipeline
from .plotting import write_figure


@dataclass(frozen=True)
class RunOutput:
    result: object
    files: tuple[Path, ...]


def run_analysis(config, output_directory=None):
    if not isinstance(config, AnalysisConfig):
        config = load_config(config)
    validate_config(config)
    values = config.values
    inputs, experiment = values["inputs"], values["experiment"]
    processing, fit = values["processing"], values["fit"]

    wavelengths, absorbance = load_spectra(
        config.input_path("measurement_spectra"), input_format(config, "measurement_spectra")
    )
    epsilon_r = load_spectrum(
        config.input_path("reactant_absorptivity"), input_format(config, "reactant_absorptivity")
    )
    epsilon_p = load_spectrum(
        config.input_path("product_absorptivity"), input_format(config, "product_absorptivity")
    )
    led = load_spectrum(config.input_path("led_emission"), input_format(config, "led_emission"))
    timestamps = load_timestamps(
        config.input_path("timestamps"), input_format(config, "timestamps")
    )
    _validate_loaded_data(wavelengths, absorbance, timestamps, epsilon_r, epsilon_p)

    power_mw, power_error_mw = experiment["power_mw"], experiment["power_error_mw"]

    smoothing, baseline = processing["led_smoothing"], processing["led_baseline"]
    initial, bounds = fit["initial_quantum_yields"], fit["quantum_yield_bounds"]
    data = AnalysisInput(
        wavelengths=wavelengths,
        absorbance=absorbance,
        timestamps=timestamps,
        epsilon_r=epsilon_r,
        epsilon_p=epsilon_p,
        led=led,
        power_mw=power_mw,
        power_error_mw=power_error_mw,
        volume_ml=experiment["volume_ul"] / 1000,
        thermal_rate=experiment["thermal_back_reaction_s_1"],
        path_length_cm=experiment["path_length_cm"],
        wavelength_limits=tuple(processing["wavelength_range_nm"]),
        baseline_correct_led=baseline["enabled"],
        led_smoothing_window=smoothing["window_points"],
        led_polynomial_order=smoothing["polynomial_order"],
        baseline_exclusion_fwhm_multiplier=baseline["exclusion_fwhm_multiplier"],
        fit_method=fit["method"],
        emission_threshold_fraction=fit.get("emission_threshold_fraction", 0.01),
        initial_yields=(initial["R_to_P"], initial["P_to_R"]),
        yield_bounds=(bounds["minimum"], bounds["maximum"]),
    )
    result = run_analysis_pipeline(data)
    files = _write_outputs(config, data, result, output_directory)
    return RunOutput(result, files)


def _write_outputs(config, data, result, output_directory):
    values, output = config.values, config.values["outputs"]
    directory = Path(output_directory) if output_directory else config.resolve(output["directory"])
    stem = directory / output["stem"]
    paths = []
    if output["write_text"]:
        paths.append(stem.with_suffix(".txt"))
    if output["write_figures"]:
        paths.extend((stem.with_suffix(".png"), stem.with_suffix(".svg")))
    if output["write_json"]:
        paths.append(directory / f"{output['stem']}_results.json")
    if output["write_config"]:
        paths.append(directory / f"{output['stem']}_config.json")
    if output.get("write_detailed_data", False):
        paths.extend((directory / f"{output['stem']}_traces.tsv",
                      directory / f"{output['stem']}_spectra.tsv"))
    existing = [path for path in paths if path.exists()]
    if existing and not output["overwrite"]:
        raise FileExistsError("Output exists and overwrite is false: " + ", ".join(map(str, existing)))

    directory.mkdir(parents=True, exist_ok=True)
    wavelength = values["experiment"]["irradiation_wavelength_nm"]
    if output["write_text"]:
        write_results(stem.with_suffix(".txt"), result, data, wavelength)
    if output["write_figures"]:
        write_figure(stem, result, data, values["plots"]["absorbance_residual_percentile"])
    if output["write_json"]:
        summary = result_summary(result, data, wavelength)
        (directory / f"{output['stem']}_results.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
    if output["write_config"]:
        (directory / f"{output['stem']}_config.json").write_text(
            json.dumps(values, indent=2) + "\n", encoding="utf-8"
        )
    if output.get("write_detailed_data", False):
        write_detailed_data(stem, result, data)
    return tuple(paths)


def _validate_loaded_data(wavelengths, absorbance, timestamps, epsilon_r, epsilon_p):
    if absorbance.shape[1] != len(timestamps):
        raise ValueError(f"Found {absorbance.shape[1]} spectra but {len(timestamps)} timestamps")
    if np.any(np.diff(wavelengths) <= 0) or np.any(np.diff(timestamps) < 0):
        raise ValueError("Wavelengths must increase and timestamps must not decrease")
    low = max(epsilon_r[0][0], epsilon_p[0][0])
    high = min(epsilon_r[0][-1], epsilon_p[0][-1])
    if low >= high or wavelengths[-1] < low or wavelengths[0] > high:
        raise ValueError("Measurement and molar-absorptivity spectra do not overlap")
