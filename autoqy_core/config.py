"""Load and validate versioned AutoQY JSON configurations."""

from dataclasses import dataclass
import json
from pathlib import Path


class ConfigError(ValueError):
    pass


@dataclass(frozen=True)
class AnalysisConfig:
    values: dict
    base_directory: Path
    source: Path | None = None

    def resolve(self, value):
        path = Path(value)
        return path if path.is_absolute() else self.base_directory / path

    def input_path(self, name):
        return self.resolve(self.values["inputs"][name])


def load_config(path):
    path = Path(path).resolve()
    try:
        values = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ConfigError(f"Cannot load {path}: {error}") from error
    config = AnalysisConfig(values, path.parent, path)
    validate_config(config)
    return config


def input_format(config, name):
    inputs = config.values["inputs"]
    formats = inputs.get("formats", {})
    if name in formats:
        return formats[name]
    if name == "timestamps":
        return inputs.get("timestamp_format", "ahk_csv")
    return inputs.get("spectral_format", "spectragryph_tsv")


def validate_config(config):
    values = config.values
    errors = []
    if values.get("schema_version") != 1:
        errors.append("schema_version must be 1")

    required = {
        "analysis": ("id", "reactant_name", "product_name"),
        "inputs": ("measurement_spectra", "reactant_absorptivity",
                   "product_absorptivity", "timestamps"),
        "experiment": ("volume_ul", "thermal_back_reaction_s_1",
                       "irradiation_wavelength_nm", "path_length_cm"),
        "processing": ("wavelength_range_nm", "led_smoothing", "led_baseline"),
        "fit": ("method", "initial_quantum_yields", "quantum_yield_bounds"),
        "plots": ("absorbance_residual_percentile",),
        "outputs": ("directory", "stem", "write_text", "write_figures",
                    "write_json", "write_config", "overwrite"),
    }
    for section_name, keys in required.items():
        section = values.get(section_name)
        if not isinstance(section, dict):
            errors.append(f"{section_name} must be an object")
            continue
        errors.extend(f"{section_name}.{key} is required" for key in keys if key not in section)
    if errors:
        raise ConfigError("Invalid configuration:\n- " + "\n- ".join(errors))

    inputs = values["inputs"]
    experiment = values["experiment"]
    irradiation_source = experiment.get("irradiation_source", "optical_power")
    if irradiation_source not in {"optical_power", "chemical_actinometer"}:
        errors.append(
            "experiment.irradiation_source must be optical_power or chemical_actinometer"
        )
    if irradiation_source == "optical_power":
        for name in ("power_mw", "power_error_mw"):
            if name not in experiment:
                errors.append(f"experiment.{name} is required for optical_power")
        if "led_emission" not in inputs:
            errors.append("inputs.led_emission is required for optical_power")
    elif "photon_flux_mol_s" not in experiment:
        errors.append(
            "experiment.photon_flux_mol_s is required for chemical_actinometer"
        )

    formats = inputs.get("formats", {})
    if formats and not isinstance(formats, dict):
        errors.append("inputs.formats must be an object")
    for name in required["inputs"]:
        if not config.input_path(name).is_file():
            errors.append(f"inputs.{name} does not exist: {config.input_path(name)}")
        _validate_format(input_format(config, name), name, errors)
    if irradiation_source == "optical_power" and "led_emission" in inputs:
        if not config.input_path("led_emission").is_file():
            errors.append(
                f"inputs.led_emission does not exist: {config.input_path('led_emission')}"
            )
        _validate_format(input_format(config, "led_emission"), "led_emission", errors)

    uncertainty = values.get("uncertainty", {})
    if not isinstance(uncertainty, dict):
        errors.append("uncertainty must be an object")
    else:
        epsilon_uncertainty = uncertainty.get("epsilon", {})
        if not isinstance(epsilon_uncertainty, dict):
            errors.append("uncertainty.epsilon must be an object")
        else:
            epsilon_method = epsilon_uncertainty.get("method", "none")
            if epsilon_method not in {"none", "deterministic_extremes"}:
                errors.append(
                    "uncertainty.epsilon.method must be none or deterministic_extremes"
                )
            if epsilon_uncertainty.get("error_metric", "sd") not in {"sd", "sem"}:
                errors.append("uncertainty.epsilon.error_metric must be sd or sem")

    for name in ("volume_ul", "path_length_cm", "irradiation_wavelength_nm"):
        if not _positive(experiment[name]):
            errors.append(f"experiment.{name} must be positive")
    if irradiation_source == "optical_power" and all(
            name in experiment for name in ("power_mw", "power_error_mw")):
        if not _positive(experiment["power_mw"]):
            errors.append("experiment.power_mw must be positive")
        if not _nonnegative(experiment["power_error_mw"]):
            errors.append("experiment.power_error_mw must be nonnegative")
        elif (_positive(experiment["power_mw"])
              and experiment["power_error_mw"] >= experiment["power_mw"]):
            errors.append("experiment.power_error_mw must be smaller than power_mw")
    if irradiation_source == "chemical_actinometer" and "photon_flux_mol_s" in experiment:
        photon_flux = experiment["photon_flux_mol_s"]
        photon_flux_error = experiment.get("photon_flux_error_mol_s", 0)
        if not _positive(photon_flux):
            errors.append("experiment.photon_flux_mol_s must be positive")
        if not _nonnegative(photon_flux_error):
            errors.append("experiment.photon_flux_error_mol_s must be nonnegative")
        elif _positive(photon_flux) and photon_flux_error >= photon_flux:
            errors.append(
                "experiment.photon_flux_error_mol_s must be smaller than photon_flux_mol_s"
            )
    if not _nonnegative(experiment["thermal_back_reaction_s_1"]):
        errors.append("experiment.thermal_back_reaction_s_1 must be nonnegative")
    if not _nonnegative(experiment.get("thermal_forward_reaction_s_1", 0)):
        errors.append("experiment.thermal_forward_reaction_s_1 must be nonnegative")

    processing = values["processing"]
    wavelength_range = processing["wavelength_range_nm"]
    if not (isinstance(wavelength_range, list) and len(wavelength_range) == 2
            and all(isinstance(value, (int, float)) for value in wavelength_range)
            and wavelength_range[0] < wavelength_range[1]):
        errors.append("processing.wavelength_range_nm must contain two increasing values")
    smoothing = processing["led_smoothing"]
    if not isinstance(smoothing, dict) or smoothing.get("method") != "savgol":
        errors.append("processing.led_smoothing.method must be savgol")
    elif smoothing.get("window_points", 0) <= smoothing.get("polynomial_order", 0):
        errors.append("LED smoothing window_points must exceed polynomial_order")
    baseline = processing["led_baseline"]
    if not isinstance(baseline, dict) or not isinstance(baseline.get("enabled"), bool):
        errors.append("processing.led_baseline.enabled must be true or false")
    elif not _positive(baseline.get("exclusion_fwhm_multiplier")):
        errors.append("LED baseline exclusion_fwhm_multiplier must be positive")

    fit = values["fit"]
    methods = {"concentrations", "emission", "regularized_concentrations",
               "ode_absorbance"}
    if fit["method"] not in methods:
        errors.append("fit.method must be concentrations, emission, "
                      "regularized_concentrations, or ode_absorbance")
    threshold = fit.get("emission_threshold_fraction", 0.01)
    if not (_positive(threshold) and threshold < 1):
        errors.append("fit.emission_threshold_fraction must be greater than 0 and less than 1")
    if not _positive(fit.get("regularization_strength", 1)):
        errors.append("fit.regularization_strength must be positive")
    if fit.get("absorbance_baseline_order", 1) not in {-1, 0, 1}:
        errors.append("fit.absorbance_baseline_order must be -1, 0, or 1")
    if not _positive(fit.get("robust_loss_scale", 0.02)):
        errors.append("fit.robust_loss_scale must be positive")
    initial = fit["initial_quantum_yields"]
    bounds = fit["quantum_yield_bounds"]
    minimum, maximum = bounds.get("minimum"), bounds.get("maximum")
    if not (_nonnegative(minimum) and _positive(maximum) and minimum < maximum):
        errors.append("quantum-yield bounds must be increasing and nonnegative")
    elif any(not minimum <= initial.get(name, -1) <= maximum for name in ("R_to_P", "P_to_R")):
        errors.append("initial quantum yields must lie within their bounds")

    percentile = values["plots"]["absorbance_residual_percentile"]
    if not (_positive(percentile) and percentile <= 100):
        errors.append("plots.absorbance_residual_percentile must be greater than 0 and at most 100")
    output = values["outputs"]
    if not output["stem"] or Path(output["stem"]).name != output["stem"]:
        errors.append("outputs.stem must be a filename stem, not a path")
    for name in ("write_text", "write_figures", "write_json", "write_config",
                 "write_detailed_data", "overwrite"):
        if name in output and not isinstance(output[name], bool):
            errors.append(f"outputs.{name} must be true or false")

    if errors:
        raise ConfigError("Invalid configuration:\n- " + "\n- ".join(errors))
    return config


def _validate_format(value, name, errors):
    spec = {"type": value} if isinstance(value, str) else value
    if not isinstance(spec, dict):
        errors.append(f"format for {name} must be a string or object")
        return
    allowed = {"ahk_csv", "simple_csv", "generic_delimited"} if name == "timestamps" else {
        "spectragryph_tsv", "generic_delimited"
    }
    if spec.get("type") not in allowed:
        errors.append(f"unsupported format for {name}: {spec.get('type')}")
    if spec.get("type") == "generic_delimited":
        if not isinstance(spec.get("delimiter", ","), str):
            errors.append(f"generic delimiter for {name} must be text")
        if not isinstance(spec.get("header", True), bool):
            errors.append(f"generic header for {name} must be true or false")


def _positive(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0


def _nonnegative(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value >= 0
