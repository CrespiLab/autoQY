"""Verify that generic CSV and Crespi-group example inputs are equivalent."""

from copy import deepcopy
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))

from autoqy_core.config import AnalysisConfig, load_config
from autoqy_core.power import run_power_analysis
from autoqy_core.runner import run_analysis


OUTPUT_FLAGS = (
    "write_text", "write_figures", "write_json", "write_config",
    "write_detailed_data",
)


def _analysis_result(path):
    config = load_config(path)
    values = deepcopy(config.values)
    for name in OUTPUT_FLAGS:
        values["outputs"][name] = False
    with TemporaryDirectory(prefix="autoqy-example-verification-") as temporary:
        return run_analysis(
            AnalysisConfig(values, config.base_directory, config.source), temporary
        ).result


def _assert_close(label, generic, crespi):
    if not np.allclose(generic, crespi, rtol=1e-11, atol=1e-13, equal_nan=True):
        difference = float(np.nanmax(np.abs(np.asarray(generic) - np.asarray(crespi))))
        raise AssertionError(f"{label} differs; maximum absolute difference = {difference:g}")


def verify_analysis_examples():
    checked = []
    for folder in sorted(ROOT.glob("Example-*")):
        generic_path = folder / "generic_inputs" / "analysis.json"
        crespi_path = folder / "crespi_group_inputs" / "analysis.json"
        if not generic_path.is_file() and not crespi_path.is_file():
            continue
        generic = _analysis_result(generic_path)
        crespi = _analysis_result(crespi_path)
        for name, generic_values, crespi_values in (
            ("quantum yields", generic.yield_fit.values, crespi.yield_fit.values),
            ("yield errors", generic.yield_errors, crespi.yield_errors),
            ("measured concentrations", generic.concentration_fit.concentrations,
             crespi.concentration_fit.concentrations),
            ("fitted concentrations", generic.yield_fit.concentrations,
             crespi.yield_fit.concentrations),
            ("extrapolated PSS", generic.extrapolated_pss, crespi.extrapolated_pss),
        ):
            _assert_close(f"{folder.name}: {name}", generic_values, crespi_values)
        checked.append(folder.name)
    return checked


def verify_power_example():
    folder = ROOT / "Example-Power"
    results = []
    for name in ("generic_inputs", "crespi_group_inputs"):
        config_path = folder / name / "power_analysis.json"
        document = json.loads(config_path.read_text(encoding="utf-8"))
        document.pop("output", None)
        for measurement in document["measurements"]:
            measurement["path"] = str((config_path.parent / measurement["path"]).resolve())
        results.append(run_power_analysis(document))
    generic, crespi = results
    _assert_close("Example-Power: power", generic.power_mw, crespi.power_mw)
    _assert_close("Example-Power: power error", generic.power_error_mw, crespi.power_error_mw)
    return folder.name


if __name__ == "__main__":
    examples = verify_analysis_examples()
    examples.append(verify_power_example())
    print(f"Equivalent generic/Crespi inputs verified for {len(examples)} examples:")
    print("\n".join(f"- {name}" for name in examples))
