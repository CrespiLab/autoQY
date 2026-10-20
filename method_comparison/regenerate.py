"""Rebuild the bundled method comparison from the current analysis pipeline.

Run ``python method_comparison/regenerate.py`` from the repository root.
Only the PNG/SVG figures, metrics.json and README.md in this directory are written.
"""

from copy import deepcopy
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from autoqy_core.config import AnalysisConfig, load_config
from autoqy_core.output import format_value_uncertainty
from autoqy_core.plotting import write_figure
from autoqy_core.runner import run_analysis


CASES = (
    ("340nm_raw_led", "340 nm, raw LED", "Example-1a_340nm_ManualPower_noblcorrLED"),
    ("340nm_baselined_led", "340 nm, pre-baselined LED", "Example-1b_340nm_ManualPower_blcorrLED"),
    ("455nm_first40", "455 nm, first 40 spectra", "Example-2_AB_455nm-100mA"),
    ("455nm_epsilon_error", "455 nm, with ε uncertainty", "Example-3_AB_455nm-EpsilonError"),
    ("395nm_epsilon_error", "395 nm, with ε uncertainty", "Example-4_395nm-EpsilonError"),
)
METHODS = (
    ("concentrations", "Concentrations"),
    ("emission", "Emission (legacy)"),
    ("ode_absorbance", "Full spectrum ODE"),
)


def _pair(values):
    return {"R_to_P": float(values[0]), "P_to_R": float(values[1])}


def _formatted(value, error):
    value_text, error_text = format_value_uncertainty(value, error, two_digit_threshold=2)
    return f"{value_text} ± {error_text}"


def regenerate():
    directory = Path(__file__).resolve().parent
    metrics = []
    guide = ["""# Compare the fitting methods

This guide runs the same five datasets with all three quantum-yield methods.
Use it to see where the methods agree and which inputs need a closer look.
Each PNG has a matching SVG. Start with the [example index](../ExampleData/README.md)
if you want to run one dataset yourself.

## What the methods do

- **Concentrations** (`concentrations`): finds reactant and product concentrations
  from each spectrum, keeps one total concentration, then fits the reaction model.
- **Emission (legacy)** (`emission`): fits the wavelengths covered by the active
  LED band and assumes the sample starts as pure reactant.
- **Full spectrum ODE** (`ode_absorbance`): fits the reaction model directly to
  the full spectral range, including the starting mixture, total concentration,
  and a small baseline correction for each spectrum.

## How to read the tables

**Φ R→P** is the forward quantum yield; **Φ P→R** is the reverse yield.
Values and errors are percentages. Each method uses the same data and error
settings within a dataset. The first three datasets include fit and power
errors. The last two also include the supplied **ε ranges**: AutoQY repeats
the fit with lower, central, and upper reference curves. These combined errors
are sensitivity ranges, not confidence intervals.

**Concentration difference** measures how far a method's fitted concentration
curves are from the Concentrations result, relative to that result's starting
total concentration. It measures agreement with a reference method, not accuracy.

**Absorbance RMSE** summarizes measured minus fitted absorbance over the same
full spectral range. Smaller means closer agreement. For the last two datasets,
parentheses give the RMSE range across the ε calculations. Initial and final
R fractions describe the fitted mixture at the first and last measurement.

Full spectrum ODE can fit baseline changes, so its smaller RMSE alone does
not prove it is the better method. The datasets do not include the individual
absorbance errors needed for a formal statistical ranking.
See [metrics.json](metrics.json) for unrounded results and separate error components.

To rebuild the tables and figures after installing AutoQY, run from the repository folder:

```text
python method_comparison/regenerate.py
```
"""]
    for stem, title, folder in CASES:
        config_path = ROOT / "ExampleData" / folder / "generic_inputs" / "analysis.json"
        source = load_config(config_path)
        epsilon_method = source.values.get("uncertainty", {}).get("epsilon", {}).get("method", "none")
        case = {
            "id": stem, "title": title,
            "configuration": config_path.relative_to(ROOT).as_posix(),
            "epsilon_uncertainty_method": epsilon_method, "methods": [],
        }
        guide.extend([
            f"## {title}\n",
            f"Input: [{folder}](../ExampleData/{folder}/generic_inputs/analysis.json). "
            f"ε uncertainty: `{epsilon_method}`.\n",
            "| Method | Φ R→P (%) | Φ P→R (%) | Concentration difference | Absorbance RMSE (ε range) | Initial R fraction | Final R fraction |\n"
            "|---|---:|---:|---:|---:|---:|---:|",
        ])
        reference = None
        for method, label in METHODS:
            print(f"Running {stem}: {label}", flush=True)
            values = deepcopy(source.values)
            values["fit"]["method"] = method
            for key in ("write_text", "write_figures", "write_json", "write_config", "write_detailed_data"):
                values["outputs"][key] = False
            output = run_analysis(
                AnalysisConfig(values, source.base_directory, source.source), directory
            )
            result, data = output.result, output.data
            trajectory = result.yield_fit.concentrations
            if reference is None:
                reference = trajectory.copy()
            difference = float(np.sqrt(np.mean((trajectory - reference) ** 2)) / reference[0].sum())
            fractions = trajectory[:, 0] / trajectory.sum(axis=1)
            epsilon = np.vstack((result.epsilon_r, result.epsilon_p))
            fitted_absorbance = trajectory @ epsilon * data.path_length_cm
            if result.yield_fit.absorbance_correction is not None:
                fitted_absorbance += result.yield_fit.absorbance_correction
            rmse = float(np.sqrt(np.mean((result.absorbance.T - fitted_absorbance) ** 2)))
            uncertainty = result.epsilon_uncertainty
            rmse_text = f"{rmse:.6f}"
            row = {
                "method": method, "label": label,
                "quantum_yield_percent": _pair(result.yield_fit.values * 100),
                "reported_error_percent": _pair(result.yield_errors * 100),
                "optimizer_power_error_percent": _pair(
                    (result.yield_errors if uncertainty is None else uncertainty.optimizer_power_errors) * 100
                ),
                "epsilon_error_percent": _pair(np.zeros(2) if uncertainty is None else uncertainty.epsilon_errors * 100),
                "epsilon_bound_combination_count": 0 if uncertainty is None else uncertainty.bound_combination_count,
                "concentration_difference_fraction": difference,
                "absorbance_rmse": rmse,
                "initial_reactant_fraction": float(fractions[0]),
                "final_reactant_fraction": float(fractions[-1]),
                "optimizer_success": bool(result.yield_fit.optimizer_success),
                "active_bounds": [bool(item) for item in result.yield_fit.active_bounds],
            }
            if uncertainty is not None:
                lower = uncertainty.absorbance_residual_rmse_minimum
                upper = uncertainty.absorbance_residual_rmse_maximum
                row["epsilon_absorbance_rmse_range"] = [float(lower), float(upper)]
                rmse_text += f" ({lower:.6f}–{upper:.6f})"
            case["methods"].append(row)
            yields = [_formatted(value * 100, error * 100) for value, error in zip(result.yield_fit.values, result.yield_errors)]
            difference_text = "reference" if method == "concentrations" else f"{difference:.2%}"
            guide.append(f"| {label} | {yields[0]} | {yields[1]} | {difference_text} | {rmse_text} | {fractions[0]:.3f} | {fractions[-1]:.4f} |")
            write_figure(directory / f"{stem}__{method}", result, data, values["plots"]["absorbance_residual_percentile"])
            print(f"  Phi = {yields}; absorbance RMSE = {rmse_text}", flush=True)
        guide.append(f"\nFigures: `{stem}__<method>.png` and `.svg`.\n")
        metrics.append(case)
        # Save completed datasets as they finish, so long runs retain their progress.
        (directory / "metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    guide.append("""## What to look for

Compare the yield differences with the reported errors, then inspect the plots.
Similar concentration curves can still give different reverse yields when little
reactant remains. Check the Emission method's pure-reactant assumption against
the recovered starting mixture.

ε ranges can widen the errors and change which methods agree. They do not cover
every choice of reaction model, baseline treatment, or LED processing. Use the
comparison as a check of those choices, alongside the warnings and residual plots.
""")
    (directory / "README.md").write_text("\n".join(guide), encoding="utf-8")


if __name__ == "__main__":
    regenerate()
