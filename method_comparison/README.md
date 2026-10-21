# Compare the fitting methods

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
describe sensitivity to the tested ε curves.

**Concentration difference** measures how far a method's fitted concentration
curves are from the Concentrations result, relative to that result's starting
total concentration. Use it to assess agreement with that reference method.

**Absorbance RMSE** summarizes measured minus fitted absorbance over the same
full spectral range. Smaller means closer agreement. For the last two datasets,
parentheses give the RMSE range across the ε calculations. Initial and final
R fractions describe the fitted mixture at the first and last measurement.

Full spectrum ODE includes baseline fitting, which can reduce RMSE. Assess
RMSE alongside the residual plots and the assumptions of each method. A formal
statistical ranking would also require individual absorbance uncertainties.
See [metrics.json](metrics.json) for unrounded results and separate error components.

To rebuild the tables and figures after installing AutoQY, run from the repository folder:

```text
python method_comparison/regenerate.py
```

## 340 nm, raw LED

Input: [Example-1a_340nm_ManualPower_noblcorrLED](../ExampleData/Example-1a_340nm_ManualPower_noblcorrLED/generic_inputs/analysis.json). ε uncertainty: `none`.

| Method | Φ R→P (%) | Φ P→R (%) | Concentration difference | Absorbance RMSE (ε range) | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| Concentrations | 19.76 ± 0.12 | 26.2 ± 0.2 | reference | 0.003112 | 1.000 | 0.0386 |
| Emission (legacy) | 17.64 ± 0.10 | 42.9 ± 0.3 | 1.31% | 0.003625 | 1.000 | 0.0378 |
| Full spectrum ODE | 19.36 ± 0.11 | 24.26 ± 0.17 | 0.44% | 0.002131 | 0.981 | 0.0366 |

Figures: `340nm_raw_led__<method>.png` and `.svg`.

## 340 nm, pre-baselined LED

Input: [Example-1b_340nm_ManualPower_blcorrLED](../ExampleData/Example-1b_340nm_ManualPower_blcorrLED/generic_inputs/analysis.json). ε uncertainty: `none`.

| Method | Φ R→P (%) | Φ P→R (%) | Concentration difference | Absorbance RMSE (ε range) | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| Concentrations | 17.68 ± 0.11 | 42.4 ± 0.5 | reference | 0.002573 | 1.000 | 0.0386 |
| Emission (legacy) | 17.2 ± 0.1 | 41.2 ± 0.3 | 0.33% | 0.002660 | 1.000 | 0.0368 |
| Full spectrum ODE | 17.2 ± 0.1 | 38.6 ± 0.3 | 0.48% | 0.001725 | 0.981 | 0.0363 |

Figures: `340nm_baselined_led__<method>.png` and `.svg`.

## 455 nm, first 40 spectra

Input: [Example-2_AB_455nm-100mA](../ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json). ε uncertainty: `none`.

| Method | Φ R→P (%) | Φ P→R (%) | Concentration difference | Absorbance RMSE (ε range) | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| Concentrations | 36.8 ± 0.4 | 49.9 ± 0.6 | reference | 0.009044 | 0.964 | 0.7771 |
| Emission (legacy) | 25.8 ± 0.5 | 51.4 ± 1.1 | 9.16% | 0.049952 | 1.000 | 0.8370 |
| Full spectrum ODE | 35.1 ± 0.5 | 50.9 ± 0.8 | 2.43% | 0.004201 | 0.999 | 0.7884 |

Figures: `455nm_first40__<method>.png` and `.svg`.

## 455 nm, with ε uncertainty

Input: [Example-3_AB_455nm-EpsilonError](../ExampleData/Example-3_AB_455nm-EpsilonError/generic_inputs/analysis.json). ε uncertainty: `deterministic_extremes`.

| Method | Φ R→P (%) | Φ P→R (%) | Concentration difference | Absorbance RMSE (ε range) | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| Concentrations | 40.2 ± 0.7 | 52 ± 3 | reference | 0.006469 (0.004998–0.008596) | 0.962 | 0.7588 |
| Emission (legacy) | 28 ± 2 | 53 ± 2 | 9.37% | 0.046974 (0.040787–0.052551) | 1.000 | 0.8225 |
| Full spectrum ODE | 37.9 ± 1.2 | 52 ± 3 | 1.83% | 0.004509 (0.004350–0.004706) | 0.986 | 0.7716 |

Figures: `455nm_epsilon_error__<method>.png` and `.svg`.

## 395 nm, with ε uncertainty

Input: [Example-4_395nm-EpsilonError](../ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json). ε uncertainty: `deterministic_extremes`.

| Method | Φ R→P (%) | Φ P→R (%) | Concentration difference | Absorbance RMSE (ε range) | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| Concentrations | 24.1 ± 0.7 | 20.3 ± 1.8 | reference | 0.005992 (0.002846–0.013346) | 0.973 | 0.2332 |
| Emission (legacy) | 25.2 ± 1.6 | 21 ± 3 | 0.87% | 0.002560 (0.002192–0.007021) | 1.000 | 0.2302 |
| Full spectrum ODE | 24.6 ± 0.8 | 20.8 ± 1.8 | 1.00% | 0.002029 (0.001975–0.002740) | 0.993 | 0.2336 |

Figures: `395nm_epsilon_error__<method>.png` and `.svg`.

## What to look for

Compare the yield differences with the reported errors, then inspect the plots.
Similar concentration curves can still give different reverse yields when little
reactant remains. Check the Emission method's pure-reactant assumption against
the recovered starting mixture.

ε ranges show how the supplied reference uncertainties affect the yields and
method agreement. Check the reaction model, baseline treatment, and LED
processing alongside the reported ranges, warnings, and residual plots.
