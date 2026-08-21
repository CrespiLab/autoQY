# AutoQY example data

Every example has the same two input tracks:

- `generic_inputs` contains plain, headered CSV files and is the recommended
  starting point for users outside the Crespi group;
- `crespi_group_inputs` preserves the original SpectraGryph, AHK, Avantes, or
  Thorlabs files and their matching configuration.

Each input folder contains its own portable `analysis.json` (or
`power_analysis.json`). The paired configurations produce numerically
equivalent results. Check all six pairs from the repository root with:

```powershell
python ExampleData/verify_equivalent_inputs.py
```

## Available examples

| Directory | Irradiation | Fitting method | Purpose |
|---|---:|---|---|
| `Example-1a_340nm_ManualPower_noblcorrLED` | 340 nm | `emission` | Raw LED without baseline correction |
| `Example-1b_340nm_ManualPower_blcorrLED` | 340 nm | `emission` | Externally baseline-corrected LED |
| `Example-2_AB_455nm-100mA` | 455 nm | `concentrations` | Legacy pure-NNLS concentration route |
| `Example-3_AB_455nm-EpsilonError` | 455 nm | `concentrations` | Wavelength-resolved ε uncertainty |
| `Example-4_395nm-EpsilonError` | 395 nm | `regularized_concentrations` | Regularized fit with ε uncertainty |
| `Example-Power` | 455 nm | — | Optical-power processing |

For a first command-line test, use the generic CSV version of the 455 nm
example:

```powershell
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json
```

To reproduce the same calculation from the original group files, replace
`generic_inputs` with `crespi_group_inputs`.

The power example follows the same convention:

```powershell
autoqy-core power ExampleData/Example-Power/generic_inputs/power_analysis.json
```

`LegacyResults` directories contain outputs from earlier AutoQY
implementations and are retained only for historical comparison. `ExpParams.txt`
and raw instrument source files live under `crespi_group_inputs`; current
analyses obtain experimental parameters from their JSON configuration.
