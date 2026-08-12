# AutoQY example data

This directory contains reproducible example analyses for AutoQY Core and a
separate example for processing optical-power measurements.

All paths inside the JSON configuration files are relative to the directory
containing that JSON file. The examples can therefore be copied or moved
without editing absolute paths.

## Available examples

| Directory | Irradiation | Fitting method | Purpose |
|---|---:|---|---|
| `Example-1a_340nm_ManualPower_noblcorrLED` | 340 nm | `emission` | Direct absorbance fitting using a raw LED spectrum without baseline correction |
| `Example-1b_340nm_ManualPower_blcorrLED` | 340 nm | `emission` | Direct absorbance fitting using an externally baseline-corrected LED spectrum |
| `Example-2_AB_455nm-100mA` | 455 nm | `concentrations` | Concentration extraction followed by kinetic quantum-yield fitting |
| `Example-Power` | 455 nm | — | Independent processing of Thorlabs optical-power measurements |

For a first command-line test, start with the 455 nm concentration example.

## Running an analysis

Run these commands from the repository root after installing AutoQY Core.

Validate the configuration:

```powershell
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

Run the analysis:

```powershell
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

The generated files are written to the `CoreResults` directory specified in
the corresponding `analysis.json`.

The 340 nm examples can be run in the same way:

```powershell
autoqy-core run ExampleData/Example-1a_340nm_ManualPower_noblcorrLED/analysis.json
autoqy-core run ExampleData/Example-1b_340nm_ManualPower_blcorrLED/analysis.json
```

## Processing power measurements

The power-processing example is independent of the quantum-yield analyses.

Run it non-interactively with the regions stored in its JSON configuration:

```powershell
autoqy-core power ExampleData/Example-Power/power_analysis.json
```

Alternatively, start the browser-based graphical interface:

```powershell
autoqy-core power-gui
```

If AutoQY Core was installed using the Windows installer and desktop shortcuts
were enabled during installation, the power GUI can also be opened using the
**AutoQY Power GUI** shortcut on the desktop.

The GUI allows the integration regions to be adjusted interactively and
exported as a reproducible `power_analysis.json`.

## Directory conventions

- `analysis.json` is the reproducible definition of a quantum-yield analysis.
- `power_analysis.json` is the reproducible definition of a power-processing
  analysis.
- `CoreResults` contains newly generated results.
- `LegacyResults` contains outputs from earlier AutoQY implementations. These
  files are retained for comparison and documentation but are not exact
  numerical regression targets.
- `ExpParams.txt` files are retained for historical reference. Current AutoQY
  Core analyses obtain the required experimental parameters from
  `analysis.json`.

Some directories contain additional processed or smoothed spectra retained
from the original workflow. Only files explicitly referenced by the JSON
configuration are used during an analysis.
