<p align="center">
  <img src="autoqy_core/assets/autoqy-logo.png" alt="AutoQY" width="760">
</p>

# AutoQY

AutoQY analyzes photoisomerization quantum yields, treats optical-power traces,
and prepares time-resolved spectra. Its calculation engine remains callable
from JSON, the terminal, and Python, while a series of local browser GUIs provide
the easiest entry points for routine analysis and data treatment.

The core scientific method is described in A. Volker, J. D. Steen and S. Crespi,
*Beilstein J. Org. Chem.* **2024**, 20, 1684–1692,
[doi:10.3762/bjoc.20.150](https://doi.org/10.3762/bjoc.20.150).

## Install on Windows

The guided installer creates a dedicated Conda environment, installs AutoQY,
and places five shortcuts inside **Desktop → AutoQY**.

1. **[Download Install-AutoQY.bat](https://github.com/CrespiLab/autoQY/releases/download/v2.1.0/Install-AutoQY.bat)**.
2. Double-click the downloaded file and follow the instructions (you may need
   to give Windows permission to run it).

The installer detects an existing `autoqy-core` environment and reuses it by
default; clean recreation remains optional.
The installation will produce a series of useful links in a folder on the Desktop.
The `v2.1.0` installer installs the stable release from `main`.

The Desktop `AutoQY` folder contains:

- **AutoQY Analysis GUI** — builds, validates, runs, and interactively plots an analysis;
- **AutoQY Power GUI** — treats Thorlabs optical-power traces;
- **AutoQY Spectral Treatment** — preprocesses spectra and calculates molar absorptivity;
- **AutoQY Analyze JSON** — retains the lightweight drag-and-drop JSON runner;
- **AutoQY Terminal** — opens a shell with the AutoQY environment activated.

All GUIs run only on the local computer. Uploaded data is held by the browser
session and local Python process.

While the following instruction is not mandatory, it is possible to check the installer non-destructively:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\Install-AutoQY.ps1 -CheckOnly
```

Requirements are Anaconda or Miniconda and a browser. Python, Git, and
the scientific dependencies are installed into the dedicated environment.

## Analysis GUI

**AutoQY Analysis GUI is the default interface for quantum-yield calculations.**
It follows the scientific organization of the legacy v1.1.1 interface while
building the current versioned `analysis.json` format. Use it to:

- select the measurement, ε, LED, and timestamp files and their formats;
- enter experiment, processing, fitting, uncertainty, and output settings;
- load an existing JSON or save a portable JSON with relative paths;
- validate every referenced file and setting before calculation;
- run the analysis and explore concentrations, fraction residuals, measured
  spectra, and wavelength-resolved absorbance residuals interactively;
- open the **Spectral Treatment** GUI when inputs need baseline correction, smoothing,
  ε averaging, or NMR-guided PSS subtraction.

Open **Desktop → AutoQY → AutoQY Analysis GUI**, or run:

```powershell
autoqy-core analysis-gui
```
The Analysis GUI creates the same reproducible JSON used by the terminal:

```powershell
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

For an end-to-end GUI example with raw Avantes spectra and ε range analysis,
follow the [395 nm tutorial](ExampleData/Example-4_395nm-EpsilonError/TUTORIAL.md).

Set `fit.method` to one of:

- `concentrations` — independent nonnegative spectral decomposition;
- `regularized_concentrations` — concentration fitting with a kinetic envelope;
- `ode_absorbance` — joint full-spectrum kinetic fitting;
- `emission` — legacy active-LED-band fitting.

New datasets should normally compare `regularized_concentrations` and
`ode_absorbance` and inspect both concentration and wavelength-resolved
residuals. Detailed schemas, assumptions, and file-format controls are in
[autoqy_core/README.md](autoqy_core/README.md).

Depending on configuration, an analysis writes TXT, PNG/SVG figures, results
JSON and input snapshots, concentration/residual TSV, and long-form measured
and reconstructed spectral TSV. Existing files are replaced only when
`outputs.overwrite` is `true` in the .json file.

The GUI can output TXT, JSON, TSV, PNG, and SVG files.


## Spectral Treatment GUI

**Spectral Treatment is the default AutoQY interface for spectral work.** Use it
to import and inspect spectra, select wavelengths, apply baseline correction, SVD and/or
Savitzky–Golay smoothing, prepare absorptivity spectra with wavelength-resolved
errors, and perform NMR-guided PSS subtraction. The terminal functions remain
available for reproducible or automated processing.

Use **Open files from folder** to retain the source directory as the default
save location. The final NMR save writes the reactant ε dataset as
`epsilon-spectra-reactant.tsv` and the NMR-derived dataset as
`epsilon-spectra-product.tsv`, and asks before replacing either existing file.
Closing the GUI browser window also closes its local terminal process.

Open **Desktop → AutoQY → AutoQY Spectral Treatment**, or run:

```powershell
autoqy-core smoother-gui
```

The GUI automatically detects both common SpectraGryph `.dat` layouts, Avantes
AvaSoft 8 `.Abs8` files, wavelength-by-row TSV, and CSV files. Processed exports
are written as TSV. Its workflow is independently configurable:

1. select the wavelength range;
2. optionally baseline each spectrum over a selected interval;
3. apply Savitzky–Golay or no wavelength smoothing;
4. optionally apply SVD rank reduction to ordered time-series data;
5. calculate and export wavelength-resolved molar absorptivity and uncertainty;
6. optionally reconstruct a product spectrum using a PSS composition measured by NMR.

Savitzky–Golay and SVD are off by default. When SVD is enabled, its proposed
rank is at least two components when the dataset permits. SVD should not be used
for independent repeat measurements because it mixes their variation.
The GUI reports the effective smoothing window, RMS changes, and the exact percentage 
of squared singular-value weight retained by the selected rank. Original uploaded
data is never overwritten.

## Power GUI

Open **Desktop → AutoQY → AutoQY Power GUI**, or run:

```powershell
autoqy-core power-gui
```

The GUI accepts up to three Thorlabs OPM CSV traces. Drag the boundaries around
the open-beam and cuvette regions, inspect baseline corrections, calculate the
power at the cuvette, and export reproducible JSON settings and results to be used
for automation.

Power treatment remains separate from quantum-yield analysis. Transfer the
reported `power_mw` and `power_error_mw` manually into `analysis.json` if necessary; 
**AutoQY Power GUI** does not automatically modify the .json file.


### Power processing from JSON

```powershell
autoqy-core power ExampleData/Example-Power/power_analysis.json
```
This method is experimental and planned for future automation. 


## Manual installation

Python 3.12 or newer is required. A dedicated environment is recommended.

```powershell
git clone --branch main --single-branch https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python -m pip install -e ".[power-gui]"
```

The optional dependency group installs Dash and Plotly for all three browser GUIs.
Install only the headless calculation engine with `python -m pip install -e .`.

Equivalent `venv` setup:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[power-gui]"
```

Verify with:

```powershell
autoqy-core --help
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

## Python API

```python
from autoqy_core import load_config, run_analysis, run_power_analysis

power = run_power_analysis("power_analysis.json")
analysis = run_analysis(load_config("analysis.json"))
```

## Disclaimer

AutoQY is provided as a research and data-analysis tool without any guarantee 
that its results are correct, complete, or suitable for a particular purpose. 
The software developers and contributors accept no responsibility or liability for errors, 
incorrect results, data loss, inappropriate interpretation, or any consequences arising 
from the use or misuse of the software.

Users are solely responsible for verifying the quality of their input data, 
selecting appropriate analysis settings, assessing the validity of the underlying 
assumptions and models, and independently validating any results produced by AutoQY.

Results generated by AutoQY should not be treated as a substitute for scientific judgment, 
appropriate experimental controls, uncertainty analysis, or independent verification. 
Ultimate responsibility for reported data, conclusions, publications, 
and decisions based on AutoQY output remains with the user.

