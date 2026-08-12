<p align="center">
  <img src="autoqy_core/assets/autoqy-logo.png" alt="AutoQY" width="760">
</p>

# AutoQY Core

AutoQY analyzes photoisomerization quantum yields, treats optical-power traces,
and prepares time-resolved spectra. Its calculation engine remains callable
from JSON, the terminal, and Python, while the two local browser GUIs provide
the easiest entry points for routine data treatment.

## Install on Windows

The guided installer creates a dedicated Conda environment, installs AutoQY,
and places four shortcuts inside **Desktop → AutoQY**.

1. **[Download Install-AutoQY.bat](https://github.com/CrespiLab/autoQY/releases/download/v2.0.0-beta.1/Install-AutoQY.bat)**.
2. Double-click the downloaded file and follow the instructions (you may need
   to give Windows permission to run it).

The installer detects an existing `autoqy-core` environment and asks before
removing it. It never silently deletes the environment or an existing checkout.
If the installer already resides in a valid checkout, that checkout is reused.
The installer will ask for a preferred installation folder.
The installation will produce a series of useful links in a folder on the Desktop.

The Desktop `AutoQY` folder contains:

- **AutoQY Power GUI** — treats Thorlabs optical-power traces;
- **AutoQY Spectral Treatment** — preprocesses spectra and calculates molar absorptivity;
- **AutoQY Analyze JSON** — validates and runs the `analysis.json` file used by AutoQY;
- **AutoQY Terminal** — opens a shell with the AutoQY environment activated.

While the following instruction is not mandatory, it is possible to check the installer non-destructively:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\Install-AutoQY.ps1 -CheckOnly
```

Requirements are Anaconda or Miniconda and a browser. Python, Git, and
the scientific dependencies are installed into the dedicated environment.

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
**AutoQY Power GUI** does not modify the .json file automatically.

## Spectral Treatment GUI

**Spectral Treatment is the default AutoQY interface for spectral work.** Use it
to import and inspect spectra, select wavelengths, apply baseline correction or
Savitzky–Golay smoothing, prepare absorptivity spectra with wavelength-resolved
errors, and perform NMR-guided PSS subtraction. The terminal functions remain
available for reproducible or automated processing.

Open **Desktop → AutoQY → AutoQY Spectral Treatment**, or run:

```powershell
autoqy-core smoother-gui
```

The main AutoQY installer installs this GUI and its Desktop shortcut together
with the calculation engine, Power GUI, terminal, and JSON runner.

The GUI reads both common SpectraGryph `.dat` layouts, Avantes AvaSoft 8
`.Abs8` files, wavelength-by-row TSV, and CSV files. `.Abs8` uploads are
detected automatically and processed exports are written as TSV. Its workflow is
independently configurable:

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

Both GUIs run only on the local computer. Uploaded data is held by the browser
session and local Python process.

## Quantum-yield analysis

The GUI utilities prepare inputs, while quantum-yield fitting remains a
reproducible JSON workflow:

```powershell
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

Set `fit.method` to one of:

- `concentrations` — independent nonnegative spectral decomposition;
- `regularized_concentrations` — concentration fitting with a kinetic envelope;
- `ode_absorbance` — joint full-spectrum kinetic fitting;
- `emission` — legacy active-LED-band fitting.

New datasets should normally compare `regularized_concentrations` and
`ode_absorbance` and inspect both concentration and wavelength-resolved
residuals. Detailed schemas, assumptions, and file-format controls are in
[autoqy_core/README.md](autoqy_core/README.md).

Alternatively, run **AutoQY Analyze JSON**. In this way, it is possible to drag and 
drop `analysis.json` files directly in the window. The script will check their
validity and will ask to run them. An output folder with different output files 
will be generated in the folder of origin of the `analysis.json` file.

Depending on configuration, an analysis writes TXT, PNG/SVG figures, results
JSON and input snapshots, concentration/residual TSV, and long-form measured
and reconstructed spectral TSV. Existing files are replaced only when
`outputs.overwrite` is `true` in the .json file.


### Power processing from JSON

```powershell
autoqy-core power ExampleData/Example-Power/power_analysis.json
```
This method is experimental and planned for future automation. 


## Manual installation

Python 3.12 or newer is required. A dedicated environment is recommended.

```powershell
git clone --branch feature/core-extraction --single-branch https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python -m pip install -e ".[power-gui]"
```

The optional dependency group installs Dash and Plotly for both browser GUIs.
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

The scientific method is described in A. Volker, J. D. Steen and S. Crespi,
*Beilstein J. Org. Chem.* **2024**, 20, 1684–1692,
[doi:10.3762/bjoc.20.150](https://doi.org/10.3762/bjoc.20.150).
