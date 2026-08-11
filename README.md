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

1. Open [Install-AutoQY.bat](Install-AutoQY.bat) on GitHub.
2. Use GitHub's **Download raw file** button.
3. Put the BAT in the folder where AutoQY should be installed and double-click it.

The installer detects an existing `autoqy-core` environment and asks before
removing it. It never silently deletes the environment or an existing checkout.
If the installer already resides in a valid checkout, that checkout is reused.

The Desktop `AutoQY` folder contains:

- **AutoQY Power GUI** — treat Thorlabs optical-power traces;
- **AutoQY Spectral Smoother** — denoise and export spectral datasets;
- **AutoQY Analyze JSON** — validate and run an `analysis.json` file;
- **AutoQY Terminal** — open a shell with the AutoQY environment activated.

Run a non-destructive installer check with:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\Install-AutoQY.ps1 -CheckOnly
```

Requirements are Anaconda or Miniconda and a current browser. Python, Git, and
the scientific dependencies are installed into the dedicated environment.

## Power GUI

Open **Desktop → AutoQY → AutoQY Power GUI**, or run:

```powershell
autoqy-core power-gui
```

The GUI accepts up to three Thorlabs OPM CSV traces. Drag the boundaries around
the open-beam and cuvette regions, inspect baseline corrections, calculate the
power at the cuvette, and export reproducible JSON settings and results.

Power treatment remains separate from quantum-yield analysis. Transfer the
reported `power_mw` and `power_error_mw` into `analysis.json`; AutoQY does not
modify that file automatically.

## Spectral Smoother GUI

Open **Desktop → AutoQY → AutoQY Spectral Smoother**, or run:

```powershell
autoqy-core smoother-gui
```

The GUI reads both common SpectraGryph `.dat` layouts, wavelength-by-row TSV,
and CSV files. Its processing order is visible and independently configurable:

1. select the wavelength range;
2. optionally zero each spectrum over a baseline interval;
3. apply Savitzky–Golay, Whittaker–Eilers, or no wavelength smoothing;
4. optionally apply SVD rank reduction;
5. inspect the processed spectra and removed residual, then export.

Savitzky–Golay is initially selected with a window expressed in nm. SVD starts
at two components for switching datasets but can be disabled. The GUI reports
the effective smoothing window, RMS changes, and the exact percentage of
squared singular-value weight retained by the selected rank. Original uploaded
data is never overwritten.

Both GUIs run only on the local computer. Uploaded data is held by the browser
session and local Python process; it is not sent to an AutoQY service.

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

### Power processing from JSON

```powershell
autoqy-core power ExampleData/Example-Power/power_analysis.json
```

### Outputs

Depending on configuration, an analysis writes TXT, PNG/SVG figures, results
JSON and input snapshots, concentration/residual TSV, and long-form measured
and reconstructed spectral TSV. Existing files are replaced only when
`outputs.overwrite` is `true`.

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
