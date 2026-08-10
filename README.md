<p align="center">
  <img src="autoqy_core/assets/autoqy-logo.png" alt="AutoQY" width="760">
</p>

# AutoQY Core

AutoQY Core is the GUI-independent calculation engine for fitting bidirectional
photoisomerization quantum yields from time-resolved absorption spectra.

The analysis and power-calibration workflows are deliberately separate. Each
has a versioned JSON input and can be run from the terminal, Python automation,
or a future graphical interface.

## Installation

### Prerequisites

Install these before downloading AutoQY Core:

- Python 3.12 or newer;
- `pip` (included with current Python and Conda installations);
- Git, if the repository will be downloaded with `git clone`;
- a current web browser, if the power-treatment GUI will be used.

NumPy, SciPy, pandas, and Matplotlib do **not** need to be installed manually.
They are installed by `pip` with AutoQY Core. The optional power-GUI install
also installs Dash and Plotly. Spectragryph and the Thorlabs Optical Power
Monitor software are not runtime requirements; AutoQY only reads their exported
files.

Using a dedicated Python environment is strongly recommended. Do not install
the project into Conda's `base` environment unless that is intentional.

### Option A: Conda (recommended on Windows)

Open **Anaconda PowerShell Prompt** and run:

```powershell
conda create --name autoqy-core python=3.12 pip
conda activate autoqy-core
conda install git

cd C:\Users\YOUR_USERNAME\Documents\GitHub
git clone https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
```

Replace `YOUR_USERNAME` and the parent folder if the repository should be kept
elsewhere. The final `cd AutoQY-Core` is important: the installation command
must be run from the folder containing `pyproject.toml`.

Install the calculation core only:

```powershell
python -m pip install -e .
```

Or install the core **and** the browser-based power GUI:

```powershell
python -m pip install -e ".[power-gui]"
```

The second command includes the complete core, so it is not necessary to run
both commands.

### Option B: Python `venv`

First install Python 3.12+ and Git using the normal installer for the operating
system. Then clone the repository and create an environment inside it.

Windows PowerShell:

```powershell
cd C:\Users\YOUR_USERNAME\Documents\GitHub
git clone https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[power-gui]"
```

macOS or Linux:

```bash
cd ~/Documents/GitHub
git clone https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[power-gui]"
```

### Verify the installation

With the environment activated and the terminal still in `AutoQY-Core`, run:

```powershell
autoqy-core --help
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

If `autoqy-core` is not found, use the equivalent module form:

```powershell
python -m autoqy_core --help
```

## Quantum-yield analysis

```bash
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

Set `fit.method` to `concentrations`, `regularized_concentrations`,
`ode_absorbance`, or the legacy `emission` method in `analysis.json`. New
analyses should normally compare `regularized_concentrations` and
`ode_absorbance`; see `autoqy_core/README.md` for their assumptions and the
known sensitivity limitation of `emission`. Use `--output-directory PATH` to
redirect a run without editing the configuration.

## Power calibration

```bash
autoqy-core power ExampleData/Example-2_AB_455nm-100mA/Power/power_analysis.json
```

The standalone result contains `power_mw` and `power_error_mw`. Those values can
then be entered into `analysis.json`; the analysis runner never modifies that
file or implicitly runs power processing.

### Browser interface

```bash
autoqy-core power-gui
```

This starts a local server at `http://127.0.0.1:8050` and opens it in the
default browser. Upload one to three Thorlabs OPM CSV files, select the six
off/on regions by dragging the twelve boundaries, inspect both baseline
corrections, and calculate the individual and averaged power. `Export JSON`
saves the selected regions and numerical results. Uploaded data stays in the
browser session and the local Python process; this command does not publish a
website or send data to an external service.

Use `autoqy-power-gui --no-open` if the browser should not open automatically.

## Outputs

An analysis can produce:

- a human-readable TXT summary;
- matching PNG and SVG figures;
- a machine-readable results JSON;
- a snapshot of the exact input configuration;
- a TSV of concentrations, fractions, and their residuals;
- a long-form TSV of measured, reconstructed, and kinetically fitted spectra.

Existing outputs are not replaced unless `outputs.overwrite` is `true`.
Quantum-yield values are rounded to the decimal place justified by their
uncertainty; the JSON retains the unrounded numerical values as well.
The summaries report both the fitted composition at the last experimental
timestamp and the extrapolated photostationary state (PSS) under continued
constant irradiation.

## Python API

```python
from autoqy_core import load_config, run_analysis, run_power_analysis

power = run_power_analysis("power_analysis.json")
config = load_config("analysis.json")
analysis = run_analysis(config)
```

See `autoqy_core/README.md` for configuration and file formats. The scientific
method is described in A. Volker, J. D. Steen and S. Crespi, Beilstein J. Org.
Chem. 2024, 20, 1684-1692, DOI: 10.3762/bjoc.20.150.
