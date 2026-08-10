# AutoQY Core

AutoQY Core is the GUI-independent calculation engine for fitting bidirectional
photoisomerization quantum yields from time-resolved absorption spectra.

The analysis and power-calibration workflows are deliberately separate. Each
has a versioned JSON input and can be run from the terminal, Python automation,
or a future graphical interface.

## Install

Python 3.12 or newer is required.

```bash
pip install -e .
```

To include the local browser interface for power treatment:

```bash
pip install -e ".[power-gui]"
```

## Quantum-yield analysis

```bash
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/analysis.json
```

Set `fit.method` to `concentrations` or `emission` in `analysis.json`. Use
`--output-directory PATH` to redirect a run without editing the configuration.

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
