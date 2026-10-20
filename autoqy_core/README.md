# Using AutoQY Core

AutoQY calculates quantum yields from a series of absorption spectra. You can
use its windows to prepare and analyze data, or run a saved analysis from a terminal.
Your data stay on your computer.

## Start with an example

After installing AutoQY, open **AutoQY Analysis GUI** from the AutoQY Desktop
folder, or run:

```text
autoqy-core analysis-gui
```

In **1 · Project**, load
`ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json`.
Check the results folder, then click **Run analysis**. Paths in this guide
start at the repository folder unless stated otherwise.

For a step-by-step introduction, use the [395 nm tutorial](../ExampleData/Example-4_395nm-EpsilonError/TUTORIAL.md).
For your own measurements, see the [beginner guide](../BEGINNERS-GUIDE.md).
The [example index](../ExampleData/README.md) explains all six examples.

## Choose a window

| Window | What you use it for | Terminal command |
|---|---|---|
| Analysis | Load settings, calculate quantum yields, and check the fit | `autoqy-core analysis-gui` |
| Spectral Treatment | Treat spectra, calculate ε, and reconstruct a pure-isomer spectrum | `autoqy-core smoother-gui` |
| Power | Process power-monitor measurements | `autoqy-core power-gui` |

**ε** (epsilon, or molar absorptivity) describes how strongly a compound absorbs
light after allowing for concentration and cuvette path length. **R** means
reactant; **P** means product. **PSS** is the mixture predicted when continued
irradiation no longer changes the composition.

On Windows, each window normally opens separately in Edge or Chrome. Closing
it stops its local server. If a separate window cannot open, AutoQY uses your
default browser. The version is shown in each window.

## What an analysis needs

For the usual LED-and-power workflow, prepare:

- absorption spectra measured during irradiation;
- reactant and product ε spectra;
- the LED emission spectrum;
- one irradiation time in seconds for each absorption spectrum;
- sample volume, path length, optical power and its error, and any thermal rates.

Save the file locations and settings in `analysis.json`. Relative file paths
start at the folder containing that JSON, not the terminal's current folder.
No changes to Python code or `ExpParams.txt` are needed.

### Save and share an experiment

In Analysis, **Save JSON** saves the editable experiment. **Run analysis** runs
it; running alone does not save that JSON. Both buttons check the inputs first.

Leave **Include all input data in JSON (recommended)** enabled to put copies of
the inputs inside the saved JSON. You can then move it to another computer
with a compatible AutoQY installation and choose a new results folder.
Turn this option off to save links to the original input files instead.
Use Browse to replace an input that was already embedded in a loaded JSON.

The output setting `write_config` separately saves a settings snapshot beside
the results. With `embed_inputs: true`, saved configurations and result JSON
files include the inputs. Analysis can reopen either file.

### Run from a terminal

Run these commands after installing AutoQY and activating its environment:

```text
autoqy-core validate analysis.json
autoqy-core run analysis.json
autoqy-core run analysis.json --output-directory my_results
autoqy-core power power_analysis.json
```

`validate` checks the configuration and input paths; `run` also reads and fits
the data. Use a new output folder for another run, or enable overwrite in the
configuration if you intend to replace existing results.

If the `autoqy-core` command is unavailable in an otherwise installed environment,
use `python -m autoqy_core` followed by the same arguments. The window shortcuts
`autoqy-analysis-gui`, `autoqy-smoother-gui`, and `autoqy-power-gui` are also available.

## Choose and compare fitting methods

Start with **Concentrations (legacy pure NNLS)**, then compare it with
**Full spectrum ODE**. NNLS means non-negative least squares: recovered
concentrations cannot be negative. ODE refers to the reaction model that
describes how concentrations change over time.

| Method in Analysis | JSON value | How it works |
|---|---|---|
| Concentrations (legacy pure NNLS) | `concentrations` | Finds concentrations from each spectrum, keeps one total concentration, then fits the reaction model. |
| Full spectrum ODE | `ode_absorbance` | Fits the reaction model directly to the full spectral range, including total concentration and the starting mixture. It can also fit a small baseline correction. |
| Emission (legacy) | `emission` | Fits absorbance in the active LED band and assumes a pure-reactant starting sample. Use it to reproduce older analyses. |

**Compare fit methods** runs all three using the same inputs and selected
uncertainty settings. It writes no result files. With **ε range** enabled,
each method includes ε uncertainty in its reported yield errors.
The comparison's RMSE values summarize the differences between measured and
fitted data using the central ε spectra. Smaller RMSE means closer agreement,
but extra baseline terms can also reduce it.

See the [method comparison](../method_comparison/README.md) for results and
figures from all five quantum-yield datasets.

### Optional fit controls

The controls below can stay at their defaults for a first run:

| JSON setting | Default | When it applies |
|---|---:|---|
| `absorbance_baseline_order` | `1` | Full spectrum ODE: `-1` disables baseline fitting, `0` fits an offset, `1` fits an offset and slope per spectrum. |
| `robust_loss_scale` | `0.02` | Full spectrum ODE: controls how strongly isolated large residuals are downweighted. |
| `emission_threshold_fraction` | `0.01` | Emission: selects wavelengths where LED intensity exceeds this fraction of its maximum. |
| `initial_quantum_yields` | `0.5` in each direction | Starting guesses for the fit, not measured yields. |
| `quantum_yield_bounds` | `0` to `1` | Allowed yields as fractions; `0.25` means 25%. |

These keys belong inside `fit`. Settings for other methods are ignored.
**Expected metastable photoproduct at PSS (%)** is an optional check against
a separately measured composition; it does not constrain the fit. Enter 77
for a mixture containing 77% product and 23% reactant.

## Include ε uncertainty

Power uncertainty is included even when ε uncertainty is off. To include
uncertainty in the reference spectra, use AutoQY ε CSV or TSV files containing
error columns for both compounds. In Analysis, choose **ε range** and **SD**
or **SEM**. In JSON, add:

```json
"uncertainty": {
  "epsilon": {
    "method": "deterministic_extremes",
    "error_metric": "sd"
  }
}
```

This is a section of a complete configuration, not a standalone JSON file.
SD describes variation between repeated samples. SEM describes uncertainty in
their average. Product ε reconstructed from NMR composition can have unequal
lower and upper ranges; AutoQY uses the saved ranges directly.

AutoQY repeats the selected fit with distinct combinations of lower, central,
and upper ε curves. It reports fit/power error, ε error, and combined error
separately. The combined error is a sensitivity range, not a confidence interval.
Omit this section or set its method to `none` to leave ε uncertainty out.

The [455 nm](../ExampleData/Example-3_AB_455nm-EpsilonError/README.md) and
[395 nm](../ExampleData/Example-4_395nm-EpsilonError/README.md) examples demonstrate this.

## Read the results

Blue is reactant and orange is product. Open circles show concentrations
recovered from the spectra; solid lines show the fitted reaction model.
With ε uncertainty enabled, error bars and shading show **ε ranges**.
The two yields are labelled **Φ R→P** and **Φ P→R**.

**Absorbance residuals** means measured minus fitted absorbance. Look for
patterns as well as large values. A small reported yield error does not rule
out a poor reference spectrum, baseline problem, or unsuitable reaction model.
Green checks mean no automatic flag; amber asks for review; red flags need
investigation. Expand **Technical details** when more detail is needed.

Each analysis example enables these outputs:

| File | What it contains |
|---|---|
| `.txt` | Readable results and experiment settings |
| `.png` and `.svg` | Four-panel result figure |
| `_results.json` | Results for further calculations |
| `_config.json` | Settings snapshot |
| `_traces.csv` | Concentrations, fractions, and residuals over time |
| `_spectra.csv` | Measured and calculated spectra and their residuals |

File names start with `outputs.stem`. `outputs.directory` is relative to the
configuration folder unless an absolute path is supplied.
The **last timestamp composition** is the fitted mixture at the end of the
measurement. **Extrapolated PSS** predicts the mixture after continued constant
irradiation. These need not be the same.

## Prepare spectra in Spectral Treatment

1. Load spectra in **1 · Data**. Check their order in **Loaded spectra**.
2. Select a wavelength range. Apply baseline correction or smoothing only when needed.
3. To calculate ε, enter each sample's concentration in mol/L and path length in cm.
4. To reconstruct a pure metastable spectrum, select stable and mixture/PSS spectra,
   enter the metastable percentage, then click **Reconstruct**.
5. In **5 · Output**, choose **Content to save**, then **Save CSV** or a figure format.

The [beginner guide](../BEGINNERS-GUIDE.md) explains each step. Keep SVD off for
independently prepared ε samples so their differences remain visible.
For sets larger than 60 spectra, the preview shows at most 60; treatment and
data exports still use the full set.

Reconstruction assumes two absorbing forms, a pure stable reference, and the
same total concentration and path length in the stable and mixture measurements.
An absorbance-only reconstruction cannot replace the product ε file for Analysis.
Add a stable ε reference to calibrate the reconstructed spectrum and include
reference error. NMR composition error defaults to 1 percentage point; enter
zero to omit it. No separate PSS measurement error is included.

**Minimal colors** uses blue for the first spectrum, orange for the last, and
grey for intermediate spectra. The full-color palette is designed for color-blind
readers. **Show original**, legend choices, and image export options control
the preview and figures. **Origin-style export** uses boxed axes, larger text,
and a 1950 × 1500 pixel PNG; SVG can be resized without losing sharpness.

Use **Wavelength slice over time** to inspect one wavelength. **Seconds per
timestamp** converts the stored coordinates to seconds: use `1` for seconds,
or `30` for coordinates `0, 1, 2, ...` measured every 30 seconds.
An optional exponential fit reports a lifetime and its error.

New uploads are added to the current set and receive the current treatment.
**Clear all spectra** clears the window, not your source files. **Displayed data**
lets you inspect and copy the plotted values; use CSV export for the complete dataset.
Negative values are preserved unless you select the option to convert them to zero.

Spectral Treatment autosaves the latest session locally. Wait for its save
confirmation before closing. Use **Save session (.aqys)** for a portable copy;
drop that file into **1 · Data** to reopen it.

## Process power measurements

Use **Power** to load up to three Thorlabs OPM traces. Select each measurement,
place the boundaries around its six LED-off/on regions, correct the baseline,
and inspect the result before calculating power at the cuvette.

The result reports repeatability SD and standard error. Choose which one to
use as `power_error_mw`. Transfer `power_mw` and `power_error_mw` to Analysis;
the two windows do not transfer them automatically.
The [power example](../ExampleData/Example-Power/README.md) includes ready-to-run
settings and explains both error choices. For terminal use, plain CSV with a
`power_mw` column is supported through `generic_csv`.

Missing or infinite power readings are excluded and reported. They remain as
gaps, so the selected region positions do not shift.

## File formats and units

Analysis lets you choose a format for each input through `inputs.formats`.
Generic CSV is the easiest starting point. Use one wavelength column and one
value column for a reference spectrum, or several value columns for a time series.
Timestamps must contain one elapsed irradiation time in seconds per spectrum.

| Data | Supported formats |
|---|---|
| Analysis spectra | Generic delimited text, SpectraGryph text, SPECORD WinASPECT `.dat`, Cary `.DSW/.BSW` |
| Spectral Treatment | Those formats plus Avantes `.Abs8`, with automatic detection |
| Timestamps | Generic delimited text or Crespi Lab AHK event logs |
| Power | Generic CSV or original Thorlabs OPM CSV |

Generic text settings can specify delimiter, header, skipped rows, and columns.
With headers, use column names. Otherwise, columns are numbered from zero.
`"value_columns": "remaining"` selects all remaining data columns.
The bundled configurations give complete examples of each input layout.

| Experiment value | Unit |
|---|---|
| `volume_ul` | µL |
| `path_length_cm` | cm |
| `power_mw`, `power_error_mw` | mW |
| `irradiation_wavelength_nm`, processing wavelength range | nm |
| `thermal_back_reaction_s_1`, `thermal_forward_reaction_s_1` | s⁻¹ |
| `photon_flux_mol_s`, `photon_flux_error_mol_s` | mol photons/s |

Thermal back reaction is P→R; thermal forward reaction is R→P and defaults to zero.
Enter a rate, not a half-life: rate = 0.693 / half-life in seconds.

### Use a chemical actinometer

If you measured photon flux directly, select chemical-actinometer mode in
Analysis. In JSON, set `experiment.irradiation_source` to `chemical_actinometer`
and supply `photon_flux_mol_s`, `photon_flux_error_mol_s`, and
`irradiation_wavelength_nm`, alongside the other experiment settings.
This mode treats the light as a single wavelength; it needs no LED file or power input.
For Example 4, 1.46 ± 0.03 mW at 395 nm corresponds to approximately
(4.8 ± 0.1) × 10⁻⁹ mol photons/s.

### Notes for scripted use

Schema version 1 uses input paths. Schema version 2 embeds each input's original
bytes with a filename, base64 data, and SHA-256 integrity check. Older builds
may not read version 2. `self_contained_document(config)` creates a portable
document from a validated `AnalysisConfig`.

GUI commands accept `--no-open` to start the local server without opening a
window. On Windows, shortcut launch logs are in `%LOCALAPPDATA%\AutoQY\logs`.
Spectral Treatment recovery is stored at
`%LOCALAPPDATA%/AutoQY/spectral-treatment/autosave.aqys` on Windows or
`~/.local/share/AutoQY/spectral-treatment/autosave.aqys` elsewhere.
Only the latest recovery session is retained per user; use `.aqys` copies to
keep separate experiments.
