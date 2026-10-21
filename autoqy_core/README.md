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
reactant; **P** means product. These labels follow the chosen reaction direction.
For a back reaction, R can be metastable and P stable. **PSS** is the mixture predicted when continued
irradiation no longer changes the composition.

On Windows, each window normally opens separately in Edge or Chrome. Closing
it stops its local server. If a separate window cannot open, AutoQY uses your
default browser. The version is shown in each window.

## What an analysis needs

**Preprocessing** previews each measurement, ε reference, and LED as soon as
it can be read. Messages identify inputs that need attention while the
available curves remain visible. Run analysis requires a valid set of inputs.

For the usual LED-and-power workflow, prepare:

- absorption spectra measured during irradiation;
- reactant and product ε spectra;
- the LED emission spectrum;
- one irradiation time in seconds for each absorption spectrum;
- sample volume, path length, optical power and its error, and any thermal rates.

Save the file locations and settings in `analysis.json`. Relative file paths
start at the folder containing that JSON. The JSON contains the experiment settings.

### Save and share an experiment

In Analysis, use **Save JSON** to save the editable experiment settings and
**Run analysis** to calculate results. Both buttons check the inputs first.

Leave **Include all input data in JSON (recommended)** enabled to put copies of
the inputs inside the saved JSON. You can then move it to another computer
with a compatible AutoQY installation and choose a new results folder.
Turn this option off to save links to the original input files instead.
Use Browse to replace an input that was already embedded in a loaded JSON.

To recover the stored files, load the JSON and click **Restore input files as
CSV…** in **1 · Project**. Choose a folder for the CSV copies. All data columns
are kept, including ε errors and timestamp events. Instrument formats are
converted to tables. Files are named by their role (for example,
`measurement_spectra.csv`); a number is added if that name already exists.
Restoration uses the data embedded in the loaded JSON. The button becomes
available when those data are present.

The output setting `write_config` separately saves a settings snapshot beside
the results. With `embed_inputs: true`, saved configurations and result JSON
files include the inputs. Analysis can reopen either file.

### Monte Carlo uncertainty

In **5 · Uncertainty**, choose **Monte Carlo** to propagate
the supplied reference ε, NMR composition, power or photon-flux, and kinetic-fit
uncertainties together. PSS reconstruction absorbance stays fixed.
**ε ranges** remains available, and loaded JSON files keep their saved method.

Use ε files containing complete replicate curves. A reconstructed isomer also
needs its matching calibration. Supply the NMR composition uncertainty when
available. If an older export omits it, Monte Carlo assumes ±1 percentage point
(1 SD) and records a visible warning in the GUI, text and JSON results. New Spectral
Treatment exports retain the reference curve and composition used for
reconstruction. AutoQY checks older exports using their saved calibration scale
and recoverable composition; the result records the checks available for each file.

Each simulation draws a shared reference once, reuses it in the reconstruction,
and refits the kinetics. This preserves the calibration relationship between
the two isomers in either reaction direction.

The GUI and figures show the Monte Carlo mean ± standard uncertainty (1 SD).
The error uses one significant digit, or two when it starts with 1 or 2; the
mean uses the same decimal place. Detailed JSON retains the nominal fit,
full-precision values, and asymmetric coverage interval. The default maximum is
2048 simulations. Any whole-number limit from 32 to 100000 is allowed.
Powers of two are optional. After at least 128 successful fits, AutoQY checks
the numerical precision of both mean yields and their SDs every 64 fits. It stops early when
twice their estimated simulation error is below half a unit in the last
reported SD digit. Batch estimates and moment-based estimates are combined;
the larger error estimate is used. This applies the numerical-tolerance
principle in [JCGM 101:2008, section 7.9](https://www.bipm.org/documents/20126/2071204/JCGM_101_2008_E.pdf)
to the reported mean and SD. Parallel workers finish the current block
before the check; a passed check prevents another block from starting.

The interval endpoints have their own precision check in **Technical details**.
An approximate interval can accompany a mean and SD that have passed their check.
If the simulation limit is reached first, the results state that the mean/SD
precision check has not passed and show the remaining simulation error in
**Technical details**. Failed fits and yield-limit hits remain visible.
The kinetic solver reuses fixed wavelength weights and baseline projections
to shorten repeated fits. A saved random seed makes the calculation reproducible.
Detailed output includes `_uncertainty.csv` with the successful simulated yields.

For a longer run, open **Monte Carlo settings** and select **Use multiple CPUs**.
Parallel processing is off by default. AutoQY detects the available CPUs and
suggests up to four workers, using roughly half the available CPUs initially.
The worker limit leaves at least one CPU free for Windows and other applications.
Each worker fits one simulation at a time and uses one numerical-library thread.
More workers use more memory; short runs can take longer because worker startup
takes time. The counter shows the selected worker count.

Input draws keep the same random sequence in either mode. Fits are saved in
simulation order, and precision checks use the same completed prefixes. You can
change worker counts and continue compatible saved runs. Small differences in
floating-point results can occur between numerical-library thread settings.

A live counter shows how many simulations have finished and how many were reused.
Each completed simulation is saved automatically in `.autoqy-monte-carlo`
inside the output folder. To continue a run, keep the same output folder,
inputs, fit settings, and random seed, raise **Maximum simulations**, and
press **Run analysis**. For example, increasing 256 to 4096 reuses the
first 256 simulations. You can close and reopen the GUI between runs.
An interrupted run continues from its last completed simulation.

Changes to the input data, uncertainty inputs, fit method, or seed start
a separate saved run. Coverage and precision settings can reuse the same
draws. Changes to sampling or fitting code start a separate saved run. Updates to
messages or stopping checks keep compatible saved simulations.
To update existing result files, enable **Allow replacing existing results** or
choose a new **Output stem** while keeping the output folder.


Older reconstructed ε exports may lack the saved NMR composition error.
AutoQY continues with an assumed ±1 percentage point (1 SD) and shows a warning.
To replace this assumption, reopen the stable
and PSS absorbance spectra in **Spectral Treatment → 4 · Reconstruction**,
add the stable ε reference, enter **Composition error (± percentage points)**,
and reconstruct and save the ε CSV again. Load the new CSV in Analysis.
An explicitly saved composition error, including zero, is used as supplied.
Other missing calibration or reconstruction details still require a complete export.
**ε ranges** remains available for files that contain the lower and upper ε curves.

The interval endpoint tolerance defaults to 0.005 in yield fraction (0.5
percentage points). This separate check uses batch estimates and binomial
order-statistic brackets, requiring twice the estimated endpoint sampling
error to fall below the tolerance. It reports the precision of the coverage
interval; automatic stopping follows the mean and SD precision check.

The [NIST Uncertainty Machine manual](https://uncertainty.nist.gov/NISTUncertaintyMachine-UserManual.pdf)
describes propagation through a joint input distribution. AutoQY uses this model:

- Reference ε has nonnegative truncated Gaussian marginals. A Gaussian copula
  estimates wavelength correlations from complete replicate deviations. SD or
  SEM parameterizes the untruncated measurement model.
- The NMR fraction has a Gaussian distribution restricted to 0–1. Its supplied
  error represents one standard uncertainty. Each draw uses nonnegative
  constrained subtraction to reconstruct the other isomer.
- Kinetic-fit uncertainty uses a bootstrap of centred whole-spectrum residuals
  across timestamps, assuming those residual spectra are exchangeable.
- Power or photon-flux errors represent one standard uncertainty. These errors,
  NMR composition errors, reference repeatability, and kinetic residuals are
  assumed independent, with the shared calibration preserved explicitly.

The bootstrap supplies the fit contribution. PSS reconstruction absorbance,
concentrations, and path lengths stay at their supplied values. These scope and
sampling assumptions are saved with the results. Wavelength correlations are
approximate with a small replicate set. The plotted intervals apply separately
at each timestamp.

Add these settings to the analysis configuration:

```json
"uncertainty": {
  "epsilon": {
    "method": "monte_carlo",
    "error_metric": "sd",
    "samples": 2048,
    "workers": 1,
    "seed": 0,
    "coverage_probability": 0.95,
    "endpoint_tolerance": 0.005
  }
}
```

For ε ranges, set `"method": "deterministic_extremes"`.

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

New analyses start with **Full spectrum ODE**, which fits the reaction model
directly to the spectra. Loaded JSON files keep their saved method. Compare
with **Concentrations (legacy pure NNLS)** to check agreement. NNLS means
non-negative least squares: recovered concentrations cannot be negative.
ODE refers to the reaction model describing how concentrations change over time.

| Method in Analysis | JSON value | How it works |
|---|---|---|
| Concentrations (legacy pure NNLS) | `concentrations` | Finds concentrations from each spectrum, keeps one total concentration, then fits the reaction model. |
| Full spectrum ODE | `ode_absorbance` | Fits the reaction model directly to the full spectral range, including total concentration and the starting mixture. It can also fit a small baseline correction. |
| Emission (legacy) | `emission` | Fits absorbance in the active LED band and assumes a pure-reactant starting sample. Use it to reproduce older analyses. |

**Compare fit methods** runs all three using the same inputs and selected
uncertainty settings and displays the comparison in the window. With **ε ranges** enabled,
each method includes ε uncertainty in its reported yield errors.
Diagnostics also check the combined error and every tested ε combination.
They flag yield limits, unstable fits, and failed fits even when the central
fit looks precise. A yield limit can truncate the reported sensitivity range.
The comparison's RMSE values summarize the differences between measured and
fitted data using the central ε spectra. Smaller RMSE means closer agreement,
and depends on the baseline terms fitted by each method.

See the [method comparison](../method_comparison/README.md) for results and
figures from all five quantum-yield datasets.

### Optional fit controls

The controls below can stay at their defaults for a first run:

| JSON setting | Default | When it applies |
|---|---:|---|
| `absorbance_baseline_order` | `1` | Full spectrum ODE: `-1` disables baseline fitting, `0` fits an offset, `1` fits an offset and slope per spectrum. |
| `robust_loss_scale` | `0.02` | Full spectrum ODE: controls how strongly isolated large residuals are downweighted. |
| `emission_threshold_fraction` | `0.01` | Emission: selects wavelengths where LED intensity exceeds this fraction of its maximum. |
| `initial_quantum_yields` | `0.5` in each direction | Initial guesses supplied to the optimizer. |
| `quantum_yield_bounds` | `0` to `1` | Allowed yields as fractions; `0.25` means 25%. |

These keys belong inside `fit`. Settings for other methods are ignored.
**Expected metastable photoproduct at PSS (%)** is an optional check against
a separately measured composition. AutoQY compares it with the fitted PSS. Enter 77
for a mixture containing 77% product and 23% reactant. The check uses the
species assigned to P, including in a back-reaction analysis.

## Include ε uncertainty

Power uncertainty is included even when ε uncertainty is off. To include
uncertainty in the reference spectra, use AutoQY ε CSV or TSV files containing
error columns for both compounds. In Analysis, choose **ε ranges** and **SD**
or **SEM**. In JSON, add:

```json
"uncertainty": {
  "epsilon": {
    "method": "deterministic_extremes",
    "error_metric": "sd"
  }
}
```

Add this section to the complete analysis configuration.
SD describes variation between repeated samples. SEM describes uncertainty in
their average. Product ε reconstructed from NMR composition can have unequal
lower and upper ranges; AutoQY uses the saved ranges directly.

AutoQY repeats the selected fit with distinct combinations of lower, central,
and upper ε curves. It reports fit/power error, ε error, and combined error
separately. The combined error describes sensitivity to the tested ε curves.
Use Monte Carlo to calculate a coverage interval under its stated sampling model.
Omit this section or set its method to `none` to leave ε uncertainty out.

The [455 nm](../ExampleData/Example-3_AB_455nm-EpsilonError/README.md) and
[395 nm](../ExampleData/Example-4_395nm-EpsilonError/README.md) examples demonstrate this.

## Read the results

Blue is reactant and orange is product. Open circles show concentrations
recovered from the spectra; solid lines show the fitted reaction model.
**ε ranges** shows variation across the tested ε curves. **Monte Carlo**
reports mean ± standard uncertainty. In the concentration plot, shaded bands
show pointwise fit intervals and error bars show the intervals for recovered
concentrations (95% coverage by default). Detailed data includes these intervals.
The two yields are labelled **Φ R→P** and **Φ P→R**.
The shaded band around the fraction residual shows a pointwise Monte Carlo
coverage interval (95% by default), or the variation across the tested ε curves
when using **ε ranges**.

**Absorbance residuals** means measured minus fitted absorbance. Look for
patterns as well as large values. Assess the reported yield errors together
with the reference spectra, baseline, and suitability of the reaction model.
Green checks show passed diagnostics; amber asks for review; red flags need
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
irradiation. Their difference indicates how far the measurement ended from the predicted PSS.

## Prepare spectra in Spectral Treatment

1. Load spectra in **1 · Data**. Check their order in **Loaded spectra**.
2. Select a wavelength range. Apply baseline correction or smoothing only when needed.
3. To calculate ε, enter each sample's concentration in mol/L and path length in cm.
4. To reconstruct a pure metastable spectrum, select stable and mixture/PSS spectra,
   enter the metastable percentage, then click **Reconstruct**.
5. In **5 · Output**, choose **Content to save**, then **Save CSV** or a figure format.

The **Spectral data** information icon explains loading; the main Output
information icon explains saving. Open and Save As dialogs share the latest
folder you choose, including spectrum/session loading, ε references, CSV,
sessions, and PNG/SVG figures. Clicking an upload box uses the same local file
chooser as **Choose files**. Dragging files onto it still works.

The [beginner guide](../BEGINNERS-GUIDE.md) explains each step. Keep SVD off for
independently prepared ε samples so their differences remain visible.
For sets larger than 60 spectra, the preview shows at most 60; treatment and
data exports still use the full set.

Reconstruction assumes two absorbing forms, a pure stable reference, and the
same total concentration and path length in the stable and mixture measurements.
Analysis uses ε spectra. Add a stable ε reference to calibrate a reconstructed
absorbance spectrum and include
reference error. NMR composition error defaults to 1 percentage point; enter
zero to omit it. PSS reconstruction absorbance stays fixed.

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
**Clear all spectra** clears the loaded set; source files remain on disk. **Displayed data**
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
use as `power_error_mw`. Copy `power_mw` and `power_error_mw` into Analysis.
The [power example](../ExampleData/Example-Power/README.md) includes ready-to-run
settings and explains both error choices. For terminal use, plain CSV with a
`power_mw` column is supported through `generic_csv`.

Missing or infinite power readings are excluded and reported. They remain as
gaps to preserve the positions of the selected regions.

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
For a first-order reaction, lifetime τ = 1 / rate, in seconds.
Analysis shows this lifetime for each positive thermal rate.

### Use a chemical actinometer

If you measured photon flux directly, select chemical-actinometer mode in
Analysis. In JSON, set `experiment.irradiation_source` to `chemical_actinometer`
and supply `photon_flux_mol_s`, `photon_flux_error_mol_s`, and
`irradiation_wavelength_nm`, alongside the other experiment settings.
This mode uses the supplied photon flux at a single irradiation wavelength.
For Example 4, 1.46 ± 0.03 mW at 395 nm corresponds to approximately
(4.8 ± 0.1) × 10⁻⁹ mol photons/s.

### Notes for scripted use

Schema version 1 uses input paths. Schema version 2 embeds each input's original
bytes with a filename, base64 data, and SHA-256 integrity check. Version 2
requires an AutoQY build with embedded-input support. `self_contained_document(config)` creates a portable
document from a validated `AnalysisConfig`.

GUI commands accept `--no-open` to start the local server without opening a
window. On Windows, shortcut launch logs are in `%LOCALAPPDATA%\AutoQY\logs`.
Spectral Treatment recovery is stored at
`%LOCALAPPDATA%/AutoQY/spectral-treatment/autosave.aqys` on Windows or
`~/.local/share/AutoQY/spectral-treatment/autosave.aqys` elsewhere.
Only the latest recovery session is retained per user; use `.aqys` copies to
keep separate experiments.
