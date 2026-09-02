# AutoQY Core configuration

One `analysis.json` contains every input path, unit-bearing experimental value,
processing choice, fit setting, plot setting, and output instruction required
for a reproducible quantum-yield run. Relative paths are resolved from the JSON
file location.

Power-monitor processing is intentionally separate; see
`ExampleData/Example-Power/generic_inputs/power_analysis.json`.

## Minimal files for one analysis

A run needs the installed `autoqy_core` package and six experiment files:

- `analysis.json`;
- time-series absorbance spectra;
- reactant molar absorptivity;
- product molar absorptivity;
- LED emission spectrum;
- timestamps.

No `ExpParams.txt`, terminal prompts, or Python source edits are required.

## Terminal interface

```bash
python -m autoqy_core validate analysis.json
python -m autoqy_core run analysis.json
python -m autoqy_core power power_analysis.json
python -m autoqy_core analysis-gui
python -m autoqy_core smoother-gui
python -m autoqy_core power-gui
```

The equivalent installed commands are `autoqy-core`, `autoqy-analysis-gui`,
`autoqy-smoother-gui`, and `autoqy-power-gui`.

## Local GUI behavior

Analysis, Spectral Treatment, and Power display the version read from
`pyproject.toml`. On Windows, each GUI opens in a dedicated Edge or Chrome app
window with a temporary isolated browser profile, rather than becoming another
tab in the user's main browser. The default browser is used when app-window
mode is unavailable. All pages and uploaded data remain on the local computer.

The GUI window and Python server monitor each other. Closing the app window
stops its server, and closing the terminal or Python process closes the browser
process tree. A normal browser tab that cannot close itself instead displays a
message telling the user to close the page and restart the GUI. Pass
`--no-open` to run a GUI server without opening a window.

### Analysis GUI

Run `autoqy-core analysis-gui` to build a complete `analysis.json`, load an
existing configuration, or run an analysis. Its six control panels cover the
project location and files, experiment, fit, uncertainty, and output settings.
Native file and folder selectors are available for local paths, and the live
JSON preview shows the configuration being constructed.

Both **Save JSON** and **Run analysis** validate every referenced file and
setting. There is no separate validation button. **Run analysis** does not save
the editable `analysis.json`; **Save JSON** does. The
`outputs.write_config` option independently controls the configuration snapshot
written beside analysis results. Portable JSON can store input paths relative
to the saved JSON file.

After a run, the GUI reports the two quantum yields and fit status, lists the
generated files, and provides interactive tabs for concentrations, fraction
residuals, preprocessing, reference/reconstruction compatibility, and
wavelength-resolved absorbance residuals. Automatic diagnostics are grouped as
green, amber, or red checks. **Compare fit methods** runs regularized
concentrations, full-spectrum ODE, and legacy concentration fitting with the
same nominal inputs; it writes no result files and disables epsilon uncertainty
for the comparison. **Open Spectral Treatment** starts that GUI in its own
dedicated window.

## Fitting methods

`"method": "regularized_concentrations"` is the recommended concentration
route. It fits all spectra together with a single conserved total concentration.
Each timestamp retains an independent reactant fraction, but those fractions
have a soft exponential envelope with a free starting fraction and free plateau.
The exponential is only a regularizer, not the quantum-yield model; the quantum
yields are still obtained from the full photochemical ODE.
`regularization_strength` controls the soft constraint. Larger values follow
the envelope more closely; the default is `1.0`.

`"method": "ode_absorbance"` jointly fits the quantum yields, total
concentration, and initial composition directly to the full measured spectral
range. The concentration trajectory follows the photochemical ODE. An optional
small baseline correction is fitted independently to each spectrum, and a
robust loss reduces the influence of isolated bad wavelengths. This method is
slower but avoids deriving the result from independently clipped concentration
points. `absorbance_baseline_order` accepts `-1` (off), `0` (offset), or `1`
(offset and slope); the default is `1`. `robust_loss_scale` defaults to `0.02`.

`"method": "concentrations"` is the legacy pure-NNLS concentration route. It
recovers every time point independently by non-negative spectral decomposition,
then fits the photochemical ODE to the resulting traces. It is fast and
transparent, but a mismatch between experimental and reference spectra can pin
several consecutive points to exactly 100% reactant or product. It has no
temporal information during spectral decomposition.

`"method": "emission"` is retained to reproduce the older direct-absorbance
calculation. It fits only wavelengths inside the active LED-emission band and
assumes the first spectrum contains pure reactant. When the absorbance area and
spectral shape sampled by the LED change little during irradiation—for example,
when reactant and product absorptivities are similar in that band—the fit is
poorly conditioned and can report large quantum-yield errors or concentrations
that disagree with the full-spectrum methods. `emission_threshold_fraction`
defines the active band as a fraction of the processed LED maximum.

All four methods use the same rate equations, path length, power uncertainty,
two thermal directions, bounds, outputs, and plotting pipeline. The thermal
part of the reactant derivative is
`d[R]/dt = … + k_P→R[P] - k_R→P[R]`; the product derivative is its negative.
`thermal_forward_reaction_s_1` is optional and defaults to zero, so existing
configurations retain their prior behavior. A practical configuration
containing every method-specific control is:

```json
"fit": {
  "method": "regularized_concentrations",
  "regularization_strength": 1.0,
  "absorbance_baseline_order": 1,
  "robust_loss_scale": 0.02,
  "emission_threshold_fraction": 0.01,
  "initial_quantum_yields": {"R_to_P": 0.5, "P_to_R": 0.5},
  "quantum_yield_bounds": {"minimum": 0.0, "maximum": 1.0}
}
```

Controls unused by the selected method are ignored. For a new dataset, compare
`regularized_concentrations` with `ode_absorbance` and inspect both the
concentration and wavelength-resolved residuals. Agreement is stronger evidence
than either result alone.

### Bundled azobenzene comparison

The table reports concentration-trajectory RMSE relative to `concentrations`,
divided by the initial total concentration. It is a regression check for the
bundled data, not a general accuracy estimate.

| Example | `emission` | `regularized_concentrations` | `ode_absorbance` |
|---|---:|---:|---:|
| 340 nm, raw LED | 1.31% | 0.93% | 0.44% |
| 340 nm, pre-baselined LED | 0.33% | 0.76% | 0.48% |
| 455 nm, first 40 spectra | 9.16% | 1.53% | 2.43% |

Both new methods therefore reproduce the established concentration trajectory
within 2.5% normalized RMSE in all three examples. The older `emission` result
also agrees for the 340 nm examples but not for the 455 nm example, illustrating
why agreement should be checked rather than assumed.

## Input formats

Each input has its own entry under `inputs.formats`, allowing mixed formats in
one experiment.

### Generic CSV (recommended)

Use `generic_delimited` with a comma delimiter and a header row. Measurement
CSV files place wavelength in the first column and one spectrum in every
remaining column; single-spectrum files use one value column.

```json
{"type": "generic_delimited", "delimiter": ",", "header": true}
```

### SpectraGryph / Crespi group

```json
{"type": "spectragryph_tsv"}
```

This reads a tab-separated header, uses the first column as wavelength in nm,
and ignores a column named `Wavenumbers [1/cm]`.

### Other generic delimiters

```json
{
  "type": "generic_delimited",
  "delimiter": ",",
  "header": false,
  "skip_rows": 0,
  "wavelength_column": 0,
  "value_columns": [1]
}
```

With a header, columns may be named instead of numbered. For measurement data,
`"value_columns": "remaining"` selects every non-wavelength, non-ignored
column.

### Crespi Lab AHK timestamps

```json
{"type": "ahk_csv"}
```

The core derives cumulative irradiation time from `LEDon` and `LEDoff` events.

### Generic timestamps

```json
{
  "type": "generic_delimited",
  "delimiter": ",",
  "header": true,
  "time_column": "irradiation_seconds"
}
```

The selected column must contain elapsed irradiation time in seconds and have
one value per measured spectrum.

## Units

| Field | Unit |
|---|---|
| `volume_ul` | microlitres |
| `power_mw`, `power_error_mw` | mW |
| `thermal_back_reaction_s_1` | s^-1 |
| `thermal_forward_reaction_s_1` | s^-1 |
| `irradiation_wavelength_nm` | nm |
| `path_length_cm` | cm |
| `wavelength_range_nm` | nm |

Volume is converted internally from microlitres to millilitres.

## Optional molar-absorptivity uncertainty

The established error-less calculation remains the default. To propagate the
wavelength-resolved errors exported by Spectral Treatment, point
`inputs.reactant_absorptivity` and `inputs.product_absorptivity` to their AutoQY
epsilon CSV files (legacy TSV remains supported) and add:

```json
"uncertainty": {
  "epsilon": {
    "method": "deterministic_extremes",
    "error_metric": "sd"
  }
}
```

`sd` is the default for independently prepared samples; `sem` is also
available. The product input may instead be an NMR-derived product epsilon CSV,
whose asymmetric non-negative bounds are used directly. The core evaluates the
distinct low/mean/high reactant-product epsilon-bound combinations through the
complete selected fit and reports optimizer-plus-power, epsilon-only, and
combined conservative error envelopes separately. Omitting this section, or
setting `method` to `none`, preserves the original calculation and input formats.

Future import parsers can target the same CSV contract with normalization and
concentration metadata left empty or zero and wavelength error set to zero; the
uncertainty loader only requires the wavelength, nominal epsilon, and selected
error columns.

## Detailed outputs

When `outputs.write_detailed_data` is true, the traces CSV contains measured
and fitted concentrations, fractions, and separate reactant/product residuals.
The spectra CSV contains original absorbance, concentration reconstruction,
kinetic reconstruction, and both residual matrices in long form.

The TXT and results JSON distinguish two endpoint values:

- **Composition at the last timestamp** is the fitted composition at the final
  experimental measurement time.
- **Extrapolated PSS** is the composition predicted by the fitted model when
  constant irradiation is continued until `dC/dt = 0`. It includes the
  configured thermal rates in both directions.

## Spectral Treatment GUI

Run `autoqy-core smoother-gui` (or `autoqy-smoother-gui`) to load a
SpectraGryph `.dat`, Avantes `.Abs8`, wavelength-by-row TSV, or CSV data. File
types are detected automatically. **Open files from folder** also makes the
source directory the default export location, and a loading indicator remains
visible while large groups of spectra are parsed.

The main workflow selects a wavelength range and optionally applies a baseline
interval, Savitzky–Golay smoothing, and SVD rank reduction. Baseline and
Savitzky–Golay operate on each spectrum independently. SVD mixes columns and is
therefore intended for ordered time-series spectra, not independent replicate
solutions. The suggested rank retains at least 99.5% of squared singular-value
weight, but it remains an operator decision. Uploaded data is never modified.

Processed absorbance can be exported without concentration values. When every
solution concentration and path length is supplied, the same workflow exports
individual molar-absorptivity spectra and their wavelength-resolved statistics.
The save folder and CSV name are editable, missing `.csv` extensions are added,
and existing files require confirmation before replacement. Negative values are
preserved by default; the output option can clamp saved absorbance and epsilon
values to zero without changing the preview or source data.

The plot legend can be hidden when many spectra are loaded. Its labels can also
be edited one per line; these edits affect only the plot and never change CSV
labels or source filenames. **Minimal colors** highlights the initial spectrum
in blue and the final spectrum in orange, with intermediate spectra in grey as
in the Analysis GUI. PNG and SVG save buttons open a Save As dialog. Saved
images contain only the plots by default; enable **Save title + legend** when
those should be included.

The optional NMR-guided PSS workflow loads a dataset whose first spectrum is
reactant and last spectrum is the final PSS. It supports its own baseline and
Savitzky–Golay settings, propagates the entered PSS composition uncertainty,
and exports reactant and reconstructed product epsilon CSV files. The primary
product epsilon is clipped at zero by default, while a raw audit column is
always retained; **Keep negative values in product epsilon** changes the primary
export behavior.

## Power-treatment GUI

Install `.[gui]`, then run `autoqy-core power-gui`. Load up to three
Thorlabs OPM CSV traces and select the active measurement from the dropdown.
Drag the twelve vertical boundaries around the six background/signal regions,
choose the baseline polynomial degree, apply the correction, and inspect the
raw and corrected plots before calculating power at the cuvette.

The result can report repeatability standard deviation or standard error. The
GUI shows both the active-measurement result and the combined result, and its
JSON export includes the selected regions so processing can be repeated from
the terminal. Power treatment remains separate from quantum-yield analysis;
transfer `power_mw` and `power_error_mw` into `analysis.json` or the Analysis
GUI when needed.

For portable scripted input, set `"format": "generic_csv"` in
`power_analysis.json` and provide one `power_mw` column. The legacy
`thorlabs_opm_csv` format remains available for original instrument exports.
