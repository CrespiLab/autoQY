# AutoQY Core configuration

One `analysis.json` contains every input path, unit-bearing experimental value,
processing choice, fit setting, plot setting, and output instruction required
for a reproducible quantum-yield run. Relative paths are resolved from the JSON
file location.

Power-monitor processing is intentionally separate and uses
`Power/power_analysis.json`.

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
python -m autoqy_core power-gui
```

## Fitting methods

`"method": "concentrations"` is the legacy-compatible concentration route. It
recovers every time point independently by non-negative spectral decomposition,
then fits the photochemical ODE to the resulting traces. It is fast and
transparent, but a mismatch between experimental and reference spectra can pin
several consecutive points to exactly 100% reactant or product. It has no
temporal information during spectral decomposition.

`"method": "regularized_concentrations"` fits all spectra together with a
single conserved total concentration. Each timestamp retains an independent
reactant fraction, but those fractions have a soft exponential envelope with a
free starting fraction and free plateau. The exponential is only a
regularizer, not the quantum-yield model; the quantum yields are still obtained
from the full photochemical ODE. `regularization_strength` controls the soft
constraint. Larger values follow the envelope more closely; the default is
`1.0`.

`"method": "ode_absorbance"` jointly fits the quantum yields, total
concentration, and initial composition directly to the full measured spectral
range. The concentration trajectory follows the photochemical ODE. An optional
small baseline correction is fitted independently to each spectrum, and a
robust loss reduces the influence of isolated bad wavelengths. This method is
slower but avoids deriving the result from independently clipped concentration
points. `absorbance_baseline_order` accepts `-1` (off), `0` (offset), or `1`
(offset and slope); the default is `1`. `robust_loss_scale` defaults to `0.02`.

`"method": "emission"` is retained to reproduce the older direct-absorbance
calculation. It fits only wavelengths inside the active LED-emission band and
assumes the first spectrum contains pure reactant. When the absorbance area and
spectral shape sampled by the LED change little during irradiation—for example,
when reactant and product absorptivities are similar in that band—the fit is
poorly conditioned and can report large quantum-yield errors or concentrations
that disagree with the full-spectrum methods. `emission_threshold_fraction`
defines the active band as a fraction of the processed LED maximum.

All four methods use the same rate equations, path length, power uncertainty,
thermal back-reaction, bounds, outputs, and plotting pipeline. A practical
configuration containing every method-specific control is:

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

### Spectragryph

```json
{"type": "spectragryph_tsv"}
```

This reads a tab-separated header, uses the first column as wavelength in nm,
and ignores a column named `Wavenumbers [1/cm]`.

### Generic spectra

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
| `irradiation_wavelength_nm` | nm |
| `path_length_cm` | cm |
| `wavelength_range_nm` | nm |

Volume is converted internally from microlitres to millilitres.

## Optional molar-absorptivity uncertainty

The established error-less calculation remains the default. To propagate the
wavelength-resolved errors exported by Spectral Treatment, point
`inputs.reactant_absorptivity` and `inputs.product_absorptivity` to their AutoQY
epsilon TSV files and add:

```json
"uncertainty": {
  "epsilon": {
    "method": "deterministic_extremes",
    "error_metric": "sd"
  }
}
```

`sd` is the default for independently prepared samples; `sem` is also
available. The product input may instead be an NMR-derived product epsilon TSV,
whose asymmetric non-negative bounds are used directly. The core evaluates the
distinct low/mean/high reactant-product epsilon-bound combinations through the
complete selected fit and reports optimizer-plus-power, epsilon-only, and
combined conservative error envelopes separately. Omitting this section, or
setting `method` to `none`, preserves the original calculation and input formats.

Future import parsers can target the same TSV contract with normalization and
concentration metadata left empty or zero and wavelength error set to zero; the
uncertainty loader only requires the wavelength, nominal epsilon, and selected
error columns.

## Detailed outputs

When `outputs.write_detailed_data` is true, the traces TSV contains measured
and fitted concentrations, fractions, and separate reactant/product residuals.
The spectra TSV contains original absorbance, concentration reconstruction,
kinetic reconstruction, and both residual matrices in long form.

The TXT and results JSON distinguish two endpoint values:

- **Composition at the last timestamp** is the fitted composition at the final
  experimental measurement time.
- **Extrapolated PSS** is the composition predicted by the fitted model when
  constant irradiation is continued until `dC/dt = 0`. It includes the
  configured thermal back-reaction.

## Power-treatment GUI

Install `.[power-gui]`, then run `autoqy-core power-gui`. The local browser
interface retains the historical power-processing sequence: load a Thorlabs
trace, drag twelve boundaries around the six background/signal regions, inspect
the open-beam and cuvette baseline corrections, and calculate the power. Up to
three traces can be processed and averaged. Exported JSON includes the selected
regions so the calibration can be repeated from the terminal.

A future full analysis GUI should edit JSON-compatible values, call
`validate_config`, and then call `run_analysis`. Power processing remains a
separate workflow and should only offer its returned `power_mw` and
`power_error_mw` for transfer into the analysis configuration.

## Spectral smoother GUI

Run `autoqy-core smoother-gui` (or `autoqy-smoother-gui`) to load a
SpectraGryph matrix, TSV, or CSV dataset and inspect its singular values,
spectral components, and measurement weights. The proposed rank is the minimum
that retains 99.5% of squared singular-value weight and is only a starting
point; the operator must check that weak chemistry or baseline drift has not
been discarded. Optional baselining subtracts each spectrum's mean over a
selected wavelength interval before decomposition. The exported file contains
the chosen low-rank reconstruction and leaves the original file unchanged.
Primary wavelength smoothing precedes optional SVD reduction. Savitzky–Golay
is initially selected with a physical window in nm; Whittaker–Eilers and no
spectral smoothing are alternatives. SVD can be disabled independently. The
preview distinguishes the input, spectrally smoothed data, final SVD result,
and total removed residual.
