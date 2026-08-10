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

`"method": "concentrations"` first recovers the two concentration traces by
non-negative spectral decomposition and fits the kinetic model to those traces.

`"method": "emission"` fits the kinetic model directly to measured absorbance
within the active LED-emission band. `emission_threshold_fraction` is the
fraction of the processed LED maximum used to define that band.

Both methods use the same rate equations, path length, power uncertainty,
thermal back-reaction, bounds, outputs, and plotting pipeline.

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

## Detailed outputs

When `outputs.write_detailed_data` is true, the traces TSV contains measured
and fitted concentrations, fractions, and separate reactant/product residuals.
The spectra TSV contains original absorbance, concentration reconstruction,
kinetic reconstruction, and both residual matrices in long form.

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
