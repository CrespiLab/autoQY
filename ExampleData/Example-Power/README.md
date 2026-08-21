# Power processing example

`generic_inputs/power_analysis.json` processes plain CSV files containing a
single `power_mw` column. `crespi_group_inputs/power_analysis.json` processes
the three original Thorlabs OPM CSV files. Both produce equivalent results.

Each trace records six explicit half-open index regions: LED off/on/off without
the cuvette assembly, followed by LED off/on/off with the jacket, cuvette, and
solvent. A cubic baseline is fitted to the two off regions surrounding each on
region. The browser GUI retains draggable lines but can export their final
indices, making the processing reproducible from `power_analysis.json`.

From the repository root, install the optional GUI and start it with:

```bash
pip install -e ".[power-gui]"
autoqy-core power-gui
```

Upload the three CSV files from `crespi_group_inputs` to reproduce the graphical
workflow. Each supplied `power_analysis.json` contains the already selected
regions and remains the non-interactive example.

Two uncertainty quantities are reported:

- `repeatability_sd_mw`: conservative root-mean-square variation within the
  corrected on regions;
- `standard_error_mw`: propagated within-region standard errors combined with
  the variation between replicate power estimates.

`uncertainty_output` selects which quantity is returned as `power_error_mw`.
The example retains `repeatability_sd` for compatibility with the historical
reported uncertainty.
