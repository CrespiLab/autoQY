# Power processing example

`power_analysis.json` processes the three Thorlabs OPM CSV files independently
from the quantum-yield analysis.

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

Upload the three CSV files from this directory to reproduce the graphical
workflow. The supplied `power_analysis.json` contains the already selected
regions and remains the non-interactive example.

Two uncertainty quantities are reported:

- `repeatability_sd_mw`: conservative root-mean-square variation within the
  corrected on regions;
- `standard_error_mw`: propagated within-region standard errors combined with
  the variation between replicate power estimates.

`uncertainty_output` selects which quantity is returned as `power_error_mw`.
The example retains `repeatability_sd` for compatibility with the historical
reported uncertainty.
