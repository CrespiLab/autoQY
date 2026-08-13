# 455 nm epsilon-error example

This self-contained example propagates wavelength-resolved uncertainty from
both the reactant and product molar-absorptivity spectra using deterministic
lower, nominal, and upper bounds.

From the repository root, run:

```text
autoqy-core validate ExampleData/Example-3_AB_455nm-EpsilonError/analysis.json
autoqy-core run ExampleData/Example-3_AB_455nm-EpsilonError/analysis.json
```

Outputs are written to the local `results` folder. The reactant TSV contains
replicate-derived standard deviations; the NMR-derived product TSV contains
asymmetric wavelength-resolved lower and upper bounds.
