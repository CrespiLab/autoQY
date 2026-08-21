# 455 nm ε uncertainty example

This self-contained example propagates wavelength-resolved uncertainty from
both the reactant and product molar-absorptivity spectra using the ε range
method with lower, nominal, and upper curves.

From the repository root, run:

```text
autoqy-core validate ExampleData/Example-3_AB_455nm-EpsilonError/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-3_AB_455nm-EpsilonError/generic_inputs/analysis.json
```

Outputs are written to the input folder's local `results` directory. The
generic reactant CSV contains replicate-derived standard deviations; the
NMR-derived product CSV contains asymmetric wavelength-resolved bounds. The
same data are preserved as legacy TSV/SpectraGryph files in
`crespi_group_inputs`.
