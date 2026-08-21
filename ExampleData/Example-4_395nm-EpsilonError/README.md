# 395 nm ε uncertainty example

This self-contained 395 nm dataset propagates the
wavelength-resolved errors of both reactant and product molar absorptivity with
the ε range method.

See [TUTORIAL.md](TUTORIAL.md) for the complete GUI workflow, beginning with
the three raw Avantes measurements in `crespi_group_inputs/source_data`.

From the repository root, run:

```text
autoqy-core validate ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json
```

Outputs are written to the local `results` folder. The configuration uses the
regularized-concentrations fit and standard deviation (`sd`) columns from the
two molar-absorptivity CSV files. The original group files and a matching
configuration remain in `crespi_group_inputs`.

The reactant input contains negative mean ε values in its low-signal
region. AutoQY reports this and constrains those values to zero before fitting.
