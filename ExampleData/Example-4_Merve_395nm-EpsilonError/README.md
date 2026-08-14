# Merve 395 nm epsilon-error example

This self-contained example uses Merve's 395 nm dataset and propagates the
wavelength-resolved errors of both reactant and product molar absorptivity with
the deterministic-extremes method.

From the repository root, run:

```text
autoqy-core validate ExampleData/Example-4_Merve_395nm-EpsilonError/analysis.json
autoqy-core run ExampleData/Example-4_Merve_395nm-EpsilonError/analysis.json
```

Outputs are written to the local `results` folder. The configuration uses the
regularized-concentrations fit and standard deviation (`sd`) columns from the
two molar-absorptivity TSV files.

The reactant input contains negative mean epsilon values in its low-signal
region. AutoQY reports this and constrains those values to zero before fitting.
