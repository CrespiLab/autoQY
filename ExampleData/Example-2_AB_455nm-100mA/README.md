# 455 nm: a first concentration analysis

Start here for a simple quantum-yield analysis. It uses the first 40 absorption spectra recorded during 455 nm irradiation.

## Run the prepared example

After installing AutoQY, open **Analysis**, load
`generic_inputs/analysis.json` from this example folder, check the results
folder, and click **Run analysis**. The prepared input files are ready to use.

Or run these commands from the repository folder:

```text
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json
```

`generic_inputs` contains CSV files. `crespi_group_inputs` contains the original
instrument files and a matching configuration. Replace `generic_inputs` with
`crespi_group_inputs` in the commands to run that version. Both give the same results.
Keep the input files beside their JSON when moving this example.

## Settings and expected result

- Method: **Concentrations** (`concentrations`).
- ε uncertainty: **off; reported errors still include fit and power uncertainty**.
- Forward yield, Φ R→P: about **36.8 ± 0.4%**.
- Reverse yield, Φ P→R: about **49.9 ± 0.6%**.

These values use the included files and settings. Changing the treatment or
method can change the result. With ε uncertainty on, the reported error also
includes the effect of the ε ranges. It is not a confidence interval.

## Find the saved files

Results go into `generic_inputs/results` (or `crespi_group_inputs/results`
for the original-file version). File names start with `Results_AB_455nm_first40_ODEMethod_Concentrations`.
The run writes text results, PNG/SVG figures, result and configuration JSON,
and two CSV tables. See [the output guide](../../autoqy_core/README.md#read-the-results).

The included configuration **stops if its output files already exist**.
For another run, choose a new results folder. From a terminal, add
`--output-directory my_results` to the run command, using a folder name you
have not used before.

## What to check

A warning about initial product is expected: the reference spectra give about 3.6% product at the start. Check the starting mixture when interpreting the fit. The Emission method assumes pure reactant and gives a different result; use Compare fit methods to inspect this.

Compare the concentration curves and the **Absorbance residuals** plot as well
as the two yields. The [method comparison](../../method_comparison/README.md)
shows this dataset with all three methods.

For your own data, follow the [beginner guide](../../BEGINNERS-GUIDE.md).

`LegacyResults` contains historical outputs. They are kept for reference and
are not the expected results of the current configuration.
