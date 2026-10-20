# 340 nm with the raw LED spectrum

This example calculates quantum yields using the raw LED spectrum. LED baseline correction is off in the supplied settings.

## Run the prepared example

After installing AutoQY, open **Analysis**, load
`generic_inputs/analysis.json` from this example folder, check the results
folder, and click **Run analysis**. The prepared input files are ready to use.

Or run these commands from the repository folder:

```text
autoqy-core validate ExampleData/Example-1a_340nm_ManualPower_noblcorrLED/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-1a_340nm_ManualPower_noblcorrLED/generic_inputs/analysis.json
```

`generic_inputs` contains CSV files. `crespi_group_inputs` contains the original
instrument files and a matching configuration. Replace `generic_inputs` with
`crespi_group_inputs` in the commands to run that version. Both give the same results.
Keep the input files beside their JSON when moving this example.

## Settings and expected result

- Method: **Emission (legacy)** (`emission`).
- ε uncertainty: **off; reported errors still include fit and power uncertainty**.
- Forward yield, Φ R→P: about **17.64 ± 0.10%**.
- Reverse yield, Φ P→R: about **42.9 ± 0.3%**.

These values use the included files and settings. Changing the treatment or
method can change the result. With ε uncertainty on, the reported error also
includes the effect of the ε ranges. It is not a confidence interval.

## Find the saved files

Results go into `generic_inputs/results` (or `crespi_group_inputs/results`
for the original-file version). File names start with `Results_340nm_Emission`.
The run writes text results, PNG/SVG figures, result and configuration JSON,
and two CSV tables. See [the output guide](../../autoqy_core/README.md#read-the-results).

The included configuration **allows replacement of existing results**.
For another run, choose a new results folder. From a terminal, add
`--output-directory my_results` to the run command, using a folder name you
have not used before.

## What to check

Compare it with [Example 1b](../Example-1b_340nm_ManualPower_blcorrLED/README.md), which uses an LED spectrum corrected before import. The two examples also use different treated absorption data, so their difference is not a test of LED correction alone.

Compare the concentration curves and the **Absorbance residuals** plot as well
as the two yields. The [method comparison](../../method_comparison/README.md)
shows this dataset with all three methods.

For your own data, follow the [beginner guide](../../BEGINNERS-GUIDE.md).

`LegacyResults` contains historical outputs. They are kept for reference and
are not the expected results of the current configuration.
