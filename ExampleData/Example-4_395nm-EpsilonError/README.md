# 395 nm with ε uncertainty

Use this example to practise the full workflow: prepare the reactant reference, reconstruct product ε from an NMR mixture composition, and calculate quantum yields. ε means molar absorptivity.

For the complete preparation workflow, follow the [395 nm tutorial](TUTORIAL.md).

## Run the prepared example

After installing AutoQY, open **Analysis**, load
`generic_inputs/analysis.json` from this example folder, check the results
folder, and click **Run analysis**. The prepared input files are ready to use.

Or run these commands from the repository folder:

```text
autoqy-core validate ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json
```

`generic_inputs` contains CSV files. `crespi_group_inputs` contains the original
instrument files and a matching configuration. Replace `generic_inputs` with
`crespi_group_inputs` in the commands to run that version. Both give the same results.
Keep the input files beside their JSON when moving this example.

## Settings and expected result

- Method: **Concentrations** (`concentrations`).
- ε uncertainty: **on, using SD and the saved product ε ranges**.
- Forward yield, Φ R→P: about **24.1 ± 0.7%**.
- Reverse yield, Φ P→R: about **20.3 ± 1.8%**.

These values use the included files and settings. Changing the treatment or
method can change the result. With ε uncertainty on, the reported error also
includes the effect of the ε ranges. It is not a confidence interval.

## Find the saved files

Results go into `generic_inputs/results` (or `crespi_group_inputs/results`
for the original-file version). File names start with `Results_395nm_EpsilonError`.
The run writes text results, PNG/SVG figures, result and configuration JSON,
and two CSV tables. See [the output guide](../../autoqy_core/README.md#read-the-results).

The included configuration **allows replacement of existing results**.
For another run, choose a new results folder. From a terminal, add
`--output-directory my_results` to the run command, using a folder name you
have not used before.

## What to check

The reactant file has negative ε values in a low-signal region. AutoQY reports them and sets them to zero. A warning about the initial product fraction (about 2.7%) is also expected. The separately measured PSS mixture contains 23% reactant and 77% product; enter 77 in the GUI product-percentage field. This value checks the result and does not constrain the fit.

Compare the concentration curves and the **Absorbance residuals** plot as well
as the two yields. The [method comparison](../../method_comparison/README.md)
shows this dataset with all three methods.

For your own data, follow the [beginner guide](../../BEGINNERS-GUIDE.md).
