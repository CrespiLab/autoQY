# 455 nm with ε uncertainty

This example includes uncertainty in both reference spectra. ε (epsilon, or molar absorptivity) describes how strongly a compound absorbs light after allowing for concentration and path length.

## Run the prepared example

After installing AutoQY, open **Analysis**, load
`generic_inputs/analysis.json` from this example folder, check the results
folder, and click **Run analysis**. The prepared input files are ready to use.

Or run these commands from the repository folder:

```text
autoqy-core validate ExampleData/Example-3_AB_455nm-EpsilonError/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-3_AB_455nm-EpsilonError/generic_inputs/analysis.json
```

`generic_inputs` contains CSV files. `crespi_group_inputs` contains the original
instrument files and a matching configuration. Replace `generic_inputs` with
`crespi_group_inputs` in the commands to run that version. Both give the same results.
Keep the input files beside their JSON when moving this example.

## Settings and expected result

- Method: **Concentrations** (`concentrations`).
- ε uncertainty: **on, using SD and the saved product ε ranges**.
- Forward yield, Φ R→P: about **40.2 ± 0.7%**.
- Reverse yield, Φ P→R: about **52 ± 3%**.

These values use the included files and settings. Changing the treatment or
method can change the result. With ε uncertainty on, the reported error also
includes the sensitivity to the tested ε curves.

## Find the saved files

Results go into `generic_inputs/results` (or `crespi_group_inputs/results`
for the original-file version). File names start with `Results_AB_455nm_first40_EpsilonError`.
The run writes text results, PNG/SVG figures, result and configuration JSON,
and two CSV tables. See [the output guide](../../autoqy_core/README.md#read-the-results).

The included configuration **allows replacement of existing results**.
For another run, choose a new results folder. From a terminal, add
`--output-directory my_results` to the run command with a new folder name.

## What to check

The reactant file contains SD values from repeated reference spectra. The product file was reconstructed using an NMR mixture composition and has unequal lower and upper ε ranges. AutoQY reports negative reactant ε values in a low-signal region and sets them to zero. A warning about the initial product fraction (about 3.8%) is also expected. Check the corresponding spectra and starting composition when reviewing these warnings.

Compare the concentration curves and the **Absorbance residuals** plot as well
as the two yields. The [method comparison](../../method_comparison/README.md)
shows this dataset with all three methods.

For your own data, follow the [beginner guide](../../BEGINNERS-GUIDE.md).
