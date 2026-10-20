# Try AutoQY with example data

These examples include the input files and saved settings needed to repeat an
analysis. Install AutoQY first. You can use the Analysis window or a terminal.

## Choose an example

| Example | What you will learn | Method in the saved settings |
|---|---|---|
| [1a · 340 nm, raw LED](Example-1a_340nm_ManualPower_noblcorrLED/README.md) | Reproduce an older analysis with a raw LED spectrum | Emission (legacy) |
| [1b · 340 nm, corrected LED](Example-1b_340nm_ManualPower_blcorrLED/README.md) | Run a similar analysis with an LED spectrum corrected before import | Emission (legacy) |
| [2 · 455 nm](Example-2_AB_455nm-100mA/README.md) | **Start here** for a simple concentration analysis | Concentrations |
| [3 · 455 nm with ε uncertainty](Example-3_AB_455nm-EpsilonError/README.md) | Include errors in the reference spectra | Concentrations |
| [4 · 395 nm with ε uncertainty](Example-4_395nm-EpsilonError/README.md) | Follow the [full preparation tutorial](Example-4_395nm-EpsilonError/TUTORIAL.md), from raw spectra to quantum yields | Concentrations |
| [Power](Example-Power/README.md) | Calculate optical power and its error from three monitor traces | Power processing |

ε (epsilon, or molar absorptivity) describes how strongly a compound absorbs
light after allowing for concentration and path length. A quantum yield counts
molecules changed per absorbed photon. AutoQY reports it in percent.

## Run your first analysis

Open **AutoQY Analysis GUI**, load
`Example-2_AB_455nm-100mA/generic_inputs/analysis.json` from this folder,
check the results folder, and click **Run analysis**.

Or run these commands from the repository folder:

```text
autoqy-core validate ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-2_AB_455nm-100mA/generic_inputs/analysis.json
```

If the installed command is unavailable, use `python -m autoqy_core` followed
by the same arguments. [The core guide](../autoqy_core/README.md) explains the
commands, plots, and saved files.

## Two ways to provide the same data

Every example has two input folders:

- **`generic_inputs`**: CSV files with column headings. Use these first.
- **`crespi_group_inputs`**: the original instrument files, with matching settings.

Each folder contains `analysis.json`, or `power_analysis.json` for Power.
Both versions give the same results. Replace `generic_inputs` with
`crespi_group_inputs` in a command to try the original files.
Keep the referenced input files with their JSON when moving an example;
the bundled JSON files contain file links, not embedded copies of the data.

Analysis outputs go into the selected input folder's `results` subfolder.
Some examples allow overwrite and others stop if results already exist;
each example README states which. To keep an earlier run, choose another
results folder or add `--output-directory my_results` to the run command.

Power is a separate calculation. Its result is written beside its JSON as
`PowerProcessing_results.json`. Follow the [power instructions](Example-Power/README.md)
to transfer its values to Analysis.

## Check your results

Each example README gives expected yields and explains any expected warnings.
The [method comparison](../method_comparison/README.md) compares all three
quantum-yield methods on the five datasets, including ε ranges where supplied.
Its result figures are regenerated from the current code.

`LegacyResults` and the power example's historical TXT file are older outputs
kept for reference. They are not the expected results of the current settings.
Current settings come from JSON; the original `ExpParams.txt` files are not used.

To check that all six pairs of input folders give the same numerical results,
run from the repository folder:

```text
python ExampleData/verify_equivalent_inputs.py
```

For your own measurements, continue with the [beginner guide](../BEGINNERS-GUIDE.md).
