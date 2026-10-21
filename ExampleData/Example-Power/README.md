# Calculate power from monitor measurements

This example uses three measurements to calculate the optical power at the
cuvette and its error. Use these values as inputs to a quantum-yield analysis.

## Run the prepared example

After installing AutoQY, run this command from the repository folder:

```text
autoqy-core power ExampleData/Example-Power/generic_inputs/power_analysis.json
```

The result is saved as `PowerProcessing_results.json` in `generic_inputs`.
The configuration stops if that file already exists. To keep it and run again,
add `--output` with a new absolute file path, or change the `output` setting
in a copy of the configuration. Relative output paths start at the configuration
folder, including paths passed through `--output`.

`generic_inputs` contains CSV files with one `power_mw` column.
`crespi_group_inputs` contains the original Thorlabs OPM files and equivalent
settings. Replace `generic_inputs` with `crespi_group_inputs` in the command
to try those files. Keep all three measurements beside their JSON.

## Inspect the measurements in the Power window

Open **AutoQY Power GUI**, or run:

```text
autoqy-core power-gui
```

Load the three CSV files from `crespi_group_inputs`. The prepared terminal
settings give the region boundaries. Select each trace in the window and
enter those boundaries before calculating power.

Each trace contains LED off/on/off measurements without the cuvette assembly,
then off/on/off with the jacket, cuvette, and solvent. A cubic curve estimates
the background around each on region. Check the corrected traces before using
the result. The supplied JSON files already contain these region choices.

## Choose the error and use the result

The result gives two error measures:

- **Repeatability SD** (`repeatability_sd_mw`): describes variation in the
  corrected signal. This is the choice saved in this example.
- **Standard error** (`standard_error_mw`): estimates uncertainty in the average,
  including differences between the repeated power measurements.

`uncertainty_output` chooses which one is reported as `power_error_mw`.
Copy `power_mw` and `power_error_mw` to the Experiment fields in Analysis.
Copy these values into the corresponding fields in Analysis.

With the supplied settings, expect about **2.067 mW**, with repeatability SD
**0.0226 mW** and standard error **0.00621 mW**. The terminal rounds the
selected result to **2.07 ± 0.02 mW**. For unrounded values, read the generated JSON. The older
`CorrectlyProcessed_Results_PowerProcessing.txt` is a historical record;
its error calculation may differ from the current one.

For a quantum-yield run using this power measurement, try the
[455 nm example](../Example-2_AB_455nm-100mA/README.md).
