# Tutorial: the 395 nm example

This example shows how to calculate stable and metastable ε, then use them
to calculate quantum yields.

Here, the **reactant** is the stable isomer and the **product** is the metastable
isomer. ε means molar absorptivity: how strongly the compound absorbs light
after allowing for concentration and cuvette path length.

All paths below start in `ExampleData/Example-4_395nm-EpsilonError`.

## Want to try the analysis first?

The `generic_inputs` folder already contains the prepared files.

1. Open the **AutoQY** folder on your Desktop and double-click **AutoQY Analysis GUI**.
2. In **1 · Project**, load `generic_inputs/analysis.json`.
3. Check that **JSON base folder** is the `generic_inputs` folder.
4. Check the output folder before running. The included settings allow overwriting.
5. Click **Run analysis**.

To prepare the spectra yourself, start below. Save your new files in a separate
folder, such as `my_395nm_inputs`, so you do not replace the included files.

## 1. Open Spectral Treatment

In AutoQY Analysis, open **1 · Project** and click **Open Spectral Treatment**.
A second window opens.

## 2. Calculate stable-isomer ε

In Spectral Treatment:

1. Open **1 · Data** and load these three files from
   `crespi_group_inputs/source_data`:

   - `2107362U1_0001A.Abs8`
   - `2107362U1_0001B.Abs8`
   - `2107362U1_0001C.Abs8`
2. In **2 · Range**, enter **250** and **700** nm.
3. Open **Preprocess spectra**. Set the baseline to **600–650 nm** and turn
   **Baseline** on. Press Enter after typing numbers, or click outside the field.
4. Choose **SavGol**, with a **5 nm** window and order **3**.
5. Leave **SVD off**. These are separate sample preparations; keep the
   differences between them so AutoQY can estimate the error.
6. Open **3 · Beer–Lambert** and enter:

| File ending | Concentration (mol/L) | Path length (cm) |
| --- | ---: | ---: |
| `0001A.Abs8` | `6.82e-5` | `1` |
| `0001B.Abs8` | `6.87e-5` | `1` |
| `0001C.Abs8` | `6.80e-5` | `1` |

7. Check that the three ε curves are similar. If one is far from the others,
   check the file order and concentrations before continuing.
8. In **5 · Output**, choose **Processed dataset / calculated ε** under
   **Content to save**. Click **Save CSV**.
9. In Save As, choose your new folder and name the file
   `reactant_absorptivity.csv`.

This file contains the three ε curves, their average, SD, and SEM.
**SD** describes differences between samples. **SEM** is the error estimated
for the average; it is smaller.

Original absorbance curves are hidden by default after treatment.
Use **Show original** if you want to compare them with the treated curves.

## 3. Calculate metastable ε from the PSS mixture

The final mixture contains **23% stable isomer** and **77% metastable isomer**,
measured by NMR. Enter the metastable percentage in Reconstruction.

1. Make sure you saved the ε CSV from step 2.
2. Click **Clear all spectra**. Load
   `crespi_group_inputs/source_data/reactant_pss.dat` into **1 · Data**.
3. Keep the same settings: **250–700 nm**, baseline **600–650 nm**,
   **SavGol 5 nm / order 3**, and **SVD off**.
4. Open **4 · Reconstruction → Reconstruct spectrum**.
5. Select the first spectrum as **Stable spectrum** and the last
   as **PSS spectrum**.
6. Drop your saved `reactant_absorptivity.csv` into
   **Stable ε reference (optional)**.
   You can use `generic_inputs/reactant_absorptivity.csv` instead.
7. Enter **77** for **Metastable (%)**.
8. Leave composition uncertainty at **1** percentage point and choose
   **Replicate SD**. Use zero only if you want to leave the NMR error out.
9. Click **Reconstruct**.

The initial and PSS absorbance curves come from the main window.
The dropped ε file provides the scale and error from the stable-isomer samples.
Both UV measurements must have the same total concentration and path length.

For this mixture, the subtraction is:

`metastable = (PSS − 0.23 × stable) / 0.77`

The plot shows stable, PSS, and metastable ε. The error bands use the stable ε
error and the NMR composition error, with PSS absorbance held fixed.
The reconstruction can produce unequal upper and lower errors.

### Check and save

Small negative values near zero can come from noise. Large negative peaks
need checking. Values below **−500 M⁻¹ cm⁻¹** stop export.

1. In **5 · Output**, choose **Stable + reconstructed ε (pair)**.
2. Click **Save CSV**. Two Save As dialogs open.
3. Save the stable file as `reactant_absorptivity.csv` and the metastable file
   as `product_absorptivity.csv` in your new folder.
   The stable file may replace the one saved in step 2; confirm only if intended.

Negative values are kept by default. To set negative values to zero in the
main exported columns, tick **Convert negative absorbance and ε values to 0**.
The metastable CSV also keeps the original calculated values in a separate column.

For a figure, open **Save figures** and click **Save PNG** or **Save SVG**.
The same content choice is used. **Image export options** lets you hide the
title, legend, or axis labels. **Minimal colors** uses blue, orange, black,
and grey. The reconstruction plot shows the calibrated ε curves.

If you only want reconstructed absorbance, leave the ε reference empty.
Save that result as an absorbance spectrum. Analysis requires the calibrated ε file.

To stop and continue later, click **Save session (.aqys)**.
Drop the saved session into **1 · Data** to reopen it.

## 4. Set up Analysis

Return to **AutoQY Analysis**.

If using the included prepared files, load `generic_inputs/analysis.json`.
If using your new ε files, put copies of these files from `generic_inputs`
in your new folder too:

- `measurement_spectra.csv`
- `led_emission.csv`
- `timestamps.csv`

Then load the example JSON and update its base folder and input paths.
Check the settings below.

### Project and files

Set **JSON base folder** to the folder containing the five input CSVs.
In **2 · Identity and data**, select:

| Input | File |
| --- | --- |
| Measurement spectra | `measurement_spectra.csv` |
| Reactant molar absorptivity | `reactant_absorptivity.csv` |
| Product molar absorptivity | `product_absorptivity.csv` |
| LED emission | `led_emission.csv` |
| Irradiation timestamps | `timestamps.csv` |

Use **Generic CSV (recommended)** for the spectral and timestamp formats.
Use `example_395nm_epsilon_error` as the analysis ID. In **Output stem**, enter
`Results_395nm_EpsilonError`. Saved result file names will start with this text.

### LED settings

In the LED **Processing** section, use:

- wavelength range: **250–800 nm**;
- smoothing window: **12 points**;
- polynomial order: **3**;
- baseline correction: **on**;
- baseline exclusion multiplier: **10**.

These settings apply to the LED file. Prepare absorbance spectra in Spectral Treatment.
The LED window is in points; the Spectral Treatment window was in nm.

### Experiment

In **3 · Experiment**, enter:

| Field | Value |
| --- | ---: |
| Sample volume | **1995 µL** |
| Path length | **1 cm** |
| Power | **1.46 mW** |
| Power error | **0.03 mW** |
| Irradiation wavelength | **395 nm** |
| Thermal P→R rate | **6.3e-5 s⁻¹** |
| Thermal R→P rate | **0 s⁻¹** |

Leave chemical-actinometer mode off for this example.
AutoQY uses the LED spectrum and measured power to calculate photon flux.

### Fit

In **4 · Fit**, choose **Concentrations (legacy pure NNLS)**.

Leave the initial yields at **0.5** and the bounds at **0–1** under
**Expert optimizer settings**. These are starting guesses.

Set **Expected metastable photoproduct at PSS (%)** to **77%** (equivalent to
23% remaining stable reactant). Analysis and Reconstruction both ask for the
metastable percentage. AutoQY compares this optional reference with the fitted PSS.

### Uncertainty

In **5 · Uncertainty**, select **ε ranges** and **Standard deviation** (SD).

AutoQY repeats the calculation using the lower, central, and upper ε values.
For the metastable form, it uses the upper and lower limits saved during
Reconstruction. The results show fit/power error, ε error, and combined error.

## 5. Run and check the results

1. In **6 · Output**, set **Results directory** to `results`.
2. Enable text, figures, JSON, configuration, and detailed-data exports.
   Enable overwrite only if you want to replace an earlier run.
3. Click **Save JSON** to save your settings.
4. Click **Run analysis** to calculate and save the selected results.
5. Read any amber or red warnings. Their explanation is shown beside the flag;
   expand **Technical details** if you need more information.
6. Check **Concentrations**, **Endpoint reconstruction**, and
   **Absorbance residuals**. The calculated spectra should follow the measured
   spectra; large differences or clear patterns in the residuals need checking.

**Compare fit methods** is an extra check. It keeps the ε uncertainty settings
selected above and displays the comparison in the window. Compare **Concentrations (legacy pure
NNLS)** with **Full spectrum ODE**. Each method's errors include its own ε range
calculations. If the yields disagree beyond those ranges, investigate the
concentration curves, reference spectra, and absorbance residuals.

See [the regenerated method comparison](../../method_comparison/README.md) for
the included 395 nm inputs and settings. Newly prepared ε files or different
treatment settings can change the values; use the included files for that check.

## Optional: run from a terminal

From the repository folder:

```text
autoqy-core validate ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json
autoqy-core run ExampleData/Example-4_395nm-EpsilonError/generic_inputs/analysis.json
```
