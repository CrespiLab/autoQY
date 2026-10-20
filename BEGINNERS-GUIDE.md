# AutoQY beginner guide

Use this guide for your own data. To practise with the included files, follow the
[395 nm tutorial](ExampleData/Example-4_395nm-EpsilonError/TUTORIAL.md).

AutoQY has two windows:

- **Spectral Treatment**: clean spectra, calculate ε, reconstruct a spectrum, and save data or figures.
- **Analysis**: calculate quantum yields from your measurements.

Click an **ⓘ** button for help with the control beside it. Press Escape to close it.

## A few terms used in the buttons

- **Absorbance**: the value measured by your UV–Vis instrument.
- **ε** (epsilon, or molar absorptivity): absorbance divided by concentration and cuvette path length. This lets you compare how strongly compounds absorb light.
- **Reactant / stable isomer (R)**: the compound you start with.
- **Product / metastable isomer (P)**: the form made by irradiation.
- **PSS**: the mixture reached when continued irradiation no longer changes its composition.
- **Quantum yield (Φ)**: how many molecules change per absorbed photon (a particle of light). A value of 0.25 means 25%.
- **SD**: how much repeated measurements differ from each other.
- **SEM**: the error estimated for their average. It is smaller than SD; do not choose it just to get a smaller error.

## 1. Open AutoQY

After installing, open the **AutoQY** folder on your Desktop.
Double-click **AutoQY Analysis GUI**. In **1 · Project**, click
**Open Spectral Treatment**. A second window opens.

## 2. Load and check your spectra

In Spectral Treatment, open **1 · Data**. Drop your files into the upload box,
click the box to choose files, or use **Choose files**.

AutoQY reads SpectraGryph text and SPECORD binary `.dat`, Cary `.DSW/.BSW`,
Avantes `.Abs8`, CSV, and TSV files.

Files loaded together are sorted by name: `sample2` comes before `sample10`.
Loading more files adds them to the current set. **Clear all spectra** removes
the current set from the window; it does not delete your files.

Open **Loaded spectra**:

- For a time series, put the spectra in measurement order.
- For separate samples, check which spectrum belongs to each concentration.
- Use the move or remove buttons to fix mistakes.

With more than 60 spectra, the plot shows only 60 spread across the series,
including the first and last. All spectra are still used for treatment, time
slices, and CSV export. Only the first and last legend controls are shown.
Sort a large series before loading it.

## 3. Choose a range and clean the spectra

Open **2 · Range**. Enter the wavelength limits you want to keep.
For example, enter `250` and `700` for 250–700 nm. Press Enter after typing a
number, or click outside the field.

Open **Preprocess spectra**. You do not have to turn everything on.

### Baseline

Use this if the spectrum sits above or below zero where there should be no absorption.

1. Choose a wavelength interval with no real absorption.
2. Enter its start and end.
3. Turn **Baseline** on.

For some samples, 600–650 nm works. Check your own spectrum first.
Do not use an interval containing a real absorption band: you would subtract
part of the signal.

### Smoothing

Use **SavGol** to reduce small, jagged changes caused by noise.
Try a **5 nm** window and polynomial order **3** as a starting point.

Compare smoothing on and off. If a peak moves or changes shape, use a smaller
window or turn smoothing off. AutoQY shows how many detector points the window uses;
the polynomial order must be smaller than that number.

### SVD

Leave **SVD off** when comparing separately prepared samples to calculate ε.
SVD can hide differences between those samples.

For a time series, SVD can reduce noise. If you try it, compare the result with
SVD off. The suggested component count is not a guarantee that the result is right.

Original spectra are hidden by default in the treated-absorbance plot.
Use **Show original** to compare them with the treated curves. This also affects
the saved absorbance figure. It does not change your source files.

## 4. Calculate ε

You need the concentration and cuvette path length for each sample.

1. Load the spectra of the separately prepared samples.
2. Set their range and treatment as above. Leave SVD off.
3. Open **3 · Beer–Lambert**.
4. Enter each concentration in **mol/L** and each path length in **cm**.

For example, `6.82e-5` means 0.0000682 mol/L. A standard 1 cm cuvette needs `1`
in the path-length field.

AutoQY calculates ε = absorbance / (concentration × path length).
With several samples, it also shows the average, SD, and SEM.

Check the individual ε curves before saving. If one is very different, check
its concentration, dilution, cuvette, and baseline. Do not hide the difference
by adding more smoothing.

Concentration fields are not shown for sets larger than 60 spectra.
Calculate ε from a smaller, separate set of known-concentration samples.

## 5. Reconstruct the pure metastable spectrum

If you can measure the pure metastable form, calculate its ε as in step 4.
You do not need reconstruction.

Use this when you have a pure stable-isomer spectrum and a PSS mixture, but
cannot measure the pure metastable form.

You need:

- the initial stable and PSS absorbance spectra;
- the percentage of metastable isomer in the PSS, measured separately, for example by NMR.

Both UV–Vis measurements must use the same total concentration and path length.
This calculation assumes only two absorbing forms. Do not scale each peak to
the same height before using it.

1. Load the UV–Vis series in **1 · Data**.
2. Set the range, baseline, and smoothing in the main controls.
3. Open **4 · Reconstruction → Reconstruct spectrum**.
4. Select **Stable spectrum** and **PSS spectrum**.
   The first and last spectra are selected by default; check them.
5. Enter **Metastable (%)**.
   If NMR gives 23% stable and 77% metastable, enter **77**, not 23.
6. Enter the composition error in percentage points.
   The default `1` means ±1 percentage point, not 1% of the entered value.
   Use `0` to leave this error out.
7. Click **Reconstruct**.

Without an ε reference, the result is an absorbance spectrum.
Use **Show reconstructed isomer** to show or hide it.

### If you also need metastable ε and its error

Before clicking **Reconstruct**, drop the stable-isomer ε CSV into
**Stable ε reference (optional)**. Use an AutoQY export from your repeated
stable-isomer measurements.

The dropped file supplies the known ε values and their errors.
The initial and PSS absorbance curves still come from the main window.

Choose **Replicate SD** or **SEM of the mean**. The reference file must contain
the error type you choose. AutoQY uses that error and the NMR composition error
to calculate the metastable ε limits. It does not add a PSS measurement error.

Both UV curves use the same scale: the stable absorbance peak is matched to
the reference ε peak. Only wavelengths shared by the UV and ε files are used.
The plot shows ε curves and error bands, not normalized absorbance curves.

Clear the ε reference to return to absorbance-only reconstruction.

### Check the result

Small negative values near zero can come from noise. Large negative peaks are
a reason to check the spectra, baseline, and PSS percentage. ε values below
−500 M⁻¹ cm⁻¹ stop export.

A very small metastable percentage makes the subtraction sensitive to noise.
Changing the data, spectrum choices, or composition requires another click on
**Reconstruct**. Changing the main treatment updates the result automatically.

## 6. Look at one wavelength over time

1. Load a time-ordered series.
2. Open **Wavelength slice over time** below the main plot.
3. Enter a wavelength, for example `450` nm.
4. Set **Seconds per timestamp**.

Leave the multiplier at `1` if the timestamps are already in seconds.
If the points are numbered 0, 1, 2… and you measured every 30 seconds, enter `30`.

**Fit exponential decay** fits a simple exponential and reports a lifetime, τ.
This is not a quantum-yield calculation. A lifetime is not a half-life:
for a first-order process, half-life = 0.693 × τ.

Check that the fitted line follows the data. If the measurement ends before
one lifetime, the lifetime and final plateau can be poorly determined.

## 7. Change names and colors

Above the main plot:

- **Legend options → Enable legend editing** lets you edit names in **Loaded spectra**.
- **Hide all** hides legend entries, not the curves. Use the **Show** controls to restore the entries you want.
- **Axis names** changes the main and time-slice axis titles.
- **Minimal colors** uses blue, orange, black, and grey. It also covers mean ε, reconstructed curves, and fits.

Without minimal colors, AutoQY uses a color-blind-friendly palette.
Colors repeat after seven spectra.

Edited names are used in saved CSV column headings. Your original files are
not renamed.

## 8. Read or copy the numbers

Open **Displayed data** below the plots. Choose the plot you want to inspect,
then scroll through its rows and columns.

- Check rows or columns and click **Copy selection**.
- If both are checked, only the cells where they meet are copied.
- For a block of cells, click the first cell and Shift-click the last.
  Use **Copy selection** or Ctrl+C.
- **Copy all** copies the whole displayed sheet, including rows below the screen.
- Paste into Excel with Ctrl+V.

The sheet cannot be edited. It follows the plotted curves, so hidden curves
and spectra left out of a large-series preview are not shown here.

## 9. Save data, figures, or a session

CSV saves the numbers in a file you can open in Excel. PNG and SVG save pictures.

Open **5 · Output**. Choose **Content to save** once:

| Choice | What is saved |
| --- | --- |
| Processed dataset / calculated ε | Treated absorbance; also ε and its errors when all concentrations and path lengths are entered |
| Reconstructed metastable spectrum | The reconstructed absorbance or ε result |
| Wavelength slice | The chosen time trace; its CSV also includes the fit, if enabled |
| Stable + reconstructed ε (pair) | Two ε CSVs; the figure shows the ε reconstruction |

Click **Save CSV**, or open **Save figures** and click **Save PNG** or **Save SVG**.
Choose the folder and name in Save As. The CSV pair asks for two file names.
The dialog asks before replacing an existing file.

CSV export includes the full treated dataset, even when the plot shows only
60 spectra. It can also include errors and other values not drawn in the figure.

Negative values are kept by default. **Convert negative absorbance and ε values
to 0** changes the saved CSV, not the plot. The reconstructed ε CSV keeps a
separate column with the original calculated values.

There is one **Image export options** menu for all figures. It controls titles,
legends, axis labels, the grid, and whether the y axis starts at zero.
**Origin-style export** makes larger text and thicker lines on a white background.
PNG is a fixed-size image; SVG stays sharp when resized.

### Save your work so you can continue later

Click **Save session (.aqys)** and choose a name in Save As.
This saves the spectra and settings, not just the plot.

To reopen it, drop one `.aqys` file into **1 · Data**, or choose it with
**Choose files**. It replaces the current session. Do not drop a session
and spectra together.

AutoQY also saves the latest state automatically. Wait for the autosave message
before closing. Clearing the spectra clears the autosaved state, but does not
delete a session file you saved yourself.

## 10. Set up a quantum-yield analysis

Return to **AutoQY Analysis**. Work through its sections in order.

### 1 · Project

Set **JSON base folder** to your experiment folder.
A JSON file records the settings and, with **Include all input data in JSON**
enabled, the complete input data. You can send that single file to another person
without also sending the original spectra, reference files, or timestamps.

If you already have `analysis.json`, load it instead of filling everything again.
Keep related input files in the same experiment folder.

### 2 · Identity and data

Give the analysis a name, then choose:

1. Measurement spectra: the UV–Vis series recorded during irradiation.
2. Reactant ε file.
3. Product ε file.
4. LED emission file, unless using an actinometer.
5. Irradiation timestamps.

The times must match the spectra in order. Use the correct file-format setting;
**Generic CSV** is suitable for the AutoQY CSV exports.

You need both ε files. An absorbance-only reconstruction is not a replacement
for the product ε file.

For the LED file, open **Processing**. Check the wavelength range, smoothing,
and baseline against your measured LED curve. These controls treat the LED,
not the absorbance spectra. The smoothing window here is in detector points,
not nm. The tutorial gives settings for its example LED.

### 3 · Experiment

Enter your measured values, with the units shown beside each field:

- sample volume in µL;
- cuvette path length in cm;
- power and power error in mW;
- irradiation wavelength in nm;
- thermal reaction rates in s⁻¹.

Thermal rates describe reactions without light. For a first-order reaction,
rate = 0.693 / half-life, with the half-life in seconds.
Do not enter a half-life into a rate field.

With power input, AutoQY uses the full LED spectrum to calculate the number of
photons. The wavelength field does not replace the LED file.

If a chemical actinometer measured your photon flux, tick **Use photon flux
from a chemical actinometer**. Enter flux and its error in **mol photons/s**,
and the irradiation wavelength. This mode uses one wavelength and does not
need an LED file.

### 4 · Fit

Start with **Regularized concentrations**. This finds the amounts of reactant
and product at each time, discourages jumps between nearby points, then fits
the reaction rates to calculate quantum yields.

Compare it with **Full-spectrum ODE absorbance**, which fits the changing
spectra directly to the reaction model.

The two legacy methods are kept for older analyses. You do not need them to
follow the tutorial.

Leave **Expert optimizer settings** at their defaults to start:
initial yields `0.5` and bounds `0–1`. These are starting guesses, not results.
Under **Method-specific controls**, the regularization strength starts at `1`.
If you measured PSS by NMR, enter the **metastable photoproduct** percentage in
**Expected metastable photoproduct at PSS (%)**. For 23% stable reactant and
77% metastable product, enter **77**. This uses the same percentage as Spectral
Treatment and is an optional diagnostic reference, not a fitting constraint.

**Compare fit methods** compares results without ε errors and does not save
result files. It does not replace your final analysis with uncertainty enabled.

### 5 · Uncertainty

Choose **ε range** if both ε files contain AutoQY error values.
Choose **Standard deviation (SD)** to use the differences between repeated
samples. Use SEM only if you want the error on their average.

AutoQY repeats the calculation with lower, central, and upper ε values.
Reconstructed product files can have different upper and lower errors.
The results list the fit/power error, ε error, and combined error separately.

### 6 · Output and 7 · Analyze

Choose where to save the analysis results and which files to write.
These Analysis settings are separate from Spectral Treatment's Save As buttons.

Leave **Include all input data in JSON (recommended)** enabled for a portable
experiment. Click **Save JSON**, then **Run analysis**. To recover or share the
experiment, load that JSON on any computer with this AutoQY build, choose a local
results directory, and run. The original input files are not needed. Disable the
option if you prefer the older file-link format.
Running does not save `analysis.json` for you.
Only enable overwrite when you want to replace an earlier run.

## 11. Check the results before using them

A small reported error does not prove the analysis is right.

Check these plots:

- **Concentrations**: do the amounts change sensibly over time?
- **Preprocessing**: does the LED emit where the compounds absorb?
- **Endpoint reconstruction**: does the calculated final spectrum match the measured one?
- **Absorbance residuals**: these are measured minus fitted absorbance.
  Clear bands or patterns can point to wrong references, a bad baseline, or another absorbing compound.

The regularized method uses a simple exponential to keep concentration points
from jumping around. That exponential is not the photochemical model.
Small patterns in its **Fraction residual** plot can occur for this reason.
Check them alongside the other plots, not on their own.

Compare the two main fit methods and the fitted PSS with your NMR result.
If they disagree badly, or small changes in treatment cause large changes in
quantum yield, investigate before reporting the result.
