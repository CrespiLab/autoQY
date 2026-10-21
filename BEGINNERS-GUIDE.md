# AutoQY beginner guide

Use this guide for your own data. To practise with the included files, follow the
[395 nm tutorial](ExampleData/Example-4_395nm-EpsilonError/TUTORIAL.md).

AutoQY has three windows:

- **Spectral Treatment**: clean spectra, calculate ε, reconstruct a spectrum, and save data or figures.
- **Analysis**: calculate quantum yields from your measurements.
- **Power**: process power-monitor measurements before entering power in Analysis.

Click an **ⓘ** button for help with the control beside it. Press Escape to close it.

## A few terms used in the buttons

- **Absorbance**: the value measured by your UV–Vis instrument.
- **ε** (epsilon, or molar absorptivity): absorbance divided by concentration and cuvette path length. This lets you compare how strongly compounds absorb light.
- **Reactant (R)**: the isomer converted in the chosen reaction direction.
- **Product (P)**: the isomer formed in that direction.
- **PSS**: the mixture reached when continued irradiation no longer changes its composition.
- **Quantum yield (Φ)**: how many molecules change per absorbed photon (a particle of light). A value of 0.25 means 25%.
- **SD**: how much repeated measurements differ from each other.
- **SEM**: the uncertainty estimated for the average of repeated measurements.

For a back reaction, R can be the metastable isomer and P the stable isomer.
An analysis can start from a PSS containing both.

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
The file chooser and all Save As dialogs open in the latest folder you chose,
including when loading an ε reference or saving a figure.
Loading more files adds them to the current set. **Clear all spectra** clears
the loaded set. Source files remain in their original folder.

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

Open **Preprocess spectra** and enable the treatments your data need.

### Baseline

Use this if the spectrum sits above or below zero where there should be no absorption.

1. Choose a wavelength interval with no real absorption.
2. Enter its start and end.
3. Turn **Baseline** on.

For some samples, 600–650 nm works. Check your own spectrum first.
Choose a region outside the absorption bands to preserve the signal.

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
SVD off. Check the suggested component count against the spectral shapes and residuals.

Original spectra are hidden by default in the treated-absorbance plot.
Use **Show original** to compare them with the treated curves. This also affects
the saved absorbance figure. Source files stay unchanged.

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
its concentration, dilution, cuvette, and baseline before choosing the smoothing settings.

Concentration fields are available for sets of up to 60 spectra.
Use a separate set of known-concentration samples to calculate ε.

## 5. Reconstruct the pure metastable spectrum

For a measured pure metastable form, calculate ε directly as in step 4.

Use this when you have a pure stable-isomer spectrum and a PSS mixture, but
cannot measure the pure metastable form.

You need:

- the initial stable and PSS absorbance spectra;
- the percentage of metastable isomer in the PSS, measured separately, for example by NMR.

Both UV–Vis measurements must use the same total concentration and path length.
The calculation assumes two absorbing forms and uses the relative sizes of
the measured absorbance curves. Preserve those sizes during treatment.

1. Load the UV–Vis series in **1 · Data**.
2. Set the range, baseline, and smoothing in the main controls.
3. Open **4 · Reconstruction → Reconstruct spectrum**.
4. Select **Stable spectrum** and **PSS spectrum**.
   The first and last spectra are selected by default; check them.
5. Enter **Metastable (%)**.
   For a mixture with 23% stable and 77% metastable isomer, enter **77**.
6. Enter the composition error in percentage points.
   The default `1` means one standard uncertainty of 1 percentage point.
   For a 77% composition, that corresponds to 76–78% at ±1 SD.
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
to calculate the metastable ε limits. The PSS absorbance is held fixed.

Both UV curves use the same scale: the stable absorbance peak is matched to
the reference ε peak. Only wavelengths shared by the UV and ε files are used.
The plot shows the calibrated ε curves and their error bands.

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
For a first-order process, half-life = 0.693 × τ.
Use Analysis to calculate quantum yields from the full reaction model.

Check that the fitted line follows the data. If the measurement ends before
one lifetime, the lifetime and final plateau can be poorly determined.

## 7. Change names and colors

Above the main plot:

- **Legend options → Enable legend editing** lets you edit names in **Loaded spectra**.
- **Hide all** hides the legend entries while keeping the curves visible. Use **Show** to restore individual entries.
- **Axis names** changes the main and time-slice axis titles.
- **Minimal colors** uses blue, orange, black, and grey. It also covers mean ε, reconstructed curves, and fits.

Without minimal colors, AutoQY uses a color-blind-friendly palette.
Colors repeat after seven spectra.

Edited names become the saved CSV column headings. Source filenames stay unchanged.

## 8. Read or copy the numbers

Open **Displayed data** below the plots. Choose the plot you want to inspect,
then scroll through its rows and columns.

- Check rows or columns and click **Copy selection**.
- If both are checked, only the cells where they meet are copied.
- For a block of cells, click the first cell and Shift-click the last.
  Use **Copy selection** or Ctrl+C.
- **Copy all** copies the whole displayed sheet, including rows below the screen.
- Paste into Excel with Ctrl+V.

The sheet is a read-only view of the plotted curves. CSV export contains the
complete treated dataset, including spectra outside a large-series preview.

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
60 spectra. It also includes the available error columns and calculation data.

Negative values are kept by default. **Convert negative absorbance and ε values
to 0** applies to the saved CSV. The plot keeps its current values, and the reconstructed ε CSV keeps a
separate column with the original calculated values.

There is one **Image export options** menu for all figures. It controls titles,
legends, axis labels, the grid, and whether the y axis starts at zero.
**Origin-style export** makes larger text and thicker lines on a white background.
PNG is a fixed-size image; SVG stays sharp when resized.

### Save your work so you can continue later

Click **Save session (.aqys)** and choose a name in Save As.
The session contains the spectra and settings.

To reopen it, drop one `.aqys` file into **1 · Data**, or choose it with
**Choose files**. Load the session on its own; it replaces the current session.

AutoQY also saves the latest state automatically. Wait for the autosave message
before closing. Clearing the spectra clears the autosaved state.
Your manually saved session files remain available.

## 10. Set up a quantum-yield analysis

Return to **AutoQY Analysis**. Work through its sections in order.

### 1 · Project

Set **JSON base folder** to your experiment folder.
A JSON file records the settings and, with **Include all input data in JSON**
enabled, the complete input data. You can send that single file to another person
without also sending the original spectra, reference files, or timestamps.

Load an existing `analysis.json` to restore its saved settings and data.
Keep related input files in the same experiment folder.

To get CSV copies of the data stored in a JSON, load it and click **Restore input
files as CSV…** here. Choose where to save them. ε error columns and timestamp
events are kept. Existing files are kept too; new copies get numbered names.
The button becomes available when the JSON contains the input data itself.

### 2 · Identity and data

Give the analysis a name, then choose:

1. Measurement spectra: the UV–Vis series recorded during irradiation.
2. Reactant ε file.
3. Product ε file.
4. LED emission file, unless using an actinometer.
5. Irradiation timestamps.

The times must match the spectra in order. Use the correct file-format setting;
**Generic CSV** is suitable for the AutoQY CSV exports.

Analysis requires an ε spectrum for each isomer. Add a stable ε reference in
Spectral Treatment to calibrate a reconstructed absorbance spectrum.

For the LED file, open **Processing**. Check the wavelength range, smoothing,
and baseline against your measured LED curve. These controls apply to the LED.
The smoothing window is in detector points. Prepare the absorbance spectra in
Spectral Treatment, where the smoothing window is in nm.

### 3 · Experiment

Enter your measured values, with the units shown beside each field:

- sample volume in µL;
- cuvette path length in cm;
- power and power error in mW;
- irradiation wavelength in nm;
- thermal reaction rates in s⁻¹.

Thermal rates describe reactions without light. For a first-order reaction,
rate = 0.693 / half-life, with the half-life in seconds.

With power input, AutoQY uses the full LED spectrum to calculate the number of
photons. The wavelength field records the nominal irradiation wavelength.

If a chemical actinometer measured your photon flux, tick **Use photon flux
from a chemical actinometer**. Enter flux and its error in **mol photons/s**,
and the irradiation wavelength. This mode models monochromatic irradiation
using the supplied photon flux.

### 4 · Fit

New analyses use **Full spectrum ODE** by default. It fits the changing
spectra directly to the reaction model, including the starting mixture and
total concentration. Loading a JSON keeps that analysis's saved method.

Compare it with **Concentrations (legacy pure NNLS)**, which finds reactant
and product amounts independently at each time, then fits the reaction model. Check the concentration curves and spectral residuals
alongside the reported yields. **Emission (legacy)** is also
available for reproducing older analyses.

Leave **Expert optimizer settings** at their defaults to start:
initial yield guesses `0.5` and allowed yield bounds `0–1`.
If you measured PSS by NMR, enter the **metastable photoproduct** percentage in
**Expected metastable photoproduct at PSS (%)**. For 23% stable reactant and
77% metastable product, enter **77**. This uses the same percentage as Spectral
Treatment. AutoQY compares the fitted PSS with this optional reference value.
The check uses the percentage of the species assigned to P in Analysis.

**Compare fit methods** uses the same inputs and uncertainty settings as your
analysis. The comparison is displayed in the window.
The displayed yield errors therefore include the configured ε uncertainty;
leave that option off only when you want to compare fit/power errors alone.

### 5 · Uncertainty

Choose **ε ranges** if both ε files contain AutoQY error values.
Choose **Standard deviation (SD)** to use the differences between repeated
samples. Use SEM only if you want the error on their average.

AutoQY repeats the calculation with lower, central, and upper ε values.
Reconstructed product files can have different upper and lower errors.
The results list the fit/power error, ε error, and combined error separately.

You can keep **ε ranges**, or choose **Monte Carlo**.
Monte Carlo uses the uncertainties already supplied and keeps the PSS absorbance
used for reconstruction fixed. It uses
replicate ε data and, for a reconstructed isomer, the matching reconstruction
with its NMR composition uncertainty. Re-export from Spectral Treatment if an
older file is missing those details.

Monte Carlo shows the mean ± standard uncertainty (1 SD). The error uses
one significant digit, or two when it starts with 1 or 2, and the mean uses
the same decimal place. The default maximum is 2048 simulations. AutoQY checks
the precision of the mean and error as it runs and can stop early. When the
check passes, the remaining simulation error is small enough for the displayed
precision. If it reaches the limit first, **Technical details** shows how much
simulation error remains.

The coverage interval has a separate endpoint check in **Technical details**.
Its limits can remain approximate even when the mean and error pass their check.
Failed fits and yields reaching their allowed limits also appear in the results.

A live counter shows how many simulations have finished and how many were reused.
Each completed simulation is saved automatically in `.autoqy-monte-carlo`
inside the output folder. To continue a run, keep the same output folder,
inputs, fit settings, and random seed, raise **Maximum simulations**, and
press **Run analysis**. For example, increasing 256 to 4096 reuses the
first 256 simulations. You can close and reopen the GUI between runs.
An interrupted run continues from its last completed simulation.

Changes to the input data, uncertainty inputs, fit method, or seed start
a separate saved run. Coverage and precision settings can reuse the same
draws. Changes to sampling or fitting code start a separate saved run. Updates to
messages or stopping checks keep compatible saved simulations.
To update existing result files, enable **Allow replacing existing results** or
choose a new **Output stem** while keeping the output folder.

The saved random seed
lets you repeat the same calculation. A large interval can still reflect limited
experimental information. Check its size, stability, and fit diagnostics together.

### 6 · Output and the analysis buttons

Choose where to save the analysis results and which files to write.
These Analysis settings are separate from Spectral Treatment's Save As buttons.

Leave **Include all input data in JSON (recommended)** enabled for a portable
experiment. Click **Save JSON**, then **Run analysis**. To recover or share the
experiment, load that JSON on any computer with this AutoQY build, choose a local
results directory, and run. The JSON carries the input data with it. Disable
the option to save file links. Use **Save JSON** whenever you want to save the
editable experiment settings.
Only enable overwrite when you want to replace an earlier run.

## 11. Check the results before using them

Assess the reported errors together with the fit and residual plots.

Check these plots:

- **Concentrations**: do the amounts change sensibly over time?
- **Preprocessing**: does the LED emit where the compounds absorb? Each valid
  input appears as you load it, even while other inputs are missing or need correction.
- **Endpoint reconstruction**: does the calculated final spectrum match the measured one?
- **Absorbance residuals**: these are measured minus fitted absorbance.
  Clear bands or patterns can point to wrong references, a bad baseline, or another absorbing compound.

Concentration plots use blue for reactant and orange for product. Open circles
show concentrations recovered from the spectra; solid lines show the fitted
reaction model. **ε ranges** shows the variation across the tested ε curves.
**Monte Carlo** reports mean ± standard uncertainty; detailed data includes the coverage intervals. The PNG and SVG exports identify the
two yields with the reaction directions Φ R→P and Φ P→R.
The shaded band around the fraction residual shows a pointwise Monte Carlo
coverage interval (95% by default), or the variation across the tested ε curves
when using **ε ranges**.

Older reconstructed ε exports may lack the saved NMR composition error. If
Monte Carlo reports missing reconstruction information, reopen the stable
and PSS absorbance spectra in **Spectral Treatment → 4 · Reconstruction**,
add the stable ε reference, enter **Composition error (± percentage points)**,
and reconstruct and save the ε CSV again. Load the new CSV in Analysis.
**ε ranges** remains available for files that contain the lower and upper ε curves.

The **Fraction residual** is measured reactant fraction minus fitted fraction.
The absorbance residual map uses a color scale centered on zero, so opposite
colors indicate underprediction and overprediction. Inspect patterns as well
as the size of the residuals.

Compare **Concentrations (legacy pure NNLS)** with **Full spectrum ODE** and
compare the fitted PSS with your NMR result.
If they disagree badly, or small changes in treatment cause large changes in
quantum yield, investigate before reporting the result.
