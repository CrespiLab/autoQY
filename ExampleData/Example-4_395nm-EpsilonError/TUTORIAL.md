# Tutorial: 395 nm analysis with deterministic epsilon uncertainty

This tutorial begins with three independently prepared reactant solutions and
finishes with a quantum-yield fit that propagates wavelength-resolved errors in
both reactant and product molar absorptivity.

## Files and experimental values

| File | Purpose | Concentration |
| --- | --- | ---: |
| `source_data/2107362U1_0001A.Abs8` | Reactant measurement 1 | 6.82e-5 M |
| `source_data/2107362U1_0001B.Abs8` | Reactant measurement 2 | 6.87e-5 M |
| `source_data/2107362U1_0001C.Abs8` | Reactant measurement 3 | 6.80e-5 M |
| `source_data/reactant_pss.dat` | Pure reactant and final PSS spectra | 23% reactant at PSS |
| `led_emission.dat` | Supplied `LED395nm_baselined.dat`, under a portable name | — |
| `measurement_spectra.dat` | Time-resolved absorbance | — |
| `timestamps.csv` | Irradiation times | — |

The optical path length is 1 cm. The irradiation power is 1.46 ± 0.03 mW.
The sample volume is 1995 µL and the thermal back-reaction rate is 6.3e-5 s⁻¹.

The ready-to-run `reactant_epsilon.tsv` and `product_epsilon.tsv` are bundled
for comparison. The following steps regenerate them from the source spectra.

## 1. Open Spectral Treatment from the Analysis GUI

Start **AutoQY Analysis**, expand **1 · Project**, and select
**Open Spectral Treatment**. Spectral Treatment opens separately so its exports
can then be selected in the Analysis GUI.

## 2. Calculate reactant molar absorptivity

1. Under **1 · Data**, select **Open files from folder** and load the three
   `.Abs8` files together, in numeric order.
2. Under **2 · Range**, set the wavelength range to **250–700 nm**.
3. Expand **Preprocess spectra**, enable **Baseline**, and set its interval to
   **600–650 nm**.
4. Select **SavGol**, then set a **5 nm** window and polynomial order **3**.
5. Leave **SVD off**. These are independent preparations, not an ordered time
   series; SVD would mix their genuine between-sample variation.
6. Under **3 · Beer–Lambert**, enter the concentrations in the same file order:
   **6.82e-5**, **6.87e-5**, and **6.80e-5 mol/L**. Enter **1 cm** for every
   path length.
7. Inspect the individual epsilon traces, mean, and shaded error band. Large
   structured differences between measurements should be investigated before
   export; the band is not a substitute for inspecting the spectra.
8. Under **4 · Output**, choose the Example 4 folder and select
   **Save reactant epsilon TSV**. The file contains processed absorbance,
   individual epsilon curves, their mean, SD, SEM, and non-negative bounds.

The bundled result is `reactant_epsilon.tsv`. If the GUI proposes a different
name, select or rename the export before using it in the analysis configuration.

## 3. Derive product epsilon from the PSS composition

1. Keep the calculated reactant epsilon loaded in Spectral Treatment.
2. Expand **5 · Optional — NMR-guided PSS subtraction** and load
   `source_data/reactant_pss.dat`. Its first spectrum is the pure reactant and
   its last spectrum is the final PSS.
3. Under **Preprocess reactant and PSS**, apply the same **600–650 nm** baseline
   and **Savitzky–Golay 5 nm / order 3** settings.
4. Enter **23** for **Reactant in final PSS (%)**.
5. Enter the experimental NMR uncertainty. The bundled product file uses the
   GUI default of **1%**. Change this value when a different uncertainty is
   justified by the composition measurement.
6. Inspect the normalized reactant and PSS curves and the reconstructed product
   epsilon. The calculation is:

   `product = (PSS − 0.23 × reactant) / (1 − 0.23)`

7. Review any negative-value warning. Small negative results down to
   −500 M⁻¹ cm⁻¹ remain visible for diagnosis; values below that threshold stop
   export. By default, the primary exported product epsilon is constrained to
   zero while the raw audit column is retained.
8. Select **Save reactant + NMR-derived epsilon TSVs**. Use the product export as
   `product_epsilon.tsv` and retain the reactant export as
   `reactant_epsilon.tsv`.

The product bounds combine the reactant measurement SD with the selected NMR
composition error. They are therefore asymmetric and wavelength dependent.

## 4. Build the quantum-yield analysis

Return to **AutoQY Analysis**. You may load the bundled `analysis.json` and
inspect each section, or enter the following values manually.

### Project and data

- JSON base folder: the `Example-4_395nm-EpsilonError` directory.
- Analysis ID: `example_395nm_epsilon_error`.
- Measurement spectra: `measurement_spectra.dat`.
- Reactant molar absorptivity: `reactant_epsilon.tsv`.
- Product molar absorptivity: `product_epsilon.tsv`.
- LED emission: `led_emission.dat`.
- Irradiation timestamps: `timestamps.csv`.
- Spectral formats: **SpectraGryph / AutoQY TSV**.
- Timestamp format: **AHK CSV**.

Store input paths relative to the JSON folder so the example remains portable.

### LED processing

Use **250–800 nm**, Savitzky–Golay window **12 points**, polynomial order **3**,
and baseline correction with exclusion multiplier **10**. This processing is
only for the LED spectrum; it is distinct from treating absorbance spectra.

### Experiment

- Sample volume: **1995 µL**.
- Path length: **1 cm**.
- Power: **1.46 mW**.
- Power error: **0.03 mW**.
- Irradiation wavelength: **395 nm**.
- Thermal back-reaction: **6.3e-5 s⁻¹**.

The irradiation wavelength is metadata and a consistency marker. Photon flux
uses the complete processed LED emission spectrum rather than a single value.

### Kinetic model

Select **Regularized concentrations** with:

- Initial Φ R→P: **0.5**.
- Initial Φ P→R: **0.5**.
- Bounds: **0–1**.
- Concentration regularization strength: **1**.
- Expected reactant at PSS: **23%**.

The other method-specific values can remain at emission threshold **0.01**,
absorbance baseline order **1**, and robust-loss scale **0.02**. They are stored
for reproducibility but only affect methods that use them.

Available fitting methods are:

- **Concentrations**: independent NNLS spectral decomposition at every time,
  followed by the kinetic fit.
- **Regularized concentrations**: conserves total concentration and softly
  regularizes fractions to a free exponential envelope; used in this example.
- **Full-spectrum ODE absorbance**: jointly fits kinetics and the complete
  absorbance matrix, with optional spectral baselines.
- **Emission (legacy)**: fits absorbance only in the active LED band and assumes
  initially pure reactant.

Use **Compare fit methods** to inspect nominal-method sensitivity. Comparison
runs without epsilon uncertainty and writes no result files; it does not replace
the selected uncertainty-enabled analysis.

### Deterministic epsilon uncertainty

Under **5 · Uncertainty** choose:

- Epsilon propagation: **Deterministic wavelength bounds**.
- Repeat-spectrum error metric: **Standard deviation (SD)**.

SD represents the observed variation among the three independent preparations.
SEM is also available, but represents uncertainty in their mean and is smaller;
select it only when that interpretation matches the scientific question.

AutoQY evaluates lower, nominal, and upper wavelength-resolved epsilon curves
for both species. The reactant uses mean ± SD with non-negative bounds. The
NMR-derived product uses its asymmetric exported bounds. The result separates
optimizer/power, epsilon, and combined quantum-yield uncertainty contributions.

### Output and validation

Set the output directory to `results`, use stem
`Results_395nm_EpsilonError`, enable text, figures, JSON, configuration, and
detailed-data exports, and enable overwrite only when replacing previous runs
is intentional.

1. Select **Validate** and inspect the preprocessing plot: initial spectrum in
   blue, final spectrum in orange, intermediate spectra in translucent grey,
   and the LED emission overlaid.
2. Review amber or red diagnostics through their information popovers. A red
   diagnostic opens the panel automatically.
3. Select **Run analysis** and inspect concentrations, fraction residuals,
   reference/reconstruction plots, and the wavelength-resolved residual map.

For the bundled files, the expected result is approximately:

- Φ R→P = **25.0 ± 1.2%**.
- Φ P→R = **21 ± 3%**.
- Extrapolated PSS = **23.3% reactant**.

Small differences can result if preprocessing or uncertainty settings differ.

## 5. Terminal equivalent

From the repository root:

```text
autoqy-core validate ExampleData/Example-4_395nm-EpsilonError/analysis.json
autoqy-core run ExampleData/Example-4_395nm-EpsilonError/analysis.json
```
