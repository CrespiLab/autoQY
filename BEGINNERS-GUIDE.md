# AutoQY Spectral Treatment and Analysis GUI

## Beginner guide: from UV–Vis spectra to spectral treatment, simple kinetics, and quantum-yield analysis

AutoQY contains two main graphical interfaces:

- **Spectral Treatment** — for preparing UV–Vis spectra, calculating molar absorptivity, inspecting kinetic traces, fitting simple exponential lifetimes, and exporting figures.
- **Analysis GUI** — for combining spectral evolution, molar absorptivities, irradiation times, LED emission, optical power, and experimental parameters to calculate photoisomerization quantum yields.

A typical quantum-yield workflow is:

```text id="o290th"
Reference UV–Vis spectra
        ↓
Spectral Treatment
        ↓
Reactant ε + Product ε
        ↓
Analysis GUI
        ↓
Φ(R→P) + Φ(P→R)
```

Spectral Treatment can also be used independently for routine UV–Vis work, including:

- baseline correction;
- smoothing;
- comparing spectral series;
- calculating molar absorptivity;
- estimating uncertainty from replicate ε measurements;
- reconstructing a photoproduct spectrum from an NMR-characterized PSS;
- following absorbance at one wavelength through a kinetic series;
- fitting a simple exponential lifetime;
- preparing spectral figures;
- simplifying legends;
- exporting PNG or SVG figures.

In most cases, the complete workflow can be performed directly from the GUIs without manually editing Python code.

---

# 1. What do I need for a quantum-yield experiment?

For a two-state photochemical reaction:

```math id="5dp8j1"
R \rightleftharpoons P
```

the analysis normally requires:

| Input | Meaning |
|---|---|
| Measurement spectra | UV–Vis spectra recorded during irradiation |
| Reactant ε | Molar absorptivity spectrum of R |
| Product ε | Molar absorptivity spectrum of P |
| LED spectrum | Measured emission spectrum of the irradiation LED |
| Irradiation timestamps | Irradiation time corresponding to every spectrum |

You will also need:

- sample volume;
- optical path length;
- irradiation power;
- uncertainty in irradiation power;
- thermal R → P rate, if relevant;
- thermal P → R rate, if relevant.

With optical-power input, the irradiation wavelength is mainly metadata and a
consistency marker. The photon-flux calculation uses the **complete processed
LED emission spectrum**, rather than treating the LED as perfectly
monochromatic.

If the photon flux was measured with a chemical actinometer, enable **Use
photon flux from a chemical actinometer** instead. Enter the flux in mol
photons/s and its uncertainty. An LED file is then unnecessary, and the
irradiation wavelength is mandatory because AutoQY uses it as the exact
monochromatic wavelength.

---

# 2. Starting AutoQY

After installation, open:

```text id="7vfmmg"
AutoQY Analysis
```

Inside the Analysis GUI, expand:

```text id="946y1g"
1 · Project
```

and click:

```text id="gfncp9"
Open Spectral Treatment
```

Spectral Treatment opens in a separate window.

---

# 3. Spectral Treatment: loading spectra

Open:

```text id="pes0i0"
1 · Data → Spectral data
```

You can:

- drag and drop spectra;
- choose individual files;
- use **Open files from folder**.

Supported formats include:

- SpectraGryph text `.dat`;
- Analytik Jena SPECORD WinASPECT binary `.dat`;
- Agilent/Varian Cary 50 and 60 `.DSW` and `.BSW`;
- Avantes `.Abs8`;
- TSV;
- CSV.

Multiple spectra can be loaded together.

## Check the spectrum order

For molar-absorptivity measurements, the concentrations entered later should correspond to the correct spectra.

For kinetic measurements, the spectra should normally be arranged in chronological order.

Expand:

```text id="qfkwcj"
Loaded spectra: order, legend, removal
```

Here you can:

- move individual spectra up or down;
- remove spectra;
- choose which traces appear in the legend;
- rename displayed legend labels.

Changing the displayed name does **not** rename the original file.

## Large kinetic datasets

When more than 60 spectra are loaded, AutoQY keeps the browser responsive by
showing at most 60 evenly spaced spectra in the interactive plot. The first and
last spectra are always included. This is only a preview limit: baseline
correction, smoothing, SVD, wavelength slices, and processed-data export still
use every loaded spectrum.

For these large datasets, **Loaded spectra** shows only the first and last
legend-name and **Show** controls. Intermediate order and removal controls are
not created, and only the endpoint spectra can appear in the legend. If the
order needs correcting, arrange the source data before loading it.

---

# 4. Select the wavelength range

Open:

```text id="ow7t7b"
2 · Range → Wavelengths
```

Choose the spectral range that is useful for the experiment.

For example:

```text id="arh3g5"
250–700 nm
```

It is usually useful to exclude detector regions that contain little information or are particularly noisy.

The selected range is used for:

- displayed spectra;
- molar-absorptivity calculations;
- uncertainty calculations;
- exported data.

---

# 5. Preprocess the spectra

Inside **2 · Range**, expand:

```text id="235zlr"
Preprocess spectra
```

## Baseline correction

Enable:

```text id="hmudvw"
Baseline
```

and select a wavelength interval where the compound is expected to have little or no absorbance.

For example:

```text id="y6krpc"
600–650 nm
```

if that region is appropriate for the molecule being studied.

Number fields used for wavelength selection and preprocessing are applied when
you press **Enter** or move to another control. This prevents a large spectrum
file from being recalculated once for every digit typed.

After baseline correction, the processed-absorbance preview shows only the
corrected spectra by default. Before supplying a complete set of Beer–Lambert
values, enable **Show original** above the plot to compare the raw and corrected
curves. Leave it disabled for a cleaner absorbance figure; the same choice is
used by PNG and SVG export.

## Savitzky–Golay smoothing

Enable:

```text id="fwxfnu"
SavGol
```

A reasonable starting point for ordinary UV–Vis spectra is:

```text id="h8r0po"
Window: 5 nm
Polynomial order: 3
```

A useful check is to compare the spectrum with smoothing turned on and off. The aim is to reduce noise without noticeably changing the underlying band shape.

<details>
<summary><strong>More about baseline correction and smoothing</strong></summary>

A baseline interval should ideally correspond to a region in which the compound has negligible absorbance.

If the selected baseline interval contains a real absorption band, subtraction can distort the spectrum.

Similarly, smoothing should normally remove high-frequency noise rather than alter peak positions or broad band shapes.

If changing the Savitzky–Golay settings visibly changes the molecular absorption profile, it is worth using a smaller window or leaving smoothing off.

</details>

---

# 6. What about SVD?

For independently prepared molar-absorptivity measurements, it is usually preferable to leave SVD off.

```text id="5x7u4y"
Independent replicate solutions
        ↓
SVD usually OFF
```

```text id="qaefce"
Ordered spectral time series
        ↓
SVD may be useful
```

<details>
<summary><strong>Why is SVD usually avoided for independent ε replicates?</strong></summary>

Independent preparations contain useful information about experimental reproducibility.

Applying SVD across those spectra can reduce the apparent variability between replicates.

For time-series data, SVD can instead be useful for reducing noise because the spectra belong to the same evolving experiment.

It is still useful to inspect the untreated spectra first.

</details>

---

# 7. Calculate molar absorptivity

Open:

```text id="3uvumb"
3 · Beer–Lambert → Concentrations
```

For each spectrum enter:

- concentration in mol/L;
- optical path length in cm.

For example:

| Spectrum | Concentration | Path length |
|---|---:|---:|
| 1 | 6.82e-5 M | 1 cm |
| 2 | 6.87e-5 M | 1 cm |
| 3 | 6.80e-5 M | 1 cm |

AutoQY applies:

```math id="1zybv5"
\varepsilon_i(\lambda) = \frac{A_i(\lambda)}{c_i l_i}
```

where:

- $A_i$ is absorbance;
- $c_i$ is concentration;
- $l_i$ is optical path length.

AutoQY then calculates:

- each individual ε spectrum;
- mean ε;
- standard deviation;
- standard error of the mean.

For $n>1$:

```math id="cy4d0k"
SEM = \frac{SD}{\sqrt{n}}
```

Per-spectrum concentration and path-length fields are hidden when a dataset has
more than 60 spectra. Large time series can still be processed and exported,
but molar absorptivity should be calculated from a separate set containing at
most 60 known-concentration spectra.

---

# 8. Inspect and export ε

Before exporting, inspect the individual ε curves.

Independent preparations should usually give reasonably similar spectra.

For example:

```text id="woq28u"
Sample 1: εmax = 18,000 M⁻¹ cm⁻¹
Sample 2: εmax = 18,300 M⁻¹ cm⁻¹
Sample 3: εmax = 25,000 M⁻¹ cm⁻¹
```

would be worth checking more closely.

Possible reasons for larger differences include:

- concentration uncertainty;
- dilution error;
- baseline differences;
- cuvette differences;
- aggregation;
- decomposition;
- instrumental variation.

To export, open:

```text id="0jmzhf"
4 · Output → Export processed dataset
```

and click:

```text id="jp9yu2"
Save processed CSV
```

A typical filename is:

```text id="tx4gu6"
reactant_absorptivity.csv
```

---

# 9. Use Spectral Treatment for simple kinetics

Spectral Treatment can also be used independently of the quantum-yield workflow.

Load a time-ordered spectral series and select a wavelength of interest.

For example:

```text id="4eccpq"
450 nm
```

Spectral Treatment extracts:

```math id="ykd4km"
A(450\ \mathrm{nm},t)
```

from every spectrum and plots the resulting kinetic trace.

If the selected wavelength lies between two detector points, the value is interpolated.

## Set the time axis

Use:

```text id="2k1vmc"
Seconds per timestamp
```

If the coordinates are already seconds, leave it at:

```text id="u19bpp"
1
```

If the spectra are numbered `0, 1, 2, 3...` and one spectrum was recorded every 30 s, enter:

```text id="nayb5p"
30
```

## Fit an exponential lifetime

Enable:

```text id="gfh0n8"
Fit exponential decay
```

AutoQY fits:

```math id="gmv9yp"
A(t) = A_{\infty} + \Delta A\,e^{-t/\tau}
```

and reports the lifetime with its fit uncertainty.

<details>
<summary><strong>Lifetime, half-life, and measurement duration</strong></summary>

The fitted lifetime is:

```math id="jdf1g1"
\tau = \frac{1}{k}
```

For a first-order process:

```math id="6o1elj"
t_{1/2} = \tau \ln 2
```

or:

```math id="nh90sm"
k = \frac{\ln 2}{t_{1/2}}
```

So lifetime and half-life are related, but are not the same quantity.

AutoQY also checks whether the measured time span extends beyond approximately one fitted lifetime.

A shorter measurement can still be fitted, but the plateau and lifetime may be less well constrained.

The simple exponential fit is most appropriate for traces that are reasonably close to single-exponential behaviour.

</details>

---

# 10. Use Spectral Treatment to prepare figures

Spectral Treatment can also prepare clean spectral figures directly.

For datasets containing at most 60 spectra, you can keep all curves visible
without showing every filename in the legend. Above 60 spectra, the figure uses
the evenly spaced preview subset described earlier; the complete series remains
in the processed CSV export.

Expand:

```text id="0r67ux"
Loaded spectra: order, legend, removal
```

Click:

```text id="6mdvgj"
Hide all
```

to hide the legend entries without removing the traces.

Then re-enable only the important spectra, for example:

```text id="zpbzu8"
First spectrum
Last spectrum
```

You can also rename the displayed labels, for example:

```text id="0cunoy"
E
PSS
```

## Minimal colors

Enable:

```text id="wfmcxi"
Minimal colors
```

to highlight the initial and final spectra while keeping intermediate traces visually quieter.

<details>
<summary><strong>Example workflow for a publication-style spectral figure</strong></summary>

For a 30-spectrum irradiation series:

1. Load all spectra.
2. Put them in chronological order.
3. Select the useful wavelength range.
4. Apply preprocessing if needed.
5. Enable **Minimal colors**.
6. Click **Hide all**.
7. Re-enable only the first and final spectra in the legend.
8. Rename them, for example, `E` and `PSS`.
9. Remove the title if unnecessary.
10. Enable **Origin-style export**.
11. Save as SVG or PNG.

SVG is useful when a vector format is preferred.

</details>

---

# 11. Obtain the product molar absorptivity

There are two common cases.

## Case A — Pure product is available

Prepare pure P at known concentration and repeat the same procedure used for the reactant.

Export, for example:

```text id="ksvehe"
product_absorptivity.csv
```

## Case B — The product cannot be isolated

If irradiation produces a known PSS mixture, AutoQY can reconstruct the product spectrum using an independently measured PSS composition.

<details>
<summary><strong>NMR-guided PSS subtraction</strong></summary>

Expand:

```text id="5z9n8a"
5 · Optional → NMR-guided PSS subtraction
```

Load a UV–Vis dataset in which:

```text id="dz1n12"
First spectrum = pure reactant
Last spectrum  = final PSS
```

Enter:

```text id="n3kept"
Reactant in final PSS (%)
```

For example:

```text id="8bf35r"
23
```

If $x$ is the reactant fraction at the PSS:

```math id="cfyh3c"
P = \frac{PSS - xR}{1-x}
```

For $x=0.23$:

```math id="hg7qtu"
P = \frac{PSS - 0.23R}{0.77}
```

You can also enter the NMR uncertainty.

The resulting product ε uncertainty is generally wavelength dependent and asymmetric.

Small negative reconstructed ε values near a zero baseline can arise from noise or subtraction uncertainty.

Larger negative spectral features may suggest checking the PSS composition, baseline, normalization, or whether more than two species are present.

</details>

---

# 12. Open the Analysis GUI

The Analysis GUI contains:

```text id="k2vbwm"
1 · Project
2 · Identity and data
3 · Experiment
4 · Fit
5 · Uncertainty
6 · Output
7 · Analyze
```

Working through them in order is usually the easiest approach.

---

# 13. Project

Expand:

```text id="tocxcd"
1 · Project → JSON and tools
```

Choose the folder containing your experiment.

For example:

```text id="it3xt1"
MyExperiment/
├── measurement_spectra.csv
├── reactant_absorptivity.csv
├── product_absorptivity.csv
├── led_emission.csv
└── timestamps.csv
```

`analysis.json` stores the configuration used for the analysis, including:

- input files;
- experimental parameters;
- fitting method;
- uncertainty treatment;
- output settings.

---

# 14. Identity and data

Expand:

```text id="t1fm3u"
2 · Identity and data → Experiment files
```

Choose a meaningful analysis ID, for example:

```text id="3aehoy"
compound1_395nm
```

and an output stem such as:

```text id="zbx30h"
compound1_395nm_results
```

You can also replace generic names such as `reactant` and `product` with labels such as `E` and `Z`.

Select:

1. **Measurement spectra**
2. **Reactant molar absorptivity**
3. **Product molar absorptivity**
4. **LED emission** (unless chemical-actinometer mode is used)
5. **Irradiation timestamps**

For new data, use **Generic CSV** where possible.

---

# 15. LED processing

Inside the LED section expand:

```text id="a7vqmg"
Processing
```

A reasonable starting point is:

```text id="xznvda"
Wavelength start:              250 nm
Wavelength end:                800 nm
Savitzky–Golay window:         12 points
Polynomial order:              3
Baseline correction:           ON
Baseline exclusion multiplier: 10
```

These settings apply only to the **LED emission spectrum**.

They do not preprocess the experimental absorbance spectra.

They are ignored when chemical-actinometer photon flux is selected.

---

# 16. Experimental parameters

Expand:

```text id="u8b47c"
3 · Experiment → Physical parameters
```

Enter the experimental values.

For example:

```text id="2exvvs"
Sample volume:          1995 µL
Path length:            1 cm
Power:                  1.46 mW
Power error:            0.03 mW
Irradiation wavelength: 395 nm
Thermal R→P:            0 s⁻¹
Thermal P→R:            6.3e-5 s⁻¹
```

For a chemical actinometer, select **Use photon flux from a chemical
actinometer** and enter **Photon flux (mol photons/s)** instead of power. The
LED file and LED-processing settings are then ignored. The irradiation
wavelength is required and represents the actual monochromatic wavelength.

AutoQY expects thermal rate constants in $\mathrm{s^{-1}}$.

For a first-order process:

```math id="gmsv6h"
k = \frac{\ln 2}{t_{1/2}}
```

---

# 17. Choose the fitting method

Expand:

```text id="llqnp3"
4 · Fit → Kinetic model
```

The available methods are:

- **Regularized concentrations**
- **Full-spectrum ODE absorbance**
- **Concentrations — legacy pure NNLS**
- **Emission — legacy**

For a new experiment, it can be useful to compare:

```text id="v2m2ar"
Regularized concentrations
vs.
Full-spectrum ODE absorbance
```

## Regularized concentrations

This is the recommended concentration-based route.

Each timestamp retains an independently adjustable reactant fraction, but the fractions are softly regularized toward an exponential envelope.

> **The exponential envelope is a regularizer. It is not the photochemical quantum-yield model.**

The resulting concentration trajectory is then fitted using the photochemical kinetic equations to obtain the quantum yields.

<details>
<summary><strong>More about the fitting methods</strong></summary>

### Regularized concentrations

The method fits all spectra together while using a conserved total concentration.

The regularizer reduces unrealistic point-to-point fluctuations.

`Concentration regularization strength` controls how strongly the fractions are encouraged to follow the smooth envelope.

The default is:

```text id="i4s1yh"
1
```

### Full-spectrum ODE absorbance

This method directly fits the complete wavelength × time absorbance matrix to the photochemical kinetic model.

It jointly fits:

- quantum yields;
- total concentration;
- initial composition;
- spectral evolution.

Optional per-spectrum baseline corrections and robust loss can help accommodate small baseline differences or isolated problematic wavelengths.

### Concentrations — legacy pure NNLS

This method independently decomposes every spectrum into reactant and product using non-negative least squares and then fits the resulting concentration trajectory.

It is simple and fast, but does not use temporal information during spectral decomposition.

### Emission — legacy

This method uses absorbance information mainly inside the active LED-emission region.

It is retained for compatibility with older AutoQY analyses.

</details>

---

# 18. Optimizer settings and expected PSS

Under:

```text id="54a7on"
Expert optimizer settings
```

the default initial values are:

```text id="jsn1p7"
Initial Φ R→P: 0.5
Initial Φ P→R: 0.5

Lower bound: 0
Upper bound: 1
```

These are numerical starting guesses rather than expected physical values.

If you have independently measured the PSS composition, enter:

```text id="qkhvda"
Expected reactant at PSS (%)
```

This provides an additional comparison between the fit and an independent experiment.

---

# 19. Compare fit methods

Click:

```text id="7lgauh"
Compare fit methods
```

For example:

| Method | Φ R→P |
|---|---:|
| Regularized concentrations | 0.25 |
| Full-spectrum ODE | 0.24 |
| Legacy NNLS | 0.26 |

shows good agreement.

If the methods differ substantially, it can be useful to inspect the spectra, reference data, and fitting assumptions.

Long method names wrap within the comparison table. On a narrow window, scroll
the table horizontally to inspect all columns. Quantum-yield values, PSS values,
and fit-status text also wrap inside their result cards rather than being cut
off.

---

# 20. Residuals and diagnostics

Residuals can be useful diagnostics, but their interpretation depends on the fitting method.

## Regularized-concentration fraction residuals

Some structure in the fraction residuals is not necessarily problematic.

The concentration recovery itself contains a soft temporal regularizer, so these residuals should not be interpreted in exactly the same way as residuals from a direct kinetic fit.

## Full-spectrum ODE absorbance residuals

The wavelength-resolved absorbance residuals provide a more direct view of how well the spectral/kinetic model reproduces the measured spectra.

<details>
<summary><strong>How to interpret the residual plots</strong></summary>

### Fraction residuals

The regularized concentration method encourages the recovered concentration points to follow a smooth exponential envelope.

Real photochemical kinetics can differ somewhat from that envelope, especially when:

- both photochemical directions are active;
- thermal reactions compete with photochemistry;
- absorption changes as composition changes.

Small structured fraction residuals can therefore occur even when the overall analysis is reasonable.

### Wavelength-resolved absorbance residuals

A coherent wavelength-dependent residual may suggest:

- differences between the reference and experimental spectra;
- an additional absorbing species;
- photodecomposition;
- spectral drift;
- concentration changes;
- baseline differences;
- limitations of the two-state model.

It is usually most useful to consider the residuals together with:

- method comparison;
- endpoint reconstruction;
- expected PSS;
- concentration evolution.

</details>

---

# 21. ε uncertainty

Expand:

```text id="6xpjgx"
5 · Uncertainty
```

If your ε reference files were produced from replicate measurements in Spectral Treatment, select:

```text id="0qv3q8"
ε range
```

You can choose:

- **Standard deviation**
- **Standard error**

For independently prepared samples, SD describes the observed variability between preparations.

SEM describes uncertainty in the estimated mean:

```math id="f67jvv"
SEM = \frac{SD}{\sqrt{n}}
```

The appropriate choice depends on which type of uncertainty you want to represent.

---

# 22. Output and run the analysis

Expand:

```text id="egh2t8"
6 · Output
```

Available outputs include:

- TXT summary;
- PNG figures;
- SVG figures;
- result JSON;
- configuration snapshot;
- detailed CSV data.

Then click:

```text id="7ow8dk"
Save JSON
```

to save the editable:

```text id="8wsv44"
analysis.json
```

and finally:

```text id="6qam4p"
Run analysis
```

A convenient reproducible workflow is:

```text id="7dx7kt"
Save JSON
    ↓
Run analysis
```

---

# 23. Inspect the results

After running the analysis, inspect:

- **Concentrations**
- **Fraction residual**
- **Preprocessing**
- **Endpoint reconstruction**
- **Absorbance residuals**

## Concentrations

For a simple R → P experiment, the reactant and product fractions will normally evolve reasonably smoothly.

## Preprocessing

Check the relationship between:

- experimental spectral evolution;
- reference ε spectra;
- LED emission.

One useful question is whether the LED emission overlaps sufficiently with the absorption spectrum.

## Endpoint reconstruction

This compares the measured endpoint with a spectrum reconstructed from the supplied reactant and product references.

Differences may arise from:

- reference-spectrum mismatch;
- baseline effects;
- decomposition;
- aggregation;
- concentration changes;
- additional species.

## Absorbance residuals

Persistent wavelength-dependent residual features can help identify spectral behaviour not fully captured by the model.

---

# 24. What should I look for in a final analysis?

A useful final check is whether the different parts of the analysis give a consistent picture.

For example:

1. Reactant ε replicate measurements are reasonably consistent.
2. Product ε has a plausible spectral shape.
3. The starting spectrum is compatible with the reference spectra.
4. Endpoint reconstruction is reasonable.
5. Concentration evolution is chemically plausible.
6. Regularized and full-spectrum ODE analyses give similar quantum yields.
7. Wavelength-resolved residuals do not contain large unexplained features.
8. Calculated PSS is reasonably consistent with an independent PSS measurement when available.
9. Small changes in preprocessing or regularization do not strongly change the result.

The different diagnostics are intended to show how sensitive the result is to the assumptions used in the analysis.

---

<details>
<summary><strong>Worked example: Example 4 — 395 nm with ε uncertainty</strong></summary>

AutoQY contains:

```text id="2cqivz"
ExampleData/Example-4_395nm-EpsilonError
```

The example uses three independently prepared reactant solutions:

```text id="q3tgqj"
6.82 × 10⁻⁵ M
6.87 × 10⁻⁵ M
6.80 × 10⁻⁵ M
```

with:

```text id="khcjwt"
Path length:             1 cm
Reactant at PSS:         23%
NMR composition error:   1%
Sample volume:           1995 µL
Power:                   1.46 ± 0.03 mW
Irradiation wavelength:  395 nm
Thermal P→R rate:        6.3 × 10⁻⁵ s⁻¹
```

The expected result is approximately:

```math id="6p5tqq"
\Phi_{R\rightarrow P} = 25.0 \pm 1.2\%
```

```math id="kmt79q"
\Phi_{P\rightarrow R} = 21 \pm 3\%
```

with an extrapolated PSS of approximately:

```math id="830jpu"
23.3\%\,R
```

Small differences can occur if preprocessing or uncertainty settings are changed.

</details>

---

# Quick workflows

<details>
<summary><strong>Make a spectral figure</strong></summary>

```text id="8qqugo"
Load spectra
        ↓
Choose wavelength range
        ↓
Baseline if useful
        ↓
SavGol if needed
        ↓
Minimal colors
        ↓
Hide all legend entries
        ↓
Enable only important traces
        ↓
Rename them
        ↓
Origin-style export
        ↓
Save SVG or PNG
```

</details>

<details>
<summary><strong>Measure a simple lifetime</strong></summary>

```text id="e61c0k"
Load time-ordered spectra
        ↓
Choose wavelength slice
        ↓
Set Seconds per timestamp
        ↓
Enable Fit exponential decay
        ↓
Read τ ± error
        ↓
Check measurement duration
```

</details>

<details>
<summary><strong>Calculate ε</strong></summary>

```text id="em9wn3"
Load independent replicate spectra
        ↓
Choose wavelength range
        ↓
Baseline if needed
        ↓
SavGol if useful
        ↓
SVD usually OFF
        ↓
Enter concentration
        ↓
Enter path length
        ↓
Inspect individual ε curves
        ↓
Save processed CSV
```

</details>

<details>
<summary><strong>Determine a quantum yield</strong></summary>

```text id="gsjo3s"
Prepare reactant ε
        ↓
Prepare or reconstruct product ε
        ↓
Open Analysis GUI
        ↓
Load measurement spectra
        ↓
Load εR and εP
        ↓
Load LED spectrum
        ↓
Load timestamps
        ↓
Enter volume, path length, power, thermal rates
        ↓
Run Regularized concentrations
        ↓
Compare with Full-spectrum ODE
        ↓
Inspect diagnostics
        ↓
Apply ε uncertainty if appropriate
        ↓
Save JSON
        ↓
Run final analysis
```

</details>

---

# A few practical points

- For independent ε measurements, keeping the replicate variability visible is usually preferable to reducing it with SVD.
- The exponential in **Regularized concentrations** is a soft regularizer, not the quantum-yield kinetic model.
- Small structured fraction residuals from the regularized method can occur and should be interpreted together with the other diagnostics.
- Wavelength-resolved residuals from the full-spectrum ODE method can be useful for identifying spectral features not captured by the model.
- Comparing more than one fitting approach is a useful way to assess how sensitive the result is to the analysis method.
- For lifetime measurements, a longer time window generally gives a better estimate of both $\tau$ and the plateau.
- For spectral figures, intermediate traces can remain visible even if only the first and last spectra are shown in the legend.
- For files above 60 spectra, the displayed intermediate traces are an evenly
  spaced preview; processing and CSV export still include the full series.
