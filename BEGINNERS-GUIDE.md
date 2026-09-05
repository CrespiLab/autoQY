# AutoQY Spectral Treatment and Analysis GUI

## Beginner tutorial: from UV–Vis spectra to publication-ready plots, lifetimes, and quantum yields

AutoQY contains two main graphical interfaces:

- **Spectral Treatment** — prepares UV–Vis spectra, calculates molar absorptivity, inspects kinetic traces, fits simple exponential lifetimes, and exports publication-ready figures.
- **Analysis GUI** — combines spectral evolution, molar absorptivities, irradiation times, LED emission, optical power, and experimental parameters to calculate photoisomerization quantum yields.

For a standard quantum-yield experiment, the workflow is:

```text
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

However, **Spectral Treatment is also useful as a standalone UV–Vis tool**.

You can use it to:

- baseline-correct spectra;
- smooth spectra;
- compare spectral series;
- calculate molar absorptivity;
- estimate uncertainty from replicate ε measurements;
- reconstruct a photoproduct spectrum from an NMR-characterized PSS;
- monitor absorbance at one wavelength during a kinetic experiment;
- fit a simple exponential lifetime;
- prepare clean spectral figures;
- hide unnecessary legend entries;
- rename displayed traces;
- export publication-ready PNG or SVG figures.

You normally do **not** need to edit Python code or JSON manually.

---

# 1. What do I need for a quantum-yield experiment?

For a two-state photochemical reaction:

$$
R \rightleftharpoons P
$$

you ultimately need:

| Input | Meaning |
|---|---|
| Measurement spectra | UV–Vis spectra recorded during irradiation |
| Reactant ε | Molar absorptivity spectrum of R |
| Product ε | Molar absorptivity spectrum of P |
| LED spectrum | Measured emission spectrum of the irradiation LED |
| Irradiation timestamps | Irradiation time corresponding to every spectrum |

You also need:

- sample volume;
- optical path length;
- irradiation power;
- uncertainty in irradiation power;
- thermal R → P rate, if relevant;
- thermal P → R rate, if relevant.

The irradiation wavelength entered in AutoQY is primarily metadata and a consistency check.

The photon-flux calculation uses the **complete processed LED emission spectrum**, rather than treating the LED as perfectly monochromatic.

---

# 2. Starting AutoQY

After installation, open:

```text
AutoQY Analysis
```

Inside the Analysis GUI, expand:

```text
1 · Project
```

and click:

```text
Open Spectral Treatment
```

Spectral Treatment opens in a separate window.

Use Spectral Treatment whenever you need to:

- prepare spectra;
- calculate ε;
- inspect a spectral series;
- perform simple kinetic analysis;
- prepare a figure.

---

# 3. Spectral Treatment: loading spectra

Open:

```text
1 · Data → Spectral data
```

You can:

- drag and drop spectra;
- choose individual files;
- use **Open files from folder**.

Supported formats include:

- SpectraGryph `.dat`;
- Avantes `.Abs8`;
- TSV;
- CSV.

Multiple spectra can be loaded together.

## Check the spectrum order

The order matters.

For molar-absorptivity measurements, the concentrations entered later must correspond to the correct spectra.

For kinetic measurements, the spectra should normally be in chronological order.

Expand:

```text
Loaded spectra: order, legend, removal
```

Here you can:

- move individual spectra up or down;
- remove spectra;
- select which traces appear in the legend;
- rename displayed legend labels.

Changing the displayed name **does not rename the original file**.

---

# 4. Select the useful wavelength range

Open:

```text
2 · Range → Wavelengths
```

Enter the spectral range that is actually useful.

For example:

```text
250–700 nm
```

Do not automatically use the complete detector range.

If your spectrometer gives unreliable values below 230 nm or above 800 nm, there is no benefit in including those regions.

The selected range is used for:

- displayed spectra;
- molar-absorptivity calculations;
- uncertainty calculations;
- exported data.

---

# 5. Baseline correction

Inside **2 · Range**, expand:

```text
Preprocess spectra
```

Enable:

```text
Baseline
```

Then choose a wavelength interval in which your compound should have essentially no absorbance.

For example:

```text
600–650 nm
```

but only if your compound genuinely does not absorb there.

If a supposedly zero-absorbance region sits around $A = 0.003$, baseline correction can remove this offset.

## Important

Do not select the baseline region blindly.

A poor baseline interval can distort the entire spectrum.

The software cannot know whether your molecule really absorbs in the selected region.

---

# 6. Savitzky–Golay smoothing

Under preprocessing, select:

```text
SavGol
```

A reasonable starting point for ordinary UV–Vis spectra is:

```text
Window: 5 nm
Polynomial order: 3
```

The program converts the selected window in nanometres into an appropriate odd number of detector points.

Use smoothing conservatively.

If the molecular absorption band changes shape substantially rather than merely becoming less noisy, the smoothing is probably too aggressive.

---

# 7. What about SVD?

Spectral Treatment also provides:

```text
SVD
```

For independently prepared molar-absorptivity measurements:

> **Leave SVD OFF.**

If you independently prepare several solutions, their differences contain real information about experimental reproducibility.

Those differences are exactly what you want when estimating uncertainty in ε.

Applying SVD across those spectra can suppress genuine replicate-to-replicate variation.

## Simple rule

```text
Independent replicate solutions
        ↓
SVD OFF
```

```text
Ordered spectral time series
        ↓
SVD may be useful
```

Even for a time series, inspect the untreated data before deciding that SVD is justified.

---

# 8. Calculating molar absorptivity

Open:

```text
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

AutoQY applies the Beer–Lambert law independently to every spectrum:

```math
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

$$
SEM=\frac{SD}{\sqrt{n}}
$$

---

# 9. Inspect the ε spectra before saving

Do not immediately save the calculated result.

Look at the individual ε curves.

Independent preparations should normally produce similar spectra.

Small random differences are expected.

Large systematic differences should be investigated.

For example:

```text
Sample 1: εmax = 18,000 M⁻¹ cm⁻¹
Sample 2: εmax = 18,300 M⁻¹ cm⁻¹
Sample 3: εmax = 25,000 M⁻¹ cm⁻¹
```

should not simply be treated as a large error bar.

Possible causes include:

- incorrect concentration;
- dilution error;
- incorrect baseline;
- dirty or mismatched cuvette;
- aggregation;
- decomposition;
- instrumental problems;
- incorrect sample.

AutoQY can calculate an SD from bad experiments. That does not make the measurements good.

---

# 10. Exporting molar absorptivity

Open:

```text
4 · Output → Export processed dataset
```

Choose a destination and click:

```text
Save processed CSV
```

When concentrations and path lengths have been supplied, the exported file contains:

- wavelength;
- processed absorbance;
- concentration;
- path length;
- individual ε spectra;
- mean ε;
- SD;
- SEM;
- lower non-negative ε bound;
- upper ε bound.

For example:

```text
reactant_absorptivity.csv
```

---

# 11. Spectral Treatment as a simple kinetic-analysis tool

Spectral Treatment can also be used independently of the quantum-yield workflow.

Suppose you record a UV–Vis spectrum repeatedly while a reaction occurs:

```text
Spectrum 0
Spectrum 1
Spectrum 2
Spectrum 3
...
Spectrum 20
```

You may simply want to ask:

> How quickly does the absorbance at 450 nm change?

You do not need the complete quantum-yield analysis for this.

---

# 12. Extracting a kinetic trace at one wavelength

Load the spectra in chronological order.

Use the wavelength-slice controls and select the wavelength of interest.

For example:

```text
450 nm
```

Spectral Treatment extracts:

$$
A(450\ \mathrm{nm},t)
$$

from every spectrum and plots absorbance against the spectrum coordinate or time.

If the selected wavelength lies between two detector points, the value is interpolated.

This is useful for:

- thermal photoswitch relaxation;
- reaction kinetics;
- degradation experiments;
- preliminary kinetic screening.

---

# 13. Setting the time axis

The wavelength-slice section contains:

```text
Seconds per timestamp
```

If the stored coordinates are already seconds, leave:

```text
1
```

For example:

```text
0, 30, 60, 90, 120...
```

already represents seconds.

If instead the spectra are simply numbered:

```text
0, 1, 2, 3, 4...
```

and one spectrum was recorded every 30 seconds, enter:

```text
30
```

The program then interprets the series as:

```text
0, 30, 60, 90, 120... s
```

---

# 14. Fitting a simple exponential lifetime

Enable:

```text
Fit exponential decay
```

AutoQY fits:

$$
A(t)=A_{\infty}+\Delta A\,e^{-t/\tau}
$$

where:

- $A_{\infty}$ is the final offset;
- $\Delta A$ is the amplitude;
- $\tau$ is the exponential lifetime.

The fitted curve is shown on the plot.

AutoQY reports:

$$
\tau \pm SE(\tau)
$$

where the uncertainty is the one-standard-error uncertainty from the fit.

## Lifetime versus half-life

The reported lifetime is:

$$
\tau=\frac{1}{k}
$$

For a first-order process:

$$
t_{1/2}=\tau\ln 2
$$

or equivalently:

$$
k=\frac{\ln 2}{t_{1/2}}
$$

Therefore:

> **Lifetime and half-life are not the same number.**

---

# 15. Is the kinetic measurement long enough?

AutoQY checks whether the measured time span extends beyond approximately one fitted lifetime.

If:

$$
t_{\mathrm{experiment}} < \tau
$$

a warning is shown.

This matters because fitting only the beginning of an exponential decay gives poor information about:

- the final plateau;
- the amplitude;
- the lifetime.

Ideally, measure substantially beyond one lifetime and preferably far enough to characterize the final plateau.

## When is this fit appropriate?

The simple exponential fit is useful for:

- first-order thermal isomerization;
- simple unimolecular relaxation;
- degradation screening;
- approximate kinetic measurements.

Do not force a single exponential onto data showing:

- clear biexponential behaviour;
- induction periods;
- multiple sequential reactions;
- oscillations;
- strong mechanistic complexity.

A converged exponential fit does not prove first-order kinetics.

---

# 16. Spectral Treatment as a publication-figure tool

Spectral Treatment can also prepare clean figures directly.

This is particularly useful for irradiation experiments containing many spectra.

Suppose you have 30 spectra showing an evolution from an initial state to a PSS.

You may want all 30 spectra visible, but you almost certainly do not want a 30-entry legend.

AutoQY separates:

```text
Which traces are plotted
```

from:

```text
Which traces appear in the legend
```

---

# 17. Cleaning up the legend

Expand:

```text
Loaded spectra: order, legend, removal
```

For a large spectral series, click:

```text
Hide all
```

This hides all trace names from the legend.

It does **not** remove the spectra from the figure.

Then selectively enable only the important entries.

For example:

```text
First spectrum
Last spectrum
```

The complete spectral evolution remains visible, while the legend contains only chemically meaningful information.

---

# 18. Renaming legend entries

Each loaded spectrum has an editable display name.

Changing this name alters the legend without changing the source file.

For example, the actual files might be:

```text
Compound1_0000s.Abs8
Compound1_0010s.Abs8
Compound1_0020s.Abs8
...
Compound1_0300s.Abs8
```

For a manuscript figure you may want only:

```text
E
PSS
```

or:

```text
Initial
PSS365
```

## Recommended workflow

For a 30-spectrum irradiation series:

1. Load all spectra.
2. Put them in chronological order.
3. Click **Hide all**.
4. Enable the legend only for spectrum 1 and spectrum 30.
5. Rename spectrum 1 to `E`.
6. Rename spectrum 30 to `PSS`.

You now have all 30 curves but only two meaningful legend entries.

---

# 19. Minimal colors

Enable:

```text
Minimal colors
```

This highlights the initial and final spectra while showing intermediate spectra in a restrained neutral style.

The purpose is to create an immediate visual hierarchy:

```text
Start → intermediate evolution → final state
```

without producing a rainbow of arbitrary colours.

This is particularly useful for photochemical irradiation series.

---

# 20. Exporting publication-ready figures

Spectral Treatment can save the plot as:

- PNG
- SVG

SVG is particularly useful for manuscripts because it is vector based and can subsequently be edited if necessary.

Expand:

```text
Image export options
```

Options include controls for:

- title;
- legend;
- grid;
- y-axis starting at zero;
- Origin-style export.

## Origin-style export

The Origin-style option produces a clean publication-oriented figure rather than simply saving a screenshot of the browser interface.

For a manuscript figure:

1. Load the spectra.
2. Select the useful wavelength range.
3. Apply preprocessing only if scientifically justified.
4. Enable **Minimal colors**.
5. Click **Hide all** for the legend.
6. Re-enable only the first and final spectra.
7. Rename them with meaningful chemical labels.
8. Remove the title if unnecessary.
9. Keep the legend if required.
10. Enable **Origin-style export**.
11. Save as SVG or PNG.

---

# 21. Obtaining the product molar absorptivity

There are two main cases.

## Case A — Pure product is available

Prepare pure P at known concentration and repeat the same procedure:

1. Load replicate product spectra.
2. Select wavelength range.
3. Baseline if appropriate.
4. Smooth only if necessary.
5. Keep SVD off for independent preparations.
6. Enter concentrations.
7. Enter path lengths.
8. Inspect the individual ε curves.
9. Export:

```text
product_absorptivity.csv
```

---

# 22. Case B — The product cannot be isolated

This situation is common for photoswitches.

Suppose irradiation produces a PSS containing:

```text
77% P
23% R
```

and NMR independently determines that 23% reactant remains.

AutoQY can reconstruct the product spectrum.

---

# 23. NMR-guided PSS subtraction

Expand:

```text
5 · Optional → NMR-guided PSS subtraction
```

Load a UV–Vis dataset in which:

```text
First spectrum = pure reactant
Last spectrum  = final PSS
```

Enter:

```text
Reactant in final PSS (%)
```

For example:

```text
23
```

If $x$ is the reactant fraction at the PSS, AutoQY reconstructs the product from:

$$
P=
\frac{PSS-xR}{1-x}
$$

For $x=0.23$:

$$
P=
\frac{PSS-0.23R}{0.77}
$$

---

# 24. Entering the NMR uncertainty

Enter:

```text
NMR error (%)
```

For example:

```text
1
```

The uncertainty in the experimentally measured PSS composition is propagated into the reconstructed product ε spectrum.

The resulting uncertainty is generally:

- wavelength dependent;
- asymmetric.

The uncertainty of the reactant ε measurement is also included.

---

# 25. Negative reconstructed ε values

Subtraction can produce negative ε values.

A true molar absorptivity cannot physically be negative.

Small negative values around a nominally zero baseline can arise from experimental noise.

A substantial negative spectral band is a warning.

Possible causes include:

- incorrect PSS composition;
- poor baseline correction;
- normalization errors;
- decomposition;
- more than two species being present.

Do not treat a strongly negative reconstructed spectrum as a successful two-component analysis.

---

# 26. Exporting the reactant and product references

Save the two reference files as, for example:

```text
reactant_absorptivity.csv
product_absorptivity.csv
```

You are now ready for the quantum-yield analysis.

---

# 27. Return to AutoQY Analysis

The Analysis GUI contains:

```text
1 · Project
2 · Identity and data
3 · Experiment
4 · Fit
5 · Uncertainty
6 · Output
7 · Analyze
```

Work through them in that order.

---

# 28. Project

Expand:

```text
1 · Project → JSON and tools
```

Choose the folder containing your experiment.

For example:

```text
MyExperiment/
├── measurement_spectra.csv
├── reactant_absorptivity.csv
├── product_absorptivity.csv
├── led_emission.csv
└── timestamps.csv
```

Think of `analysis.json` as the recipe for the analysis.

It stores:

- input files;
- experimental parameters;
- fitting method;
- preprocessing;
- uncertainty treatment;
- output settings.

---

# 29. Identity and data

Expand:

```text
2 · Identity and data → Experiment files
```

Choose a meaningful:

```text
Analysis ID
```

for example:

```text
compound1_395nm
```

and an output stem such as:

```text
compound1_395nm_results
```

You can replace generic species names such as `reactant` and `product` with chemically useful labels such as `E` and `Z`.

---

# 30. Select the input files

Select:

1. **Measurement spectra**
2. **Reactant molar absorptivity**
3. **Product molar absorptivity**
4. **LED emission**
5. **Irradiation timestamps**

For newly prepared data, use:

```text
Generic CSV (recommended)
```

where possible.

A typical measurement file looks conceptually like:

```csv
Wavelength_nm,Spectrum_1,Spectrum_2,Spectrum_3
250,...
251,...
252,...
```

---

# 31. LED processing

Inside the LED section expand:

```text
Processing
```

A reasonable starting point is:

```text
Wavelength start:              250 nm
Wavelength end:                800 nm
Savitzky–Golay window:         12 points
Polynomial order:              3
Baseline correction:           ON
Baseline exclusion multiplier: 10
```

These settings apply only to the **LED emission spectrum**.

They do not preprocess the experimental absorbance spectra.

---

# 32. Experimental parameters

Expand:

```text
3 · Experiment → Physical parameters
```

Enter the actual experimental values.

For example:

```text
Sample volume:          1995 µL
Path length:            1 cm
Power:                  1.46 mW
Power error:            0.03 mW
Irradiation wavelength: 395 nm
Thermal R→P:            0 s⁻¹
Thermal P→R:            6.3e-5 s⁻¹
```

---

# 33. Thermal rate constants

AutoQY expects thermal **rate constants**, in $\mathrm{s^{-1}}$, not half-lives.

For a first-order reaction:

$$
k=\frac{\ln 2}{t_{1/2}}
$$

If thermal conversion is negligible on the timescale of the experiment, use zero.

---

# 34. Choosing the fitting method

Expand:

```text
4 · Fit → Kinetic model
```

The fitting methods are:

- **Regularized concentrations**
- **Full-spectrum ODE absorbance**
- **Concentrations — legacy pure NNLS**
- **Emission — legacy**

For a new experiment, the most useful starting comparison is:

```text
Regularized concentrations
vs.
Full-spectrum ODE absorbance
```

---

# 35. Regularized concentrations

This is the recommended concentration-based route.

The method fits all spectra together while enforcing a conserved total concentration.

Each timestamp retains an independently adjustable reactant fraction.

However, those fractions are **softly regularized toward an exponential envelope** with:

- a free starting fraction;
- a free plateau.

> **The exponential envelope is a regularizer. It is not the photochemical quantum-yield model.**

The regularizer stabilizes the spectral decomposition and discourages noisy, physically implausible jumps between consecutive spectra.

The resulting concentration trajectory is subsequently fitted using the full photochemical kinetic equations to determine the quantum yields.

The:

```text
Concentration regularization strength
```

controls how strongly the fractions are encouraged to follow the exponential envelope.

The default is:

```text
1
```

---

# 36. Full-spectrum ODE absorbance

The:

```text
Full-spectrum ODE absorbance
```

method directly fits the complete wavelength × time absorbance matrix to the photochemical kinetic model.

It jointly fits:

- quantum yields;
- total concentration;
- initial composition;
- spectral evolution.

Optional per-spectrum baseline corrections and robust loss can handle small systematic offsets and isolated bad wavelengths.

This method is slower but provides a conceptually different route to the answer.

Agreement with the regularized-concentration method is therefore particularly useful.

---

# 37. Legacy pure-NNLS concentrations

The:

```text
Concentrations (legacy pure NNLS)
```

method independently decomposes every spectrum into reactant and product using non-negative least squares.

The resulting concentrations are then fitted to the photochemical kinetic model.

Advantages:

- simple;
- fast;
- transparent.

Disadvantages:

- no temporal information is used during spectral decomposition;
- spectral mismatch can pin points to 0% or 100%;
- noisy data can produce irregular concentration traces.

It remains useful as an independent comparison.

---

# 38. Legacy emission method

The:

```text
Emission
```

method is retained for compatibility with older AutoQY analyses.

It uses the absorbance behaviour primarily inside the active LED-emission region.

It can become poorly conditioned when reactant and product absorb similarly in the irradiated wavelength region.

For a new experiment with good full spectral information, it would normally not be the first method to rely on.

---

# 39. Optimizer settings

Expand:

```text
Expert optimizer settings
```

The defaults are:

```text
Initial Φ R→P: 0.5
Initial Φ P→R: 0.5

Lower bound: 0
Upper bound: 1
```

Normally leave these alone initially.

The values of 0.5 are only numerical starting guesses.

---

# 40. Expected PSS

If you have independently measured the PSS composition, enter:

```text
Expected reactant at PSS (%)
```

For example:

```text
23
```

If NMR gives 23% R while the fitted model predicts 60% R, the discrepancy should be investigated.

---

# 41. Compare fit methods

Before accepting the final result, click:

```text
Compare fit methods
```

For example:

| Method | Φ R→P |
|---|---:|
| Regularized concentrations | 0.25 |
| Full-spectrum ODE | 0.24 |
| Legacy NNLS | 0.26 |

would be reassuring.

By contrast:

| Method | Φ R→P |
|---|---:|
| Regularized concentrations | 0.25 |
| Full-spectrum ODE | 0.08 |
| Legacy NNLS | 0.62 |

means the result is strongly method dependent.

Do not simply choose the number you prefer.

---

# 42. How should residuals be interpreted?

Residuals are useful, but they must be interpreted according to the fitting method.

## Regularized-concentration fraction residuals

The regularized concentration method deliberately encourages the independently fitted concentration points to follow a smooth exponential envelope.

Therefore:

> **A structured fraction residual in the regularized-concentration analysis is not automatically evidence that the photochemical mechanism is wrong.**

The spectral decomposition itself has already been influenced by a temporal regularizer.

Real photochemical kinetics do not necessarily have to be a mathematically perfect single exponential, particularly when:

- both photochemical directions are active;
- thermal reactions compete with photochemistry;
- excitation conditions effectively change as composition changes.

Do not demand perfectly random fraction residuals from the regularized route.

Instead, use them to evaluate:

- how strongly the recovered trajectory differs from the final photochemical ODE fit;
- whether deviations are small or large;
- whether the result changes substantially when regularization strength is changed.

## Full-spectrum ODE absorbance residuals

The wavelength-resolved absorbance residuals from the full-spectrum ODE analysis are a stronger test of whether the assumed spectral/kinetic model reproduces the experimental data.

A coherent wavelength-dependent residual can indicate:

- incorrect reference spectra;
- an unmodelled third species;
- photodecomposition;
- spectral drift;
- concentration changes;
- baseline problems;
- limitations of the two-state model.

## Best practice

Do not judge a result from one residual plot alone.

Look for consistency between:

- Regularized concentrations
- Full-spectrum ODE absorbance
- Legacy concentration fit
- Endpoint reconstruction
- Expected PSS
- Wavelength-resolved residuals

A residual pattern appearing only in the regularized-concentration fraction plot is less concerning than a spectral feature that remains unexplained in the full-spectrum ODE residual map.

---

# 43. ε uncertainty

Expand:

```text
5 · Uncertainty
```

If your ε reference files were produced from replicate measurements in Spectral Treatment, select:

```text
ε range
```

You can then choose:

- **Standard deviation**
- **Standard error**

For independently prepared samples, SD describes the observed variability between preparations.

SEM describes uncertainty in the estimated mean:

$$
SEM=\frac{SD}{\sqrt{n}}
$$

and is therefore smaller.

Do not choose SEM simply because it produces smaller quantum-yield error bars.

---

# 44. Output settings

Expand:

```text
6 · Output
```

Choose a results folder such as:

```text
results
```

For a complete analysis, it is useful to enable:

- TXT summary;
- PNG figures;
- SVG figures;
- result JSON;
- configuration snapshot;
- detailed CSV data.

Leave overwrite disabled unless replacing an existing analysis is intentional.

---

# 45. Save the analysis configuration

Click:

```text
Save JSON
```

This saves:

```text
analysis.json
```

> **Run analysis does not automatically save the editable `analysis.json`.**

For a final reproducible analysis:

```text
Save JSON
    ↓
Run analysis
```

---

# 46. Run the analysis

Click:

```text
Run analysis
```

The GUI validates the configuration and performs the selected fit.

The result cards report:

```text
R → P quantum yield
P → R quantum yield
```

Do not stop after reading these two numbers.

---

# 47. Preprocessing and fit diagnostics

Open:

```text
Preprocessing and fit diagnostics
```

AutoQY uses approximately:

```text
Green  = no automatic threshold exceeded
Amber  = inspect carefully
Red    = invalid or strongly unstable
```

A green result does not prove that the experiment is scientifically correct.

Diagnostics are aids, not replacements for chemical judgement.

---

# 48. Concentration plot

Inspect:

```text
Concentrations
```

For a simple R → P experiment you would normally expect reasonably smooth behaviour:

```text
R decreases
P increases
```

A strongly irregular trajectory deserves investigation.

Remember that the concentration methods differ in how much temporal information they use.

Do not expect legacy NNLS and the regularized method to behave identically.

---

# 49. Fraction residuals

Inspect:

```text
Fraction residual
```

Use this plot to compare the recovered concentration trajectory with the subsequent photochemical kinetic fit.

Large residuals are worth investigating.

Small structured residuals, particularly for the regularized-concentration route, are not automatically evidence of a wrong chemical mechanism because the concentration recovery itself is softly biased toward an exponential envelope.

Interpret this plot together with:

- method comparison;
- full-spectrum residuals;
- PSS agreement;
- endpoint reconstruction.

---

# 50. Preprocessing plot

Inspect:

```text
Preprocessing
```

This allows you to view together:

- experimental spectral evolution;
- reference ε spectra;
- LED emission.

Ask:

> Does the LED actually overlap a region where the relevant species absorbs?

A mathematically converged calculation cannot rescue an experimentally inappropriate irradiation wavelength.

---

# 51. Endpoint reconstruction

Inspect:

```text
Endpoint reconstruction
```

This asks whether the measured endpoint can be represented using the supplied reactant and product spectra.

Poor reconstruction may indicate:

- incorrect ε spectra;
- baseline problems;
- decomposition;
- aggregation;
- changing sample concentration;
- formation of a third species.

---

# 52. Absorbance residuals

Inspect:

```text
Absorbance residuals
```

This wavelength-resolved residual map is particularly important for the full-spectrum ODE method.

Ideally, there should be no strong coherent spectral feature that the model systematically fails to reproduce.

For example, a persistent residual band around 520 nm could indicate an absorbing component missing from the assumed R + P model.

Possible explanations include:

- photoproduct formation;
- intermediate formation;
- spectral drift;
- incorrect molar absorptivity;
- decomposition;
- failure of the two-state model.

A coherent wavelength-resolved residual is generally more chemically informative than a small systematic deviation in the regularized concentration trace.

---

# 53. What makes a quantum-yield result convincing?

A defensible result should ideally satisfy most of the following:

1. Reactant ε replicate measurements agree reasonably well.
2. Product ε is physically plausible.
3. The starting experimental spectrum is compatible with the reference spectra.
4. Endpoint reconstruction is reasonable.
5. Concentration evolution is chemically sensible.
6. Regularized and full-spectrum ODE analyses give compatible quantum yields.
7. Legacy NNLS does not reveal a severe inconsistency.
8. Wavelength-resolved ODE residuals do not show a major unexplained spectral component.
9. Calculated PSS agrees with independently measured PSS when available.
10. Quantum yields are not artificially stuck on optimizer bounds.
11. Reasonable changes in preprocessing do not radically change Φ.
12. Reasonable changes in regularization strength do not qualitatively change the answer.

The objective is not to make every diagnostic plot perfectly featureless.

The objective is to show that the derived quantum yield is **robust to reasonable ways of analysing the same physical experiment**.

---

# 54. Worked example supplied with AutoQY

AutoQY contains:

```text
ExampleData/Example-4_395nm-EpsilonError
```

The example uses three independently prepared reactant solutions:

```text
6.82 × 10⁻⁵ M
6.87 × 10⁻⁵ M
6.80 × 10⁻⁵ M
```

with:

```text
Path length:             1 cm
Reactant at PSS:         23%
NMR composition error:   1%
Sample volume:           1995 µL
Power:                   1.46 ± 0.03 mW
Irradiation wavelength:  395 nm
Thermal P→R rate:        6.3 × 10⁻⁵ s⁻¹
```

The expected result is approximately:

```math
\Phi_{R\rightarrow P}=25.0\pm1.2\%
```

```math
\Phi_{P\rightarrow R}=21\pm3\%
```

with an extrapolated PSS of approximately:

```math
23.3\%\ R
```

Small differences can occur if preprocessing or uncertainty settings are changed.

---

# 55. Quick workflow: make a publication-ready UV–Vis figure

```text
Load spectra
        ↓
Choose wavelength range
        ↓
Baseline if scientifically justified
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

A typical final legend may contain only:

```text
E
PSS
```

while all intermediate spectra remain visible.

---

# 56. Quick workflow: measure a simple lifetime

```text
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
        ↓
Inspect whether a single exponential is appropriate
```

Remember:

$$
t_{1/2}=\tau\ln 2
$$

---

# 57. Quick workflow: calculate ε

```text
Load independent replicate spectra
        ↓
Choose wavelength range
        ↓
Baseline
        ↓
SavGol only if necessary
        ↓
SVD OFF
        ↓
Enter concentration
        ↓
Enter path length
        ↓
Inspect individual ε curves
        ↓
Save processed CSV
```

---

# 58. Quick workflow: determine a quantum yield

```text
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
Inspect spectral diagnostics
        ↓
Apply ε uncertainty if appropriate
        ↓
Save JSON
        ↓
Run final analysis
        ↓
Report Φ only after checking robustness
```

---

# 59. Rules worth remembering

## Rule 1

> A numerical fit does not prove that the chemical model is correct.

## Rule 2

> Do not use SVD to make independent ε replicates artificially agree.

## Rule 3

> The exponential in **Regularized concentrations** is a soft regularizer, not the quantum-yield kinetic model.

## Rule 4

> A structured fraction residual from the regularized method alone does not prove mechanistic failure.

## Rule 5

> Coherent wavelength-resolved residuals in the full-spectrum ODE analysis deserve more serious investigation.

## Rule 6

> Agreement between independent analysis approaches is stronger evidence than one apparently perfect fit.

## Rule 7

> For lifetime measurements, make sure the experiment spans enough of the decay to constrain both $\tau$ and the final plateau.

## Rule 8

> A 30-spectrum experiment does not need a 30-entry legend. Keep the intermediate curves, hide their legend entries, and label only the chemically meaningful spectra.

---

# 60. In one sentence

**Use Spectral Treatment to make the spectral data physically interpretable, use several analysis routes to test whether the quantum yield is robust, and never let a good-looking numerical fit substitute for inspection of the actual UV–Vis data.**
