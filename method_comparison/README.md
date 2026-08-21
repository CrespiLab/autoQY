# AutoQY method comparison

This directory compares four fitting methods on the three bundled azobenzene
examples. Every PNG has a matching SVG with the same filename stem.

## Methods

- `concentrations`: legacy pure-NNLS decomposition at every timestamp,
  followed by the photochemical ODE fit.
- `emission`: legacy direct-absorbance method restricted to the active LED band
  and initialized as pure reactant.
- `regularized_concentrations`: conserved total concentration with independent
  timestamp fractions softly regularized toward an exponential envelope.
- `ode_absorbance`: direct full-spectrum photochemical ODE fit with free initial
  composition and total concentration, robust residuals, and per-spectrum
  offset/slope baseline corrections.

## Metrics

`Concentration difference` is the concentration-trajectory RMSE relative to
the legacy `concentrations` result, divided by that method's initial total
concentration. It measures agreement with the historical reference trajectory; it does not
establish which trajectory is true.

`Absorbance RMSE` evaluates every fitted trajectory over the same complete
spectral range. It is descriptive rather than a reduced chi-squared statistic
because the bundled data do not provide pointwise absorbance uncertainties.
The `ode_absorbance` model also has baseline nuisance parameters, so its lower
RMSE should not be interpreted as a complexity-penalized model comparison.

## 340 nm, raw LED

| Method | Phi R->P (%) | Phi P->R (%) | Concentration difference | Absorbance RMSE | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| `concentrations` | 19.76 ± 0.12 | 26.16 ± 0.21 | reference | 0.003112 | 1.000 | 0.0386 |
| `emission` | 17.64 ± 0.10 | 42.85 ± 0.34 | 1.31% | 0.003625 | 1.000 | 0.0378 |
| `regularized_concentrations` | 20.27 ± 0.19 | 28.55 ± 0.78 | 0.93% | 0.003150 | 0.997 | 0.0409 |
| `ode_absorbance` | 19.36 ± 0.11 | 24.26 ± 0.17 | 0.44% | 0.002131 | 0.981 | 0.0366 |

All concentration trajectories agree within 1.4%, and the forward quantum
yields cluster around 18-20%. The reverse quantum yield does not: `emission`
returns 42.9%, whereas the other methods return 24.3-28.5%. The experiment ends
near 96% product, so the reverse contribution is inferred from a comparatively
small reactant population and is more sensitive to the LED spectrum and model.
The between-method spread is much larger than the formal within-method errors.

Files: `340nm_raw_led__<method>.png` and `.svg`.

## 340 nm, pre-baselined LED

| Method | Phi R->P (%) | Phi P->R (%) | Concentration difference | Absorbance RMSE | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| `concentrations` | 17.68 ± 0.11 | 42.36 ± 0.46 | reference | 0.002573 | 1.000 | 0.0386 |
| `emission` | 17.23 ± 0.10 | 41.18 ± 0.30 | 0.33% | 0.002660 | 1.000 | 0.0368 |
| `regularized_concentrations` | 17.94 ± 0.17 | 44.42 ± 1.27 | 0.76% | 0.002694 | 0.998 | 0.0398 |
| `ode_absorbance` | 17.19 ± 0.10 | 38.62 ± 0.26 | 0.48% | 0.001725 | 0.981 | 0.0363 |

The four concentration trajectories agree within 0.8%. Forward yields agree
within 0.8 percentage points, while reverse yields span 38.6-44.4%. This is
substantially better agreement than for the raw-LED configuration. Because the
two example directories also use different processed input files, the change
cannot be attributed exclusively to LED baseline treatment, but it confirms
that this part of the workflow materially affects the inferred reverse yield.

Files: `340nm_baselined_led__<method>.png` and `.svg`.

## 455 nm, first 40 spectra

| Method | Phi R->P (%) | Phi P->R (%) | Concentration difference | Absorbance RMSE | Initial R fraction | Final R fraction |
|---|---:|---:|---:|---:|---:|---:|
| `concentrations` | 36.77 ± 0.44 | 49.92 ± 0.63 | reference | 0.009044 | 0.964 | 0.7771 |
| `emission` | 25.79 ± 0.48 | 51.40 ± 1.09 | 9.16% | 0.049952 | 1.000 | 0.8370 |
| `regularized_concentrations` | 37.81 ± 0.38 | 50.08 ± 0.50 | 1.53% | 0.004609 | 0.983 | 0.7729 |
| `ode_absorbance` | 35.14 ± 0.50 | 50.91 ± 0.76 | 2.43% | 0.004201 | 0.999 | 0.7884 |

Here `emission` is the clear outlier. Its concentration trajectory differs by
9.2%, its full-spectrum absorbance RMSE is approximately 5-12 times larger than
the alternatives, and its forward yield is about 11 percentage points below
`concentrations`. Its pure-reactant initial condition also conflicts with the
0.964 initial reactant fraction recovered by independent decomposition.

The other three methods are much more consistent: their forward yields span
35.1-37.8%, reverse yields span 49.9-50.9%, and concentration trajectories stay
within 2.5% of the established route. Their different fitted initial fractions
(0.964-0.999) remain a visible source of systematic uncertainty.

Files: `455nm_first40__<method>.png` and `.svg`.

## Overall interpretation

1. Similar concentration curves do not guarantee equally well-determined
   quantum yields. This is most visible for the reverse 340 nm yield.
2. The legacy `emission` method can agree well when the LED-band data are
   informative, as in the baselined 340 nm example, but fails the 455 nm
   cross-method comparison.
3. `regularized_concentrations` remains close to the established concentration
   route in all three datasets and avoids independent boundary clipping.
4. `ode_absorbance` gives the smallest full-spectrum residual in all examples,
   but part of that improvement is expected from its fitted baseline terms.
5. Formal errors within one method do not include the full method-selection,
   reference-spectrum, or LED-processing uncertainty. The spread between
   methods is therefore scientifically relevant and should be reported or used
   as a sensitivity check.

A formal goodness-of-fit ranking should wait until absorbance and absorptivity
uncertainties are available. Those uncertainties would allow reduced
chi-squared values and statistically meaningful acceptance intervals rather
than relying on unweighted RMSE.
