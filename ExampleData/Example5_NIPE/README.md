# Example 5 — NIPE with degradation

## What is being measured?

F1 and F2 are two forms of the same photoswitch. In this example, the solution
starts mostly as F2. Irradiation with the 365 nm LED converts F2 into F1, while
the reverse F1-to-F2 photoreaction can also occur:

`F2 ⇌ F1`

In an ideal experiment these would be the only two species, so the total
concentration `[F1] + [F2]` would stay constant. Here it does not: the sample
also changes through degradation or another unmodelled process. This is why a
normal two-species kinetic fit cannot follow the entire experiment well.

## What are the files?

- `M1.dat`, `M2.dat`, and `M3.dat` are three repeat measurements. Each contains
  a sequence of absorption spectra collected during irradiation.
- `F1_MeCN_Epsilon.dat` and `F2_MeCN_Epsilon.dat` are the reference molar
  absorptivity spectra used to estimate how much F1 and F2 is present.
- `LED365nm_baselined.dat` describes the light emitted by the 365 nm LED. It is
  already baseline-corrected, so additional LED baseline correction is off.
- `Power_365nm_25P_1.csv`, `_2.csv`, and `_3.csv` are the three power-monitor
  measurements. `Details.txt` records the power and sample volume used in the
  ready-made configurations.
- `timestamps_30s.csv` assigns one spectrum every 30 seconds, from 0 to 900
  seconds.

## How to run the example

Open `analysis_M1.json`, `analysis_M2.json`, or `analysis_M3.json` in AutoQY
Analysis and select **Run analysis**. The files are already configured for
NIPE, the recorded average power, a 2.07 mL sample, and a 1 cm optical path.

After the run:

1. The main result reports the NIPE fit to the complete 0–900 s trace.
2. **Preprocessing and fit diagnostics** shows a red F2⇌F1 model check. This is
   expected for these data and warns that the complete two-species fit is not a
   faithful description of the sample.
3. Open **NIPE pre-plateau window analysis**. AutoQY compares short early
   sections, before the concentration curve becomes flat, and estimates the
   quantum yield near the start of irradiation where degradation has had less
   time to accumulate.

The early value is still called an **apparent quantum yield**. It is a better
estimate of the F2⇌F1 photochemistry before degradation dominates, but these
data cannot determine a separate quantum yield or mechanism for the degradation
product.

## Why does the M1 fit miss some points?

The kinetic curves assume that one constant pair of quantum yields describes
the whole experiment. M1 progressively departs from that assumption: the
amount represented by the two reference spectra changes and the local NIPE
yield decreases with irradiation time. Consequently, the complete-trace curve
is a compromise between incompatible early and late sections. The red model
check and the early-window result should therefore be reported together rather
than treating the smooth complete-trace curve as a valid degradation model.
