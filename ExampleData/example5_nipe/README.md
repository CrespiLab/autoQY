# Example 5 — NIPE with degradation

This example contains the original F1/F2 triplicate supplied for the NIPE
implementation. The samples begin predominantly as F2 and form F1 under the
365 nm LED.

Open `analysis_M1.json`, `analysis_M2.json`, or `analysis_M3.json` in AutoQY
Analysis. Each configuration uses the measured power and volume recorded in
`Details.txt`, selects NIPE, and displays the automatic pre-plateau window
analysis after the fit.

The timestamps are an explicit working assumption requested for this example:
one spectrum every 30 seconds, from 0 to 900 seconds. Replace
`timestamps_30s.csv` if recorded acquisition times become available. A 1 cm
path length is also assumed.

The supplied LED spectrum is already baselined, so the configuration leaves
additional LED baseline correction off. The original three power-monitor CSV
files and `Details.txt` are retained unchanged for provenance.

All three measurements are expected to trigger the red A⇌B diagnostic. The
pre-plateau estimate remains an apparent quantum yield because the available
F1 and F2 reference spectra cannot separately identify the degradation-product
yield.
