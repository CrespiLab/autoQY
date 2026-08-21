# 340 nm baselined-LED example

`generic_inputs/analysis.json` uses only plain, headered CSV files.
`crespi_group_inputs/analysis.json` selects the same direct-absorbance
(`emission`) fit using the original mixed SpectraGryph, headerless CSV, and
Crespi group AHK inputs.

The historical TXT and figure outputs are retained as archival records. 
They are now put in LegacyResults. The current core result is covered by 
the automated emission method. The paired inputs are checked by
`ExampleData/verify_equivalent_inputs.py`.
