# AutoQY 2.5.0

Release prepared: 20 October 2026.

AutoQY 2.5.0 lets you save and share a complete analysis in one JSON file, with
no need to send the original input files separately. It also brings faster analysis
startup, a more complete Spectral Treatment workflow, recoverable sessions,
clearer fit warnings, and more robust Windows installation and power-data handling.

## What changed since the previous main release

This comparison uses the locally available `main` commit `0a99217` (package version
2.4.0) and the `2.4.1` branch commit `9b58106`. No remote fetch was performed.
The development branch changes 39 files, with 4,615 added and 2,157 removed lines,
before release preparation. The self-contained JSON feature was added during
release preparation and is included in the 2.5.0 source and wheel.

### New: complete experiments in one JSON file

- **Include all input data in JSON (recommended)** is enabled by default in
  **Analysis → Output**. Click **Save JSON** to create a self-contained experiment.
- The file includes measurement spectra, reactant and product absorptivity
  references, uncertainty columns, LED emission when required, timestamps, and
  analysis settings. Exact file contents are preserved without rounding or losing
  binary-format information.
- Move or send this single JSON file, then use **Load existing analysis.json** on
  the receiving computer. Choose a local results directory and run the analysis;
  the original input files and their folder locations are no longer required.
- With the option enabled, **result JSON and configuration snapshots** also carry
  the complete experiment. Result JSON can be reopened in the Analysis GUI or
  passed to `autoqy-core run`.
- Existing path-based JSON remains supported. Disable the option to save file
  links for newly selected inputs. Already embedded data stays included unless
  replaced using **Browse**.
- Embedded inputs have integrity checks to detect damaged data. Including the
  data makes JSON files larger. Recipients need this updated 2.5.0 build or later;
  older builds cannot read the self-contained format.
- Separate TXT, CSV, and image exports are still available. They do not need to
  accompany the JSON to reproduce the analysis.

### Faster startup and clearer analysis

- Scientific solvers and static plotting are loaded when needed, instead of during
  every GUI launch. Empty analysis plots also avoid expensive startup validation,
  and the browser launcher no longer adds a fixed one-second delay.
- Local instrumented analysis-startup measurements improved from 9.45 seconds on
  a repeated baseline run to 2.93 seconds on a repeated changed run. This measures
  imports, app construction, and initial HTTP responses, excluding browser launch
  and rendering. It is not a guarantee for other PCs. A reusable startup profiler
  is included in `tools/profile_gui_startup.py`.
- The optional expected PSS input now asks for **metastable photoproduct (%)**.
  For 77% product and 23% remaining stable reactant, enter **77**. PSS diagnostics
  also display product percentages.
- Fit warnings explain the practical problem and what to check. Jacobian and
  spectral condition numbers remain available under **Technical details** rather
  than leading the warning. The numerical warning thresholds are unchanged.

### Spectral Treatment

- A unified reconstruction panel recovers the pure metastable absorbance spectrum
  from a stable reference, a mixture/PSS spectrum, and the metastable percentage.
- Calibrated reconstruction uses stable-isomer molar absorptivity references.
  Independent spectral and composition uncertainties can be propagated when the
  required errors are available; unknown uncertainties are not silently set to zero.
- Wavelength peak picking adds labelled maxima, with prominence and minimum
  separation controls, without modifying the exported numerical spectrum.
- Improved spectrum selection, ordering, removal, editable names, and consistent
  processing when additional files are loaded. Removing spectra preserves the
  remaining acquisition times.
- A full data sheet follows plotted spectra and supports cell, row, column, and
  whole-sheet copying. CSV exports preserve user-facing labels.
- Plot and export controls offer color-blind-friendly and minimal palettes, legend
  selection, and control over titles and axis labels in PNG/SVG exports.
- Large-dataset controls avoid building unnecessary per-spectrum widgets; previews
  limit the displayed spectra while keeping full data available for processing.

### Session recovery and Windows reliability

- Spectral sessions autosave locally and can be saved and reopened as `.aqys`
  files, including data and processing settings. Atomic writes protect the previous
  recovery file if saving fails; older JSON autosaves remain readable.
- Browser/process tracking and longer close/reload grace periods prevent ordinary
  page reloads or delayed heartbeats from prematurely stopping a GUI. Windowless
  Windows launchers retain diagnostics in local log files.
- GUI shortcuts use dedicated Windows GUI entry points. JavaScript assets needed
  for clipboard and recovery behavior are included in package data.
- The installer downloads its PowerShell bootstrap into a temporary folder and
  cleans up those temporary files without replacing a repository-local installer.
- Uninstallation checks for programs using AutoQY files, lets users close them or
  cancel, and preserves the installation and shortcuts if environment removal fails.

### Power processing and documentation

- Isolated non-finite power readings are excluded with a warning while preserving
  sample positions and region boundaries. Invalid baseline/measurement regions
  produce explicit errors instead of unreliable statistics.
- The beginner guide, package reference, and 395 nm example tutorial have been
  rewritten around the current workflow.
- Regression coverage now includes reconstruction, recovery, peak picking,
  invalid power samples, startup, clearer warnings, composition conversion,
  browser lifecycle, installation, and uninstallation.

## Compatibility and upgrade notes

- Python 3.12 or later is still required; declared dependency requirements are
  unchanged from the previous main branch.
- Existing schema-v1 analysis JSON remains compatible. The stored field is still
  `analysis.expected_pss_reactant_percent`; the GUI converts it using
  `reactant_percent = 100 - product_percent`. Do not replace old JSON values with
  product percentages manually.
- Self-contained JSON uses schema version 2 with lossless embedded file contents
  and integrity checks, preserving uncertainty columns and binary formats. It
  requires this 2.5.0 build or later; prior builds do not support embedded inputs.
- Expected PSS remains a diagnostic reference; it does not constrain the fit.
  The existing kinetic solvers and fit equations are unchanged. New reconstruction
  capabilities and handling of previously rejected invalid power readings are
  separate processing changes.
- Reinstall/update the package to regenerate Windows GUI shortcuts and installed
  package metadata. Merely copying a changed version file does not update an
  existing installed distribution.
- The package version, citation version, and GUI-displayed version are now 2.5.0.
  The citation date is 20 October 2026; historical 2.4.1 profiling references remain
  labelled with the version on which those measurements were made.

## Validation

The prepared source passed **180 Python tests and 39 subtests**, plus the
JavaScript clipboard checks. Nine new portability tests cover GUI save/load,
reopening result JSON, preserving binary contents and uncertainty columns,
detecting damaged embedded data, and compatibility with older file-link inputs.
The sharing test deletes the original input files, reopens the JSON from another
folder, and reproduces the same input arrays and fitted quantum yields.

A separate numerical example confirms that 77% product maps to the former
23% reactant input without changing fitted results. The rebuilt wheel was checked
against the current implementation. See `RELEASE-VALIDATION-2.5.0.md` for details.

No release has been published, no remote branch has been changed, and no new
version-specific DOI is asserted by this preparation.
