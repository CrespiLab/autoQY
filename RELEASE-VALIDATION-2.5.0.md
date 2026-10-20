# AutoQY 2.5.0 release validation

Validated on 20 October 2026 using the prepared source exported from development
commit `9b58106`, with version and documentation updates for release 2.5.0.
The reference main branch was local commit `0a99217`, package version 2.4.0.
No remote fetch or publication was performed.

- Full Python suite: **180 passed**, with **39 subtests passed**, in 85.25 seconds.
- The suite emitted 76 warnings: Dash DataTable deprecation notices and known
  negative-absorptivity warnings from example reference data. No test failed.
- JavaScript clipboard tests passed for full sheets, row/column intersections,
  ranges, quoting, and numerical precision.
- Package metadata, citation version, and GUI version helper all report **2.5.0**.
  Existing GUI-layout tests also verified version displays in all three GUIs.
- Citation release date is **2026-10-20**. The existing DOI was preserved; no new
  version-specific DOI was invented.
- A wheel was built successfully. Its metadata version, GUI entry points, and
  clipboard/recovery JavaScript assets were inspected and verified.
- Active version references were scanned. Remaining 2.4.1 mentions describe the
  development branch and historical startup measurements.
- The beginner guide and 395 nm tutorial now agree with the GUI's metastable
  photoproduct percentage input; example JSON retains its compatible reactant field.

Python 3.12.10; NumPy 2.5.3; SciPy 1.18.1; pandas 3.0.5; Matplotlib 3.11.2;
Dash 3.4.0; Plotly 6.9.0. No dependency requirements were changed.

## Self-contained JSON feature

Nine new portability tests and four subtests passed. Coverage includes GUI save
and load callbacks, reopening from a different folder after deleting every original
input file, identical input arrays and fitted yields, preserved uncertainty
columns, self-contained result/configuration outputs, lossless binary data,
checksum rejection, chemical-actinometer mode without an LED file, and legacy
path-based configurations. A Windows read-only-directory attribute issue in the
test fixture was corrected by copying file contents without source-directory
attributes; the deleted-originals test then passed.

The wheel was rebuilt and its three changed runtime modules were compared
byte-for-byte with the prepared source. The shared Example 4 JSON retains the
original uncertainty columns and deterministic-extremes configuration. It validates
as schema version 2 without external input paths.

The Windows installer was not executed against a live installation, and the wheel
was not installed into the user's existing environment. Startup measurements in
the release notes are the earlier instrumented observations, not a new older-PC test.

## Prepared files and scope

The prepared source copy is `AutoQY-Core-2.5.0` in the active workspace. The changes
have also been applied directly to
`C:\Users\stecr\OneDrive\Documenti\GitHub\AutoQY-Core` on 20 October 2026,
after checking that the original checkout was clean at commit `9b58106`.
The checkout remains on branch `2.4.1`, with package version 2.5.0 and uncommitted
release changes. No commit or push was performed.

`2.4.1-to-2.5.0.patch` contains the incremental release/version/documentation edits
against commit `9b58106`; apply it to a matching clean 2.4.1 checkout. The source
ZIP and wheel already contain the complete 2.5.0 code. `main-to-2.4.1.diff` records
the accumulated development differences from the locally available old main.

The citation date assumes the proposed 20 October release date. Adjust it if
publication occurs on a later date. No remote changes, tags, or GitHub releases
were created.
