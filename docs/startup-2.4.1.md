# AutoQY 2.4.1 startup and analysis GUI changes

Base commit: `8810b40` (local branch `2.4.1`).

## Findings and fixes

The package initializer eagerly imported the pipeline, solvers, output writers,
Matplotlib, and runner for every entry point. Consequently, even opening a GUI
loaded unrelated scientific modules. The analysis GUI also imported the runner
and spectral solvers before they were needed. The runner loaded Matplotlib even
when static figure export was disabled. First-time Matplotlib font discovery and
the large import tree added work before a window could appear.

The package now resolves its existing public exports on first access. Analysis
loads the runner and spectral solvers when checking or processing data; the
runner loads static plotting only when exporting figures. Other GUIs also benefit
from the lighter package initializer, although their startup times were not benchmarked.

Five blank analysis plots previously built and validated Plotly figure/template
objects at startup. They now use equivalent figure dictionaries; actual result
plots retain the existing Plotly implementation. The browser launcher also no
longer waits an extra fixed second after the HTTP server has bound its socket.

The profile showed substantial import-related filesystem work. Loading fewer
modules should particularly help machines with slow storage or CPUs, but no older
PC or antivirus-specific measurements were available. The launcher creates a new
isolated browser profile each time; this remains a possible separate contributor
to window-opening time and was left unchanged to preserve window/process cleanup.

## Measured startup

Windows, Python 3.12.10, existing local dependencies, same source folder:

| Measurement | Baseline first run | Baseline repeat | Changed first run | Changed repeat |
| --- | ---: | ---: | ---: | ---: |
| Import analysis GUI | 11.052 s | 7.102 s | 2.066 s | 1.615 s |
| Construct application | 5.222 s | 2.258 s | 1.761 s | 1.242 s |
| Initial HTTP requests | 0.296 s | 0.091 s | 0.125 s | 0.073 s |
| Total | 16.571 s | 9.451 s | 3.952 s | 2.930 s |
| Loaded modules | 2062 | 2061 | 1108 | 1108 |
| Matplotlib pyplot loaded | Yes | Yes | No | No |

These are cProfile-instrumented observations, not controlled cold-disk benchmarks.
Cache warming, OneDrive, machine load, and profiling overhead affect timings. One
changed measurement overlapped a test run. The three initial HTTP endpoints all
returned 200. These timings exclude interpreter startup, Conda activation, browser
launch/rendering, and browser-triggered callbacks. Numerical imports are deferred,
so their remaining cost is paid on the first operation that needs them.

To repeat on an affected PC, activate the AutoQY environment and run from this
repository in separate fresh processes:

```powershell
python tools/profile_gui_startup.py --gui analysis --output startup-1.json
python tools/profile_gui_startup.py --gui analysis --output startup-2.json
```

The tool writes timing/environment JSON and a cumulative profile text file. It
does not open a browser or start a persistent server. `--gui power` and `--gui smoother`
are also supported. Use the same dependencies, filesystem location, and cache
conditions when comparing versions. Browser launch should be timed separately
using the normal desktop shortcut.

## Composition and warnings

The analysis form accepts expected metastable photoproduct percent, from 0 to 100.
Schema-v1 JSON retains the reactant-percent field: entering 77 saves 23; loading
an existing value of 23 displays 77. Blank values remain blank. Invalid percentages
are rejected server-side. Expected PSS is a diagnostic reference, not a fit input
or concentration constraint. Spectral Treatment already used metastable percent,
including the complement required by its older reconstruction calculation.

PSS diagnostics now compare product percentages. Existing solver equations and
warning thresholds are unchanged. Fit and spectral warnings explain their meaning
and next steps in visible text. Condition numbers are available under collapsed
Technical details. Method-comparison warnings use the same plain terminology.

## Validation

Full suite: **171 tests passed, 35 subtests passed**. Reported warnings were existing
Dash DataTable deprecation notices and reference-data negative-absorptivity warnings.
Tests cover fresh-process lazy imports, public exports, initial HTTP responses,
browser scheduling, percentage endpoints/blank/invalid values, old JSON loading,
save/load round trips, visible labels, and fit-warning boundaries/failure cases.

Example 4 was run through both the old configuration and the new GUI conversion;
quantum yields and extrapolated PSS agreed. A separate static-export smoke test
successfully generated nonempty PNG and SVG files after deferring Matplotlib.
Its fitted yields were 0.2501211870 and 0.2112826573 (fractions).

The tests used the local Python 3.12.10 runtime with NumPy 2.5.3, SciPy 1.18.1,
pandas 3.0.5, Matplotlib 3.11.2, Dash 3.4.0, and Plotly 6.9.0. No dependency
versions were changed. The package remains version 2.4.1; no installer was rebuilt.
