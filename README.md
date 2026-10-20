<p align="center">
  <img src="autoqy_core/assets/autoqy-logo.png" alt="AutoQY" width="760">
</p>

# AutoQY

AutoQY calculates photoisomerization quantum yields and helps prepare optical-power and spectral data.

**Version 2.5.0** — see [what changed since the previous main release](RELEASE-NOTES-2.5.0.md).

If you want to use the program on Windows, follow the quickstart below. 
<details>
<summary><strong>Advanced: manual installation on Windows, Linux, or macOS — click to expand</strong></summary>

Use this method if you prefer to install and start AutoQY from a terminal. You
need Conda and an internet connection. Manual installation does not create the
Windows Desktop shortcuts described above.

### Windows (Anaconda Prompt or Miniconda Prompt)

```powershell
conda create --name autoqy-core python=3.12 pip git
conda activate autoqy-core
git clone --branch main --single-branch https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python -m pip install --upgrade pip
python -m pip install -e ".[gui]"
```

### Linux

```bash
conda create --name autoqy-core python=3.12 pip git
conda activate autoqy-core
git clone --branch main --single-branch https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python -m pip install --upgrade pip
python -m pip install -e '.[gui]'
```

### macOS

```bash
conda create --name autoqy-core python=3.12 pip git
conda activate autoqy-core
git clone --branch main --single-branch https://github.com/CrespiLab/autoQY.git AutoQY-Core
cd AutoQY-Core
python -m pip install --upgrade pip
python -m pip install -e '.[gui]'
```

After installation, start one GUI at a time:

```text
autoqy-analysis-gui
autoqy-smoother-gui
autoqy-power-gui
```

To use AutoQY again later, open a terminal and run `conda activate autoqy-core`
before starting a GUI.

</details>


<details>
<summary><strong>Windows quickstart — click to expand</strong></summary>

<br>

### 1. Install Conda first

You need **[Miniconda](https://www.anaconda.com/download/success)** or **[Anaconda](https://www.anaconda.com/download)**. Install one of them before installing AutoQY.

You also need an internet connection during installation. AutoQY installs its own Python environment automatically.

### 2. Download the installer

**[Download Install-AutoQY.bat](https://github.com/CrespiLab/autoQY/releases/latest/download/Install-AutoQY.bat)**

### 3. Choose where AutoQY will be installed

The installer creates a folder named `AutoQY-Core`.

The easiest method is:

1. Move `Install-AutoQY.bat` into the folder where you want `AutoQY-Core` to be created.
2. Double-click the BAT file.
3. When asked whether to use the current folder, answer **Yes**.

For example, if the BAT file is in `Documents\Software`, AutoQY will be installed in `Documents\Software\AutoQY-Core`.

Alternatively, leave the BAT file in Downloads. When asked whether to use the current folder, answer **No**, then copy and paste the full path of the folder you want.

### 4. Follow the installer

When asked for a Git branch, press Enter to install `main`, or enter another
branch name such as `develop` or `feature/my-change`. Accept the Windows prompt
if one appears, answer the remaining questions, and wait for the installation
to finish. It can take several minutes.

The BAT launcher downloads its PowerShell installer into a temporary folder
and removes that download when it finishes, including after a failed install.
An existing `Install-AutoQY.ps1` beside the BAT file is left untouched.

Advanced users running the PowerShell installer directly can select a branch
without an interactive prompt:

```powershell
.\Install-AutoQY.ps1 -Branch feature/my-change
```

When installation is complete, open the new **AutoQY** folder on your **Desktop**.

The GUI shortcuts run without a terminal window on Windows. After updating
an existing installation, rerun the installer to refresh these shortcuts.
Diagnostic logs are saved in `%LOCALAPPDATA%\AutoQY\logs`.


## Start here after installation

For most work, double-click:

**Desktop → AutoQY → AutoQY Analysis GUI**

The Desktop folder also contains:

- **AutoQY Spectral Treatment** — prepares spectra, calculates molar absorptivity,
  and reconstructs a pure metastable-isomer absorbance spectrum from a stable
  reference and a mixture of known composition, without requiring ε conversion.
  Save a compressed binary `.aqys` session in **Output** and drop it back into
  Spectral Treatment to resume the saved state. Autosave uses the same format.
  CSV content, names, destination, and all save buttons are managed in **Output**.
- **AutoQY Power GUI** — treats Thorlabs optical-power traces.
- **Uninstall AutoQY** — removes AutoQY and can optionally remove its Conda environment.
  If AutoQY is open, save your work and close its windows first. Uninstall waits
  until they are closed; you can type `cancel` without deleting anything.

</details>

## Follow the tutorial

If you want to click buttons right on, start with the **[395 nm step-by-step tutorial](ExampleData/Example-4_395nm-EpsilonError/TUTORIAL.md)**. It shows how to prepare spectra, build an analysis, run it, and inspect the results.

A step-by-step guide of the different GUIs and how to inspect and export the results can be found in the **[Beginners Guide](BEGINNERS-GUIDE.md)**.

More ready-to-run datasets are listed in **[ExampleData](ExampleData/README.md)**.

## Advanced use

Command-line usage, file formats, fitting methods, configuration details, and the Python API are documented in **[autoqy_core/README.md](autoqy_core/README.md)**.

The terminal and direct JSON runner remain available for legacy or automated
workflows, but they are not added to the Desktop. To see the available commands,
open Anaconda Prompt, Miniconda Prompt, PowerShell, or another terminal and run:

```powershell
conda run --name autoqy-core autoqy-core --help
```

## Citation

The scientific method is described in A. Volker, J. D. Steen and S. Crespi, *Beilstein J. Org. Chem.* **2024**, 20, 1684–1692, [doi:10.3762/bjoc.20.150](https://doi.org/10.3762/bjoc.20.150).

Software citation metadata is provided in [CITATION.cff](CITATION.cff).

## Disclaimer

AutoQY is a research tool and may produce incorrect or unsuitable results if the inputs or settings are inappropriate. Users are responsible for checking their data, settings, assumptions, uncertainties, and outputs before reporting or publishing results. The developers accept no responsibility for errors, data loss, or consequences arising from use of the software.


## Analysis composition and fit warnings

In Analysis → Fit → Method-specific controls, enter the **expected metastable
photoproduct at PSS (%)**. For example, enter **77** for 77% photoproduct and
23% remaining stable reactant. Leave it blank if unknown. This is a diagnostic
reference and does not constrain the fitted quantum yields.

Existing schema-v1 JSON files still store `analysis.expected_pss_reactant_percent`.
The GUI converts with `reactant_percent = 100 - product_percent` when saving and
converts back when loading. Do not change existing JSON values to product percent.

Fit checks now explain the problem and what to check directly in the diagnostics
panel. **Review** asks for inspection; **Red flag** means the result may be unreliable.
Jacobian and spectral condition numbers remain under **Technical details**.
The numerical checks and thresholds have not changed.

For startup profiling and validation details, see [the startup notes](docs/startup-2.4.1.md).

## Share a complete experiment as one JSON file

In Analysis → Output, leave **Include all input data in JSON (recommended)**
enabled and click **Save JSON**. The saved file includes the exact measurement
spectra, both absorptivity references (including uncertainty columns), LED data
when needed, timestamps, and analysis settings. You can move or send this one
file without the original input files. Choose **Load existing analysis.json**
on the receiving computer and select a local results directory before running.

With this option enabled, result JSON and configuration snapshots written during
a run also include the complete experiment. A result JSON can be opened with the
same load button or passed to `autoqy-core run`. TXT/CSV/image outputs are still
separate exports; they are not required to rerun the experiment.

Disable the option to save file links for newly selected inputs. Data loaded from
a self-contained file stays embedded unless you replace it using **Browse**.
Existing path-based JSON remains supported. Self-contained files use schema
version 2 and need this AutoQY 2.5.0 build or later; older builds cannot open them.
Embedded files use lossless base64 data and SHA-256 integrity checks. This preserves
binary formats and all numerical precision, at the cost of larger JSON files.
