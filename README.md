<p align="center">
  <img src="autoqy_core/assets/autoqy-logo.png" alt="AutoQY" width="760">
</p>

# AutoQY

AutoQY calculates photoisomerization quantum yields and helps prepare optical-power and spectral data.

If you just want to use the program on Windows, follow the quickstart below. You can ignore the rest for now.

## Windows quickstart

### 1. Install Conda first

You need **Miniconda or Anaconda**. Install one of them before installing AutoQY.

You also need an internet connection during installation. AutoQY installs its own Python environment, so you do not need to install Python separately.

### 2. Download the installer

**[Download Install-AutoQY.bat](https://github.com/CrespiLab/autoQY/releases/download/v2.3.0/Install-AutoQY.bat)**

### 3. Choose where AutoQY will be installed

The installer creates a folder named `AutoQY-Core`.

The easiest method is:

1. Move `Install-AutoQY.bat` into the folder where you want `AutoQY-Core` to be created.
2. Double-click the BAT file.
3. When asked whether to use the current folder, answer **Yes**.

For example, if the BAT file is in `Documents\Software`, AutoQY will be installed in `Documents\Software\AutoQY-Core`.

Alternatively, leave the BAT file in Downloads. When asked whether to use the current folder, answer **No**, then copy and paste the full path of the folder you want.

### 4. Follow the installer

Accept the Windows prompt if one appears, answer the questions, and wait for the installation to finish. It can take several minutes.

When installation is complete, open the new **AutoQY** folder on your **Desktop**.

## Start here after installation

For most work, double-click:

**Desktop → AutoQY → AutoQY Analysis GUI**

The Desktop folder also contains:

- **AutoQY Spectral Treatment** — prepares spectra and calculates molar absorptivity.
- **AutoQY Power GUI** — treats Thorlabs optical-power traces.
- **AutoQY Analyze JSON** — runs an existing analysis file.
- **AutoQY Terminal** — opens a command line with AutoQY ready.
- **Uninstall AutoQY** — removes AutoQY and can optionally remove its Conda environment.

The GUI window and its terminal belong together. Closing either one closes the other.

## Follow the tutorial

Start with the **[395 nm step-by-step tutorial](ExampleData/Example-4_395nm-EpsilonError/TUTORIAL.md)**. It shows how to prepare spectra, build an analysis, run it, and inspect the results.

More ready-to-run datasets are listed in **[ExampleData](ExampleData/README.md)**.

## The three GUIs

- **Analysis GUI:** creates or loads an analysis, runs the quantum-yield calculation, and shows the results.
- **Spectral Treatment:** corrects, smooths, and exports spectra. It can also calculate reactant and product molar absorptivity.
- **Power GUI:** processes optical-power measurements and reports the power and uncertainty to enter in an analysis.

Important: **Run analysis** does not save the editable `analysis.json`. Use **Save JSON** when you want to keep or reuse the setup.

All GUIs run locally on your computer. Your uploaded data is not sent to an AutoQY server.

## Advanced use

Command-line usage, file formats, fitting methods, configuration details, and the Python API are documented in **[autoqy_core/README.md](autoqy_core/README.md)**.

To see the available commands, open **AutoQY Terminal** and run:

```powershell
autoqy-core --help
```

## Citation

The scientific method is described in A. Volker, J. D. Steen and S. Crespi, *Beilstein J. Org. Chem.* **2024**, 20, 1684–1692, [doi:10.3762/bjoc.20.150](https://doi.org/10.3762/bjoc.20.150).

Software citation metadata is provided in [CITATION.cff](CITATION.cff).

## Disclaimer

AutoQY is a research tool and may produce incorrect or unsuitable results if the inputs or settings are inappropriate. Users are responsible for checking their data, settings, assumptions, uncertainties, and outputs before reporting or publishing results. The developers accept no responsibility for errors, data loss, or consequences arising from use of the software.
