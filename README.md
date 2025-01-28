## autoQuant
**January 28th, 2025**

autoQuant is a graphical user interface (GUI) for the calculation of the isomerization quantum yield using data recorded according to the following publication:

A. Volker, J. D. Steen, S. Crespi, A fiber-optic spectroscopic setup for isomerization quantum yield determination, Beilstein J. Org. Chem. 2024, 20, 1684–1692, DOI: 10.3762/bjoc.20.150.

### Installation

#### Create a new Python environment and install pip
Python 3.12 or higher is required
##### For example with Conda:
```bash
conda create -n autoQuant
conda activate autoQuant
conda install pip
```

#### Download the source files
##### - Clone using the URL:
```bash
git clone https://github.com/JornSteen/autoQuant.git
```
A folder called "autoQuant" is downloaded.

##### - Download ZIP folder and unpack at desired location
A folder called "autoQuant-main" is downloaded.

#### Install
##### From cloned folder
```bash
cd autoQuant
pip install .
```
##### From ZIP
```bash
cd autoQuant-main
pip install .
```

### Run
Execute main.py with Python:
#### From cloned folder
```bash
cd autoQuant
python main.py
```
#### From ZIP
```bash
cd autoQuant-main
python main.py
```
The GUI should appear after a short while.
