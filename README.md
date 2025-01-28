## autoQuant
**January 28th, 2025**

autoQuant is a graphical user interface (GUI) for the calculation of the isomerization quantum yield using data recorded according to the following publication:

A. Volker, J. D. Steen, S. Crespi, A fiber-optic spectroscopic setup for isomerization quantum yield determination, Beilstein J. Org. Chem. 2024, 20, 1684–1692, DOI: 10.3762/bjoc.20.150.

### Installation

#### Create a new Python environment and install pip
Python 3.12 or higher is required
##### Conda
```bash
conda create -n autoQuant
conda activate autoQuant
conda install pip
```
##### Linux
```bash
sudo apt install python3-venv
python3 -m venv autoQuant
source autoQuant/bin/activate
```

#### Download the source files
##### Clone using URL at desired location:
```bash
git clone https://github.com/JornSteen/autoQuant.git
```
A folder called "autoQuant" is downloaded.

##### Download ZIP folder and unpack at desired location
A folder called "autoQuant-main" is unpacked.

#### Install
```bash
cd autoQuant(-main)
pip install .
```

### Run
Execute with Python:
```bash
cd autoQuant(-main)
python main.py
```
The GUI should appear after a short while.

### Tested successfully with:
- python 3.12.8

- numpy 2.2.2
- scipy 1.15.1
- lmfit 1.3.2
- pandas 2.2.3
- pyqt 5.15.9 (PyQt5)
- matplotlib 3.10.0
