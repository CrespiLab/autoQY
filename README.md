## autoQuant
**January 28th, 2025**

autoQuant is a graphical user interface (GUI) for the calculation of the isomerization quantum yield using data recorded according to the following publication:

A. Volker, J. D. Steen, S. Crespi, A fiber-optic spectroscopic setup for isomerization quantum yield determination, Beilstein J. Org. Chem. 2024, 20, 1684–1692, DOI: 10.3762/bjoc.20.150.

### Installation
Python 3.12 or higher is required

#### Conda
The Anaconda Powershell Prompt is a good tool.
##### Create a new Python environment and install pip
```bash
(base) conda create -n autoQuant
(base) conda activate autoQuant
(autoQuant) conda install pip
```

##### Download the source files
###### Clone using URL at desired location:
```bash
(autoQuant) cd desired-location
(autoQuant) \desired-location> git clone https://github.com/JornSteen/autoQuant.git
```
A folder called "autoQuant" is downloaded.

###### Download ZIP folder and unpack at desired location
A folder called "autoQuant-main" is unpacked.

##### Install
```bash
(autoQuant) \desired-location> cd autoQuant(-main)
(autoQuant) \desired-location\autoQuant(-main)> pip install .
```

#### Linux
##### Download the source files
###### Clone using URL at desired location:
```bash
~$ cd desired-location
desired-location$ git clone https://github.com/JornSteen/autoQuant.git
```
A folder called "autoQuant" is downloaded.

###### Download ZIP folder and unpack at desired location
A folder called "autoQuant-main" is unpacked.

##### Create a new Python environment and install pip
Create the Python environment in the downloaded autoQuant(-main) folder.
```bash
~$ cd desired-location
desired-location$ sudo apt install python3-venv
desired-location$ python3 -m venv autoQuant(-main)
desired-location$ source autoQuant(-main)/bin/activate
```
If necessary, install pip:
```bash
(autoQuant(-main)) desired-location$ sudo apt install pip
```

##### Install
```bash
(autoQuant(-main)) desired-location$ cd autoQuant(-main)
(autoQuant(-main)) desired-location/autoQuant(-main)$ pip install .
```

### Run
Execute with Python:
```bash
(autoQuant(-main)) cd autoQuant(-main)
(autoQuant(-main)) autoQuant(-main)$ python main.py
```
The GUI should appear after a short while.

N.B. Currently, there is a runtime error upon launching the programme in Linux, which causes re-sizing issues for the window.

### Tested successfully with:
- python 3.12.8

- numpy 2.2.2
- scipy 1.15.1
- lmfit 1.3.2
- pandas 2.2.3
- pyqt 5.15.9 (PyQt5)
- matplotlib 3.10.0
