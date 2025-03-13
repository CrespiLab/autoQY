## autoQY
**February 20th, 2025**

autoQY is a graphical user interface (GUI) for the calculation of the isomerization quantum yield using data recorded according to the following publication:

A. Volker, J. D. Steen, S. Crespi, A fiber-optic spectroscopic setup for isomerization quantum yield determination, Beilstein J. Org. Chem. 2024, 20, 1684–1692, DOI: 10.3762/bjoc.20.150.

### Installation
Python 3.12 or higher is required

#### Conda
The Anaconda Powershell Prompt is a good tool.
##### Create a new Python environment and install pip
```bash
(base) conda create -n autoQY
(base) conda activate autoQY
(autoQY) conda install pip
```

##### Download the source files
###### Clone using URL at desired location:
```bash
(autoQY) conda install git
(autoQY) cd desired-location
(autoQY) \desired-location> git clone https://github.com/CrespiLab/autoQY.git
```
A folder called "autoQY" is downloaded.

###### Download ZIP folder and unpack at desired location
A folder called "autoQY-main" is unpacked.

##### Install
```bash
(autoQY) \desired-location> cd autoQY(-main)
(autoQY) \desired-location\autoQY(-main)> pip install .
```

#### Linux
##### Download the source files
###### Clone using URL at desired location:
```bash
~$ cd desired-location
desired-location$ git clone https://github.com/CrespiLab/autoQY.git
```
A folder called "autoQY" is downloaded.

###### Download ZIP folder and unpack at desired location
A folder called "autoQY-main" is unpacked.

##### Create a new Python environment and install pip
Create the Python environment in the downloaded autoQY(-main) folder.
```bash
~$ cd desired-location
desired-location$ sudo apt install python3-venv
desired-location$ python3 -m venv autoQY(-main)
desired-location$ source autoQY(-main)/bin/activate
```
If necessary, install pip:
```bash
(autoQY(-main)) desired-location$ sudo apt install pip
```

##### Install
```bash
(autoQY(-main)) desired-location$ cd autoQY(-main)
(autoQY(-main)) desired-location/autoQY(-main)$ pip install .
```

### Run
Execute with Python:
```bash
(autoQY(-main)) cd autoQY(-main)
(autoQY(-main)) autoQY(-main)$ python main.py
```
The GUI should appear after a short while.

### Notes
#### Linux
Currently, there is a runtime error upon launching the programme in Linux, which causes re-sizing issues for the window.

#### Windows
In some cases, upon launching autoQuant on a small screen (for example, a laptop), the font sizes are too big, leading to difficulty reading the labels. This is most likely due to a Windows scale setting that increases the size of text.
To change this setting, go to:
- Settings
- Display
- Scale and Layout
- Change the size of text, apps and other items: set to 100% (the default is probably 150%)

Start autoQY again and, hopefully, the GUI looks normal now.

### Tested successfully with:
- python 3.12.8

- numpy 2.2.2
- scipy 1.15.1
- lmfit 1.3.2
- pandas 2.2.3
- pyqt 5.15.9 (PyQt5)
- matplotlib 3.10.0
