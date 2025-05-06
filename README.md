# autoQY
**May 6th, 2025**

autoQY is a graphical user interface (GUI) for the calculation of the isomerization quantum yield using data recorded according to the following publication:

A. Volker, J. D. Steen, S. Crespi, A fiber-optic spectroscopic setup for isomerization quantum yield determination, Beilstein J. Org. Chem. 2024, 20, 1684–1692, DOI: 10.3762/bjoc.20.150.

## Installation
Python 3.12 or higher is required

### Conda
The Anaconda Powershell Prompt is a good tool.
#### Create a new Python environment and install pip
```bash
(base) conda create -n autoQY
(base) conda activate autoQY
(autoQY) conda install pip
```

#### Download the source files
##### Clone using URL at desired location:
```bash
(autoQY) conda install git
(autoQY) cd desired-location
(autoQY) \desired-location> git clone https://github.com/CrespiLab/autoQY.git
```
A folder called "autoQY" is downloaded.

#### Install
```bash
(autoQY) \desired-location> cd autoQY
(autoQY) \desired-location\autoQY> pip install -e .
```

### Linux
#### Download the source files
Clone using URL at desired location:
```bash
~$ cd desired-location
desired-location$ git clone https://github.com/CrespiLab/autoQY.git
```
A folder called "autoQY" is downloaded.

#### Create a new Python environment and install pip
Create the Python environment in the downloaded autoQY(-main) folder.
```bash
~$ cd desired-location
desired-location$ sudo apt install python3-venv
desired-location$ python3 -m venv autoQY
desired-location$ source autoQY/bin/activate
```
If necessary, install pip:
```bash
(autoQY) desired-location$ sudo apt install pip
```

#### Install
```bash
(autoQY) desired-location$ cd autoQY
(autoQY) desired-location/autoQY$ pip install -e .
```

## Run
Make sure to activate the environment and be in the install directory.
```
(base) \autoQY> conda activate autoQY
```
Or
```
/autoQY$ source autoQY/bin/activate
```

Execute with command-line prompt:
```
(autoQY) autoqy
```
The GUI should appear after a short while.

## Notes
### Linux
Currently, there is a runtime error upon launching the programme in Linux, which causes re-sizing issues for the window.

### Windows
In some cases, upon launching autoQuant on a small screen (for example, a laptop), the font sizes are too big, leading to difficulty reading the labels. This is most likely due to a Windows scale setting that increases the size of text.
To change this setting, go to:
- Settings
- Display
- Scale and Layout
- Change the size of text, apps and other items: set to 100% (the default is probably 150%)

Start autoQY again and, hopefully, the GUI looks normal now.
