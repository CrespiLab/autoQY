# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 2024

@author: Jorn Steen

Window for PowerProcessing module

"""
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit #change in baseline correction

from PyQt5 import uic, QtWidgets
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

import QY.LoadedData as LoadedData
from tools.plotting import MplCanvas
import PowerProcessing.baseline_power as BaselinePower
import QY.Results as Results

class WindowPowerProcessing(QtWidgets.QMainWindow):
    """Class for PowerProcessing module."""
    def __init__(self, parent):
        self.parent = parent
        super(WindowPowerProcessing, self).__init__()
        uic.loadUi('UIs/WindowPowerProcessing.ui', self)  # Load the UI file you provided
        
        ##!!! RE-NAME WINDOW TITLE
        
        self.x = None
        self.Power = None
        self.loaded_powerdata = [] 
        self.section_LEDon_no = None
        self.section_LEDon_yes = None
        self.line_positions = [] ## list of line positions (this is the parent window of MplCanvas in plotting.py)
        
        self.WindowPP_baselineCorrectionButton.clicked.connect(self.baseline_correction)
        self.WindowPP_calculatePowerButton.clicked.connect(self.calculate_power_at_cuvette)
        
        options = QtWidgets.QFileDialog.Options()
        self.filename_power, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                              "Load Data File", "",
                                                              "Data Files (*.csv);;All Files (*)",
                                                              options=options)

        self.load_power(self.filename_power)

    def load_power(self,file_name):
        """Load the data from a file and plot it in a new window."""
        if file_name:
            self.x, self.Power = self.load_and_validate_data(file_name)
            
            ## Save x and Power without overwriting
            self.loaded_powerdata.append((self.x, self.Power))

            ## Define filename for Results_PowerProcessing file
            path_split=file_name.split('/')[0:-1] # leave only filepath (remove name)
            path='\\'.join(path_split) # re-join into string
            end="Results_PowerProcessing" # simple name with added info
            Results.savefilename_power = f"{path}\\{end}"

            if self.Power is not None and self.x is not None:
                self.add_new_tab_PowerData(self.plot_power_data, f"Power Data {LoadedData.count}")
    
    def load_and_validate_data(self,file_name):
        """Load data and handle errors."""
        x, data_power = self.load_powerdata(file_name) # from tools.power_data.py: converts W (raw data) to mW
        if data_power is None or x is None:
            print("Failed to load data.")
        return x, data_power

    def load_powerdata(self,file_path):
        """Load data from the CSV file generated by the Thorlabs software Optical Power Monitor (OPM)."""
        try:
            data = pd.read_csv(file_path, skiprows=14)
            data_power = data["Power (W)"] * 1000     # Convert to mW
            x = np.arange(len(data_power))            # Index for the x-axis
            return x, data_power
        except Exception as e:
            print(f"Error loading data: {e}")
            return None, None
            
    def add_new_tab_PowerData(self, plot_func, title, *args):
        """Create a new tab with a plot and a navigation toolbar."""
        try:
            tab = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout()  # Use QVBoxLayout to stack toolbar and canvas vertically

            ## Create the custom MplCanvas
            canvas = MplCanvas(self, LoadedData.count)  # Pass dataset # for initialization

            ## Create the navigation toolbar for the canvas
            toolbar = NavigationToolbar(canvas, self)

            ## Add the toolbar and canvas to the layout
            layout.addWidget(toolbar)  # Add the toolbar at the top
            layout.addWidget(canvas)   # Add the canvas (plot area) below the toolbar

            tab.setLayout(layout)

            ## Call the plotting function to populate the canvas
            if LoadedData.count is not None:
                print(f'plot nr {LoadedData.count}')  # print dataset #
                plot_func(canvas, *args) # call plot_power_data func
            else:
                print('plot nr None')
                plot_func(canvas, *args) ##!!! still needed?

            ## Add the new tab to the tab widget
            self.tabWidget.addTab(tab, title)

            ## Make tabs closable
            self.tabWidget.setTabsClosable(True)

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to add new tab: {e}")
    
    def plot_power_data(self, canvas):
        """Plot the loaded data using the MplCanvas."""
        if self.x is None or self.Power is None:
            QtWidgets.QMessageBox.warning("Error", "No data loaded")
            return
        
        ##!!! remove the use of idx (don't need it anymore, I think)
            ### replaced by LoadedData.count: but this seems overly complicated
        
        # self.line_positions.extend([[0] * 12 for _ in range(idx - len(self.line_positions) + 1)])
        self.line_positions.extend([[0] * 12 for _ in range(LoadedData.count - len(self.line_positions) + 1)])
        
        ### Use the plot_power method from MplCanvas
        canvas.plot_power(self.x, self.Power, self.filename_power)
        self.line_positions[LoadedData.count] = self.get_sorted_line_positions(canvas)
    
    def get_sorted_line_positions(self,canvas):
        """Get and sort the x-positions of vertical lines."""
        line_positions = [line.get_xdata()[0] for line in canvas.lines]
        line_positions.sort()
        return line_positions
    
    def update_display(self, idx, line_index, new_x):
        if self.line_positions is not None:
            self.line_positions[idx][line_index] = new_x 
        else:
            print("self.line_positions is empty")

    ############################################################################################################
    ############################ BASELINE CORRECTION ############################
    ############################################################################################################

    def baseline_correction(self):
        """Apply baseline correction to each loaded dataset and create a new tab for each."""
        if not self.loaded_powerdata:
            QtWidgets.QMessageBox.warning(self, "Error", "No data loaded")
            return

        self.add_new_tab_PowerData(self.plot_baseline_corrected_no, f"Baseline correction {LoadedData.count}: No Jacket and No Cuvette")
        self.add_new_tab_PowerData(self.plot_baseline_corrected_yes, f"Baseline correction {LoadedData.count}: With Jacket and Cuvette with Solvent")

    def plot_baseline_corrected_no(self, canvas):
        case = 0 # no cuvette
        x, data_power = self.loaded_powerdata[0]
        line_positions = self.line_positions[LoadedData.count]
        
        if len(line_positions) != 12: # Ensure 12 positions are selected
            QtWidgets.QMessageBox.warning(self, "Error", "You need to select exactly 12 line positions.")
            return
        sections = BaselinePower.choose_sections(line_positions)
        baselined, baseline, self.section_LEDon_no = self.calculate_baseline_and_power_no(x, data_power, sections)
        
        """Plot the baseline-corrected data."""
        canvas.plot_baseline_correction(x, data_power, baseline, baselined, sections,
                                        case, self.filename_power)

    def plot_baseline_corrected_yes(self, canvas):
        case = 1 # yes cuvette
        x, data_power = self.loaded_powerdata[0]
        line_positions = self.line_positions[LoadedData.count]
        
        if len(line_positions) != 12: # Ensure 12 positions are selected
            QtWidgets.QMessageBox.warning(self, "Error", "You need to select exactly 12 line positions.")
            return

        sections = BaselinePower.choose_sections(line_positions)
        baselined, baseline, self.section_LEDon_yes = self.calculate_baseline_and_power_yes(x, data_power, sections)

        """Plot the baseline-corrected data."""
        canvas.plot_baseline_correction(x, data_power, baseline, baselined, sections,
                                        case, self.filename_power)

    def calculate_baseline_and_power_no(self, x, data_power, sections):
        """Calculate baseline and baseline-corrected power."""
        ## Baseline correction: no jacket, no cuvette
        x_masked = np.concatenate((x[sections["start_0"]:sections["end_0"]], x[sections["start_0_1"]:sections["end_0_1"]]))
        y_masked = np.concatenate((data_power[sections["start_0"]:sections["end_0"]], data_power[sections["start_0_1"]:sections["end_0_1"]]))
        x_masked = x_masked[~np.isnan(y_masked)]
        y_masked = y_masked[~np.isnan(y_masked)]

        ## Polynomial fit (n = 3)
        n = 3
        p0 = np.full(n, 0.000000001)
        popt, _ = curve_fit(BaselinePower.fit_func, x_masked, y_masked, p0=p0)
        baseline = np.polyval(popt, x)
        baselined = data_power - baseline

        section_on_no = baselined[sections["start_1"]:sections["end_1"]]

        return baselined, baseline, section_on_no

    def calculate_baseline_and_power_yes(self, x, data_power, sections):
        ## Baseline correction: jacket, cuvette with solvent
        x_masked = np.concatenate((x[sections["start_2"]:sections["end_2"]], x[sections["start_4"]:sections["end_4"]]))
        y_masked = np.concatenate((data_power[sections["start_2"]:sections["end_2"]], data_power[sections["start_4"]:sections["end_4"]]))
        x_masked = x_masked[~np.isnan(y_masked)]
        y_masked = y_masked[~np.isnan(y_masked)]
        
        ## Polynomial fit (n = 3)
        n = 3  # Degree of polynomial for baseline correction
        p0 = np.full(n + 1, 0.000000001)  # Initial guess of baseline coefficients
        
        popt, _ = curve_fit(BaselinePower.fit_func, x_masked, y_masked, p0=p0)
        baseline = np.polyval(popt, x)
        baselined = data_power - baseline

        section_on_yes = baselined[sections["start_3"]:sections["end_3"]]

        return baselined, baseline, section_on_yes


    def calculate_power_at_cuvette(self):
        """Calculate the power at the cuvette and the standard deviation."""
        if self.section_LEDon_no is None or self.section_LEDon_yes is None:
            QtWidgets.QMessageBox.warning(self, "Error", "Please perform baseline correction first.")
            return
        
        sections_noyes = [self.section_LEDon_no, self.section_LEDon_yes]
        
        if not sections_noyes:
            QtWidgets.QMessageBox.warning(self, "Error", "No baseline correction has been applied yet.")
            return

        means = []
        stds = []

        for baselined_data in sections_noyes:
            baselined_data = np.asarray(baselined_data).ravel()
            if len(baselined_data) == 0:
                continue
            mean_power = np.nanmean(baselined_data)
            std_power = np.nanstd(baselined_data)
            means.append(mean_power)
            stds.append(std_power)

        power_at_cuv = np.mean(means) # average over powers without and with cuvetteholder+cuv to get power at cuvette
        variance_at_cuv = np.sum(np.square(stds)) / len(means)
        std_at_cuv = np.sqrt(variance_at_cuv)

        LoadedData.PowersAtCuvette[LoadedData.count] = power_at_cuv # calculated power
        LoadedData.ErrorsAtCuvette[LoadedData.count] = std_at_cuv # calculated error

        self.parent.labels_power[LoadedData.count].setPlainText(f"{LoadedData.PowersAtCuvette[LoadedData.count]:.2f}") # calculated power
        self.parent.labels_error[LoadedData.count].setPlainText(f"{LoadedData.ErrorsAtCuvette[LoadedData.count]:.2f}") # calculated power

        


        
