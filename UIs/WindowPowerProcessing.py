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

import autoQuant.LoadedData as LoadedData
# import tools.load_data as LoadData
from tools.plotting import MplCanvas
import tools.baseline_power as BaselinePower

class WindowPowerProcessing(QtWidgets.QMainWindow):
    """Class for PowerProcessing module."""
    def __init__(self):
        super(WindowPowerProcessing, self).__init__()
        uic.loadUi('UIs/WindowPowerProcessing.ui', self)  # Load the UI file you provided
        
        
        self.x = None
        self.Power = None
        self.loaded_powerdata = [] 
        self.power_baselined_noyes = []
        
        self.line_positions = [] ## list of line positions (this is the parent window of MplCanvas in plotting.py)
        
        ##!!! ELABORATE THIS WINDOW WITH BUTTONS: CALCULATE POWER

        self.WindowPP_baselineCorrectionButton.clicked.connect(self.baseline_correction)
        
        
        options = QtWidgets.QFileDialog.Options()
        self.file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                              "Load Data File", "",
                                                              "Data Files (*.csv);;All Files (*)",
                                                              options=options)

        self.load_power(self.file_name)

    def load_powerdata(self,file_path):
        """Load data from a CSV file."""
        try:
            data = pd.read_csv(file_path, skiprows=14)
            Power = data["Power (W)"] * 1000     # Convert to mW
            x = np.arange(len(Power))            # Index for the x-axis
            return x, Power
        except Exception as e:
            print(f"Error loading data: {e}")
            return None, None
    
    def load_and_validate_data(self,file_name):
        """Load data and handle errors."""
        x, Power = self.load_powerdata(file_name) # from tools.power_data.py: converts W (raw data) to mW
        if Power is None or x is None:
            print("Failed to load data.")
        return x, Power
    
    def load_power(self,file_name):
        """Load the data from a file and plot it in a new window."""
        if file_name:
            # LoadedData.filename_power = file_name
            # LoadedData.x, LoadedData.Power = self.load_and_validate_data(LoadedData.filename_power)
            
            self.filename_power = file_name
            self.x, self.Power = self.load_and_validate_data(self.filename_power)
            
            ## Save x and Power without overwriting
            # LoadedData.loaded_powerdata.append((LoadedData.x, LoadedData.Power))
            self.loaded_powerdata.append((self.x, self.Power))


            if self.Power is not None and self.x is not None:
                self.add_new_tab_PowerData(self.plot_power_data, f"Power Data {LoadedData.count}", LoadedData.count)
            
            print(f"WPP-load_power===LoadedData.count:{LoadedData.count}")
            
            # LoadedData.count += 1
    
            # if LoadedData.Power is not None and LoadedData.x is not None:
            #     self.add_new_tab_PowerData(self.plot_power_data, f"Power Data {LoadedData.count+1}", LoadedData.count)
            # LoadedData.count += 1
    
    def add_new_tab_PowerData(self, plot_func, title, idx=None, *args):
        """Create a new tab with a plot and a navigation toolbar."""
        try:
            tab = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout()  # Use QVBoxLayout to stack toolbar and canvas vertically

            # Create the custom MplCanvas
            canvas = MplCanvas(self, idx)  # Pass idx for initialization

            # Create the navigation toolbar for the canvas
            toolbar = NavigationToolbar(canvas, self)

            # Add the toolbar and canvas to the layout
            layout.addWidget(toolbar)  # Add the toolbar at the top
            layout.addWidget(canvas)   # Add the canvas (plot area) below the toolbar

            tab.setLayout(layout)

            # Call the plotting function to populate the canvas
            if idx is not None:
                print(f'plot nr {idx}')  # Start from 1
                plot_func(canvas, idx, *args)  # Pass idx only if provided
            else:
                print('plot nr None')
                plot_func(canvas, *args)  # No idx passed

            # Add the new tab to the tab widget
            self.tabWidget.addTab(tab, title)

            # Make tabs closable
            self.tabWidget.setTabsClosable(True)

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to add new tab: {e}")
    
    def plot_power_data(self,canvas, idx):
        """Plot the loaded data using the MplCanvas."""
        # if LoadedData.x is None or LoadedData.Power is None:
        #     QtWidgets.QMessageBox.warning("Error", "No data loaded")
        #     return

        if self.x is None or self.Power is None:
            QtWidgets.QMessageBox.warning("Error", "No data loaded")
            return
        
        self.line_positions.extend([[0] * 12 for _ in range(idx - len(self.line_positions) + 1)])
        ### Use the plot_power method from MplCanvas
        
        # canvas.plot_power(LoadedData.x, LoadedData.Power, LoadedData.filename_power)
        canvas.plot_power(self.x, self.Power, self.filename_power)
        self.line_positions[idx] = self.get_sorted_line_positions(canvas)
    
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

        ##!!! HERE: ADD NEW TAB IN POWERDATA-WINDOW

        # Iterate over each loaded dataset using the saved idx (LoadedData.count)
        print(f"WPP-baseline_correction===len(self.loaded_powerdata):{len(self.loaded_powerdata)}")
        # for idx in range(len(self.loaded_powerdata)):
        #     ## Baseline correction: no jacket, no cuvette
        #     self.add_new_tab_PowerData(self.plot_baseline_corrected_no, f"Baseline correction: No Jacket and No Cuvette (Power Data {idx+1})", idx)
            
        #     # Baseline correction: jacket, cuvette with solvent
        #     self.add_new_tab_PowerData(self.plot_baseline_corrected_yes, f"Baseline correction: With Jacket and Cuvette with Solvent (Power Data {idx+1})", idx)

        self.add_new_tab_PowerData(self.plot_baseline_corrected_no, f"Baseline correction {LoadedData.count}: No Jacket and No Cuvette", LoadedData.count)
        print("WPP-baseline_correction===after self.plot_baseline_corrected_no")
        self.add_new_tab_PowerData(self.plot_baseline_corrected_yes, f"Baseline correction {LoadedData.count}: With Jacket and Cuvette with Solvent", LoadedData.count)

    def plot_baseline_corrected_no(self, canvas, idx):
        
        # print(f"WPP-plot_baseline_corrected_no===self.loaded_powerdata:{self.loaded_powerdata}")
        # print(f"WPP-plot_baseline_corrected_no===self.loaded_powerdata[0]:{self.loaded_powerdata[0]}")

        case = 0 # no cuvette
        # x, Power = self.loaded_powerdata[idx]
        # line_positions = self.line_positions[idx]

        # print(f"WPP-plot_baseline_corrected_no===self.line_positions:{self.line_positions}")
        # print(f"WPP-plot_baseline_corrected_no===self.line_positions[1]:{self.line_positions[1]}")
        
        x, Power = self.loaded_powerdata[0]
        # print(f"WPP-plot_baseline_corrected_no===x:{x}")
        # print(f"WPP-plot_baseline_corrected_no===Power:{Power}")

        line_positions = self.line_positions[1]
        
        if len(line_positions) != 12: # Ensure 12 positions are selected
            QtWidgets.QMessageBox.warning(self, "Error", "You need to select exactly 12 line positions.")
            return
        #if idx == 1:
            # Proceed with baseline correction
        sections = BaselinePower.choose_sections(line_positions)
        baselined, baseline, section_on_no = self.calculate_baseline_and_power_no(x, Power, sections)#baseline_correction(Power, x, sections)
        
        self.power_baselined_noyes.append(section_on_no) # power values LED ON, no cuvette, no jacket
        
        """Plot the baseline-corrected data."""
        canvas.plot_baseline_correction(x, Power, baseline, baselined, sections,
                                        case, self.filename_power)


    def plot_baseline_corrected_yes(self, canvas, idx):
    ##!!! MOVE TO WindowPowerProcessing.py
        case = 1 # yes cuvette
        x, Power = self.loaded_powerdata[0]
        line_positions = self.line_positions[1]
        
        if len(line_positions) != 12: # Ensure 12 positions are selected
            QtWidgets.QMessageBox.warning(self, "Error", "You need to select exactly 12 line positions.")
            return

        sections = BaselinePower.choose_sections(line_positions)
        baselined, baseline, section_on_yes = self.calculate_baseline_and_power_yes(x, Power, sections) #baseline_correction(Power, x, sections)

        self.power_baselined_noyes.append(section_on_yes) # power values LED ON, yes cuvette, yes jacket
        
        """Plot the baseline-corrected data."""
        canvas.plot_baseline_correction(x, Power, baseline, baselined, sections,
                                        case, self.filename_power)

    def calculate_baseline_and_power_no(self, x, Power, sections):
        """Calculate baseline and baseline-corrected power."""
        # Baseline correction: no jacket, no cuvette
        x_masked = np.concatenate((x[sections["start_0"]:sections["end_0"]], x[sections["start_0_1"]:sections["end_0_1"]]))
        y_masked = np.concatenate((Power[sections["start_0"]:sections["end_0"]], Power[sections["start_0_1"]:sections["end_0_1"]]))
        x_masked = x_masked[~np.isnan(y_masked)]
        y_masked = y_masked[~np.isnan(y_masked)]

        # Polynomial fit (n = 3)
        n = 3
        p0 = np.full(n, 0.000000001)
        popt, _ = curve_fit(BaselinePower.fit_func, x_masked, y_masked, p0=p0)
        baseline_1 = np.polyval(popt, x)
        baselined_1 = Power - baseline_1

        section_on_no = baselined_1[sections["start_1"]:sections["end_1"]]

        return baselined_1, baseline_1, section_on_no

    def calculate_baseline_and_power_yes(self, x, Power, sections):
        # Baseline correction: jacket, cuvette with solvent
        x_masked = np.concatenate((x[sections["start_2"]:sections["end_2"]], x[sections["start_4"]:sections["end_4"]]))
        y_masked = np.concatenate((Power[sections["start_2"]:sections["end_2"]], Power[sections["start_4"]:sections["end_4"]]))
        x_masked = x_masked[~np.isnan(y_masked)]
        y_masked = y_masked[~np.isnan(y_masked)]
        
        # Polynomial fit (n = 3)
        n = 3  # Degree of polynomial for baseline correction
        p0 = np.full(n + 1, 0.000000001)  # Initial guess of baseline coefficients
        
        popt, _ = curve_fit(BaselinePower.fit_func, x_masked, y_masked, p0=p0)
        baseline_2 = np.polyval(popt, x)
        baselined_2 = Power - baseline_2

        section_on_yes = baselined_2[sections["start_3"]:sections["end_3"]]
        # power_std = np.nanstd(baselined_2[sections["start_3"]:sections["end_3"]])

        return baselined_2, baseline_2, section_on_yes

