# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 2024

@author: Jorn Steen

Window for PowerProcessing module

"""

# import sys
from PyQt5 import uic, QtWidgets
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import pandas as pd
import numpy as np

import autoQuant.LoadedData as LoadedData
import tools.load_data as LoadData
from tools.plotting import MplCanvas
import tools.baseline_power as BaselinePower



class WindowPowerProcessing(QtWidgets.QMainWindow):
    """Class for PowerProcessing module."""
    def __init__(self):
        super(WindowPowerProcessing, self).__init__()
        uic.loadUi('UIs/WindowPowerProcessing.ui', self)  # Load the UI file you provided
        
        
        self.line_positions = []
        
        ##!!! ELABORATE THIS WINDOW WITH BUTTONS: BASELINE CORR; CALCULATE POWER

        # self.WindowPP_baselineCorrectionButton.clicked.connect(self.)
        
        
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
        # options = QtWidgets.QFileDialog.Options()
        # file_name, _ = QtWidgets.QFileDialog.getOpenFileName("Load Data File", "",
        #                                                      "Data Files (*.csv);;All Files (*)",
        #                                                      options=options)
        if file_name:
            LoadedData.filename_power = file_name
            LoadedData.x, LoadedData.Power = self.load_and_validate_data(LoadedData.filename_power)
            
            # Save x and Power without overwriting
            LoadedData.loaded_powerdata.append((LoadedData.x, LoadedData.Power))
    
            if LoadedData.Power is not None and LoadedData.x is not None:
                self.add_new_tab_PowerData(self.plot_power_data, f"Power Data {LoadedData.count+1}", LoadedData.count)
            LoadedData.count += 1
    
    def add_new_tab_PowerData(self, plot_func, title, idx=None, *args):
        """Create a new tab with a plot and a navigation toolbar."""
        try:
            tab = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout()  # Use QVBoxLayout to stack toolbar and canvas vertically

            # Create the custom MplCanvas
            canvas = MplCanvas(self, idx)  # Pass idx for initialization

            print(f"WPP-add_new_tab_PowerData === idx:{idx} and type:{type(idx)}")

            # Create the navigation toolbar for the canvas
            toolbar = NavigationToolbar(canvas, self)

            # Add the toolbar and canvas to the layout
            layout.addWidget(toolbar)  # Add the toolbar at the top
            layout.addWidget(canvas)   # Add the canvas (plot area) below the toolbar

            tab.setLayout(layout)

            ############
            # Call the plotting function to populate the canvas
            # idx functionality removed as it's unused in this version
            # plot_func(canvas)
            ############
            
            # Call the plotting function to populate the canvas
            if idx is not None:
                print(f'plot nr {idx + 1}')  # Start from 1
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
        if LoadedData.x is None or LoadedData.Power is None:
            QtWidgets.QMessageBox.warning("Error", "No data loaded")
            return
        
        print(f"WPP-plot_power_data === LoadedData.line_positions: {LoadedData.line_positions} and type {type(LoadedData.line_positions)}")
        
        self.line_positions.extend([[0] * 12 for _ in range(idx - len(self.line_positions) + 1)])
        ### Use the plot_power method from MplCanvas
        canvas.plot_power(LoadedData.x, LoadedData.Power, LoadedData.filename_power)
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
            print("LoadedData.line_positions is empty")

