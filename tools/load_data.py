from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import pandas as pd
import numpy as np

import autoQuant.LoadedData as LoadedData
import tools.load_data as LoadData
from tools.plotting import MplCanvas


def load_powerdata(file_path):
    """Load data from a CSV file."""
    try:
        data = pd.read_csv(file_path, skiprows=14)
        Power = data["Power (W)"] * 1000     # Convert to mW
        x = np.arange(len(Power))            # Index for the x-axis
        return x, Power
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None

def load_and_validate_data(file_name):
    """Load data and handle errors."""
    x, Power = load_powerdata(file_name) # from tools.power_data.py: converts W (raw data) to mW
    if Power is None or x is None:
        print("Failed to load data.")
    return x, Power

def load_power(file_name):
    """Load the data from a file and plot it in a new window."""
    # options = QtWidgets.QFileDialog.Options()
    # file_name, _ = QtWidgets.QFileDialog.getOpenFileName("Load Data File", "",
    #                                                      "Data Files (*.csv);;All Files (*)",
    #                                                      options=options)
    if file_name:
        LoadedData.filename_power = file_name
        LoadedData.x, LoadedData.Power = LoadData.load_and_validate_data(LoadedData.filename_power)
        
        # Save x and Power without overwriting
        LoadedData.loaded_powerdata.append((LoadedData.x, LoadedData.Power))

        

        if LoadedData.Power is not None and LoadedData.x is not None:
            add_new_window_PowerData(plot_power_data, f"Power Data {LoadedData.count+1}", LoadedData.count)
        LoadedData.count += 1
        
def add_new_window_PowerData(self, plot_func, title="New Plot", idx=None, *args):
    """Open a new standalone window with the given plot widget."""
    ##!!! ELABORATE THIS WINDOW WITH BUTTONS: BASELINE CORR; CALCULATE POWER
    try:
        # Create a new window widget
        window = QtWidgets.QWidget()
        window.setWindowTitle(title)

        # Create a layout for the window
        layout = QtWidgets.QVBoxLayout()
        window.setLayout(layout)

        # Create the custom MplCanvas
        canvas = MplCanvas(self, idx)  # Pass idx for initialization

        # Create a navigation toolbar
        toolbar = NavigationToolbar(canvas, self)

        # Add the toolbar and canvas to the layout
        layout.addWidget(toolbar)  # Toolbar at the top
        layout.addWidget(canvas)   # Canvas below the toolbar

        # Call the plotting function to populate the canvas
        if idx is not None:
            print(f'plot nr {idx + 1}')  # Start from 1
            plot_func(canvas, idx, *args)  # Pass idx only if provided
        else:
            print('plot nr None')
            plot_func(canvas, *args)  # No idx passed

        # Show the window
        window.show()

        # Keep a reference to the window to prevent it from being garbage collected
        if not hasattr(self, "_open_windows"):
            self._open_windows = []  # Create a list to store open windows
        self._open_windows.append(window)  # Store the window
    except Exception as e:
        QtWidgets.QMessageBox.critical(self, "Error", f"Failed to create new window: {e}")


def plot_power_data(canvas, idx):
    """Plot the loaded data using the MplCanvas."""
    if LoadedData.x is None or LoadedData.Power is None:
        QtWidgets.QMessageBox.warning("Error", "No data loaded")
        return
    
    LoadedData.line_positions.extend([[0] * 12 for _ in range(idx - len(LoadedData.line_positions) + 1)])
    ### Use the plot_power method from MplCanvas
    canvas.plot_power(LoadedData.x, LoadedData.Power, LoadedData.filename_power)
    LoadedData.line_positions[idx] = get_sorted_line_positions(canvas)

def get_sorted_line_positions(canvas):
    """Get and sort the x-positions of vertical lines."""
    line_positions = [line.get_xdata()[0] for line in canvas.lines]
    line_positions.sort()
    return line_positions

