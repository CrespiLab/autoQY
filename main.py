"""
Created on Mon Feb 12 11:03:16 2024
@author: Anouk
Modified by: Jorn, Aug 23rd 2024

This script is used for processing of the data measured with the 
    Optical Power Monitor (OPM) software.
- Baseline correction
- Fit to get value of power with standard deviation

The data is recorded in the following sequence:
- no Luma40 jacket, no cuvette, LED off
- no Luma40 jacket, no cuvette, LED on
- no Luma40 jacket, no cuvette, LED off
- Luma40 jacket, cuvette with solvent, LED off
- Luma40 jacket, cuvette with solvent, LED on
- Luma40 jacket, cuvette with solvent, LED off

Interactive Feature:
- Pick the x-coordinates of the areas for baseline correction interactively

TO DO:
[] add feature: pick two or more input files; 
   then perform averaging and error calculation to yield final result


GUI features to add:
[] Pop-up window to select file

"""
import sys, os
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from PyQt5 import QtWidgets, uic
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import autoQuant.Integration as Integration
import autoQuant.ExpParam as ExpParam
from tools.power_data import load_data
from tools.plotting import MplCanvas
from tools.processing import choose_sections
#from tools.style import apply_dark_theme
from scipy.optimize import curve_fit #change in baseline correction

def fit_func(x, *coeffs):  #not used
    return np.polyval(coeffs, x)

class PowerProcessingApp(QtWidgets.QMainWindow):
        ############################
        #Error Handling  ***********
        ############################
    def __init__(self):
        super(PowerProcessingApp, self).__init__()
        uic.loadUi('qt.ui', self)  # Load the UI file you provided


        ############################
        #Change style  *************
        ############################
        #self.setStyleSheet(get_stylesheet())

        self.filename = None
        self.RefPower = None
        self.x = None

        self.count = 0

        self.led_file = None
        self.spectral_data_path = './'
        self.spectral_data_file = 'output_AutoQuant'
        self.eps_file_a = None
        self.eps_file_b = None

        self.loaded_data = [] 
        self.line_positions = [] 
        self.all_corrected_power = []

        self.LEDindex_first, self.LEDindex_last = None, None
        self.emission_wavelengths = None
        self.emission_Intensity = None
        self.epsilon_A_wavelengths = None 
        self.epsilon_A_values = None
        self.epsilon_B_wavelengths = None
        self.epsilon_B_values = None
        self.SpectralData_Wavelengths = None
        self.SpectralData_Abs = None
        self.SpectralData_Index = None
        self.emission_interp = None
        self.emission_Intensity_proc = None

        self.wavelength_low = None
        self.wavelength_high = None
        self.epsilon_A_interp = None
        self.epsilon_B_interp = None

        # Button connections
        self.loadDataButton.clicked.connect(self.load_power)
        self.baselineCorrectionButton.clicked.connect(self.baseline_correction)
        self.calculatePowerButton.clicked.connect(self.calculate_total_power)

        self.radioButton.toggled.connect(self.handle_radio_selection)
        self.radioButton_2.toggled.connect(self.handle_radio_selection)
        self.radioButton_3.toggled.connect(self.handle_radio_selection)
        self.radioButton_4.toggled.connect(self.handle_radio_selection)

        # Connecting buttons to their respective methods
        ##!!! implement the SingleWavelength option as well
        self.LoadLED.clicked.connect(lambda: self.load_file("LED Emission"))
        self.LoadTrans.clicked.connect(lambda: self.load_file("Epsilons A"))
        self.LoadCis.clicked.connect(lambda: self.load_file("Epsilons B"))
        self.LoadSpectra.clicked.connect(lambda: self.load_file("Spectral Data"))

        #self.label.setText("Total Power: -- mW ± -- mW")


        ############# EXPERIMENTAL PARAMETERS #############
        self.V = 3.0                 # Volume in ml
        self.k_BA = 7.240e-7         # Thermal back reaction rate s-1

        ######## POWER ########
        self.I0_avg = 743             # Photon flux in microWatt
        self.I0_err = 4                # Error on photon fplux in microWatt
        self.LEDw = 340

        # Set text in QPlainTextEdit using setPlainText
        self.plainTextEdit.setPlainText(str(self.V))
        self.plainTextEdit_2.setPlainText(str(self.k_BA))
        self.plainTextEdit_3.setPlainText(str(self.I0_avg))
        self.plainTextEdit_4.setPlainText(str(self.I0_err))
        self.plainTextEdit_5.setPlainText(str(self.LEDw))

        # Connect the textChanged signal to the update functions
        self.plainTextEdit.textChanged.connect(self.update_V)
        self.plainTextEdit_2.textChanged.connect(self.update_k_BA)
        self.plainTextEdit_3.textChanged.connect(self.update_I0_avg)
        self.plainTextEdit_4.textChanged.connect(self.update_I0_err)
        self.plainTextEdit_5.textChanged.connect(self.update_LEDw)

        # Adding new buttons for custom plot functions from your scripts
        self.plotLEDButton.clicked.connect(self.process_LED)
        self.plotEpsilonButton.clicked.connect(self.plot_epsilon)
        self.plotQuantButton.clicked.connect(self.plot_auto_quant)

        # For storing data and baseline corrected data
        self.data = None
        self.baseline_corrected_data = None
        self.filename_LED = None

        self.loadDataButton.setEnabled(False)
        self.baselineCorrectionButton.setEnabled(False)
        self.calculatePowerButton.setEnabled(False)
        self.LoadLED.setEnabled(False)
        #self.radioButton_5.setEnabled(False)
        #self.radioButton_6.setEnabled(False)

    # Update methods for the parameters
    def update_V(self):
        try:
            self.V = float(self.plainTextEdit.toPlainText())  # Convert the input to a float
            print(f"Updated V to {self.V}")
        except ValueError:
            pass  # Handle the case where the input is not a valid number

    def update_k_BA(self):
        try:
            self.k_BA = float(self.plainTextEdit_2.toPlainText())  # Convert the input to a float
            print(f"Updated k_BA to {self.k_BA}")
        except ValueError:
            pass

    def update_I0_avg(self):
        try:
            self.I0_avg = int(self.plainTextEdit_3.toPlainText())  # Convert the input to an integer
            print(f"Updated I0_avg to {self.I0_avg}")
        except ValueError:
            pass

    def update_I0_err(self):
        try:
            self.I0_err = int(self.plainTextEdit_4.toPlainText())  # Convert the input to an integer
            print(f"Updated I0_err to {self.I0_err}")
        except ValueError:
            pass

    def update_LEDw(self):
        try:
            self.LEDw = int(self.plainTextEdit_5.toPlainText())  # Convert the input to an integer
            print(f"Updated LEDw to {self.LEDw}")
        except ValueError:
            pass

    def handle_radio_selection(self):
        if self.radioButton_3.isChecked():
            # Enable TextLabel and disable Load button
            self.plainTextEdit_3.setEnabled(True)
            self.plainTextEdit_4.setEnabled(True)
            self.loadDataButton.setEnabled(False)
            self.loadDataButton.setEnabled(False)
            self.baselineCorrectionButton.setEnabled(False)
            self.calculatePowerButton.setEnabled(False)
        
        if self.radioButton_4.isChecked():
            # Disable TextLabel and enable Load button
            self.plainTextEdit_3.setEnabled(False)
            self.plainTextEdit_4.setEnabled(False)
            self.loadDataButton.setEnabled(True)
            self.baselineCorrectionButton.setEnabled(True)
            self.calculatePowerButton.setEnabled(True)
        
        if self.radioButton.isChecked():
            self.plainTextEdit_5.setEnabled(True)
            self.LoadLED.setEnabled(False)
            #self.radioButton_5.setEnabled(False)
            #self.radioButton_6.setEnabled(False)

        if self.radioButton_2.isChecked():
            self.plainTextEdit_5.setEnabled(False)
            self.LoadLED.setEnabled(True)
            #self.radioButton_5.setEnabled(True)
            #self.radioButton_6.setEnabled(True)

    def update_display(self, idx, line_index, new_x):
        self.line_positions[idx][line_index] = new_x 

    def load_power(self):
        """Load the data from a file and plot it on a new tab."""
        options = QtWidgets.QFileDialog.Options()
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                             "Load Data File", "",
                                                             "Data Files (*.csv);;All Files (*)",
                                                             options=options)
        if file_name:
            self.filename = file_name
            self.x, self.RefPower = self.load_and_validate_data(self.filename)
            # Save x and RefPower without overwriting
            self.loaded_data.append((self.x, self.RefPower))

            if self.RefPower is not None and self.x is not None:
                self.add_new_tab(self.plot_power_data, f"Loaded Data {self.count}", self.count)
            self.count += 1

    def load_and_validate_data(self, file_name):
        """Load data and handle errors."""
        x, RefPower = load_data(file_name)
        if RefPower is None or x is None:
            print("Failed to load data.")
        return x, RefPower

    def add_new_tab(self, plot_func, title, idx=None, *args):
        """Create a new tab with a plot and a navigation toolbar."""
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout()  # Use QVBoxLayout to stack toolbar and canvas vertically

        # Create the custom MplCanvas, passing idx
        canvas = MplCanvas(self, idx)  # Now idx will be passed correctly

        # Create the navigation toolbar for the canvas
        toolbar = NavigationToolbar(canvas, self)

        # Add the toolbar and canvas to the layout
        layout.addWidget(toolbar)  # Add the toolbar at the top
        layout.addWidget(canvas)   # Add the canvas (plot area) below the toolbar

        tab.setLayout(layout)

        # Call the plotting function to plot on the canvas
        if idx is not None:
            print('plot nr',idx)
            plot_func(canvas, idx, *args)  # Pass idx only if it's provided
            #self.tabWidget.tabCloseRequested.connect(self.close_tab,idx)
        else:
            print('plot nr None')
            plot_func(canvas, *args)  # Don't pass idx if it's None

        # Add the new tab to the tab widget
        self.tabWidget.addTab(tab, title)
        # Make tabs closable
        self.tabWidget.setTabsClosable(True)

        # Connect the tab close request to a custom method to handle closing
        self.tabWidget.tabCloseRequested.connect(self.close_tab)

    def close_tab(self, index):
        """Handle the tab close request."""
        # Get the tab to be closed
        tab = self.tabWidget.widget(index)

        # Handle any cleanup for idx and loaded data
        if hasattr(tab, 'idx') and tab.idx is not None:
            idx = tab.idx
            print(f"Closing tab with idx: {idx}")

            # Remove idx and loaded data
            if idx in self.idx:
                self.idx.remove(idx)  # Remove idx from self.idx list
                print(f"Removed idx: {idx} from self.idx")
            
            if idx in self.loaded_data:
                self.loaded_data[idx] = None# Remove the corresponding data from loaded_data
                self.line_positions[idx] = None
                print(f"Removed data for idx: {idx} from loaded_data")
        
        # Remove the tab
        self.tabWidget.removeTab(index)

    def plot_power_data(self, canvas, idx):
        """Plot the loaded data using the MplCanvas."""
        if self.x is None or self.RefPower is None:
            QtWidgets.QMessageBox.warning(self, "Error", "No data loaded")
            return
        
        self.line_positions.extend([[0] * 12 for _ in range(idx - len(self.line_positions) + 1)])
        # Use the plot_power method from MplCanvas
        canvas.plot_power(self.x, self.RefPower, self.filename)
        self.line_positions[idx] = self.get_sorted_line_positions(canvas)


    ############################
    #plot_sections *************
    ############################################################################################################

    def baseline_correction(self):
        """Apply baseline correction to each loaded dataset and create a new tab for each."""
        if not self.loaded_data:
            QtWidgets.QMessageBox.warning(self, "Error", "No data loaded")
            return

        # Iterate over each loaded dataset using the saved idx (self.count)
        for idx in range(len(self.loaded_data)):
            # Baseline correction: no jacket, no cuvette
            self.add_new_tab(self.plot_baseline_corrected_no, f"Baseline correction: no jacket, no cuvette {idx}", idx)
            
            # Baseline correction: jacket, cuvette with solvent
            self.add_new_tab(self.plot_baseline_corrected_yes, f"Baseline correction: jacket, cuvette with solvent {idx}", idx)

    def calculate_baseline_and_power_no(self, x, RefPower, sections):
        """Calculate baseline and baseline-corrected power."""
        # Baseline correction: no jacket, no cuvette
        x_masked = np.concatenate((x[sections["start_0"]:sections["end_0"]], x[sections["start_0_1"]:sections["end_0_1"]]))
        y_masked = np.concatenate((RefPower[sections["start_0"]:sections["end_0"]], RefPower[sections["start_0_1"]:sections["end_0_1"]]))
        x_masked = x_masked[~np.isnan(y_masked)]
        y_masked = y_masked[~np.isnan(y_masked)]

        # Polynomial fit (n = 3)
        n = 3
        p0 = np.full(n, 0.000000001)
        popt, _ = curve_fit(fit_func, x_masked, y_masked, p0=p0)
        baseline_1 = np.polyval(popt, x)
        baselined_1 = RefPower - baseline_1

        return baselined_1, baseline_1
    
    def calculate_baseline_and_power_yes(self, x, RefPower, sections):
        # Baseline correction: jacket, cuvette with solvent
        x_masked = np.concatenate((x[sections["start_2"]:sections["end_2"]], x[sections["start_4"]:sections["end_4"]]))
        y_masked = np.concatenate((RefPower[sections["start_2"]:sections["end_2"]], RefPower[sections["start_4"]:sections["end_4"]]))
        x_masked = x_masked[~np.isnan(y_masked)]
        y_masked = y_masked[~np.isnan(y_masked)]

        # Polynomial fit (n = 6)
        n = 6
        p0 = np.full(n, 0.000000001)
        popt, _ = curve_fit(fit_func, x_masked, y_masked, p0=p0)
        baseline_2 = np.polyval(popt, x)
        baselined_2 = RefPower - baseline_2

        return baselined_2, baseline_2

    def get_sorted_line_positions(self, canvas):
        """Get and sort the x-positions of vertical lines."""
        line_positions = [line.get_xdata()[0] for line in canvas.lines]
        line_positions.sort()
        return line_positions

    def plot_baseline_corrected_no(self, canvas, idx):
        case = 0 # no cuvette
        x, RefPower = self.loaded_data[idx]
        line_positions = self.line_positions[idx]
        #print(line_positions)
        #filename = f"Baseline Corrected {idx+1}"
        # Ensure 12 positions are selected
        if len(line_positions) != 12:
            QtWidgets.QMessageBox.warning(self, "Error", "You need to select exactly 12 line positions.")
            return
        #if idx == 1:
            # Proceed with baseline correction
        sections = choose_sections(line_positions)
        baselined, baseline = self.calculate_baseline_and_power_no(x, RefPower, sections)#baseline_correction(RefPower, x, sections)
        self.all_corrected_power.append(baselined)
        """Plot the baseline-corrected data."""
        canvas.plot_baseline_correction(x, RefPower, baseline, baselined, sections, case)


    def plot_baseline_corrected_yes(self, canvas, idx):
        case = 1 # yes cuvette
        x, RefPower = self.loaded_data[idx]
        line_positions = self.line_positions[idx]
        #print(line_positions)
        #filename = f"Baseline Corrected {idx+1}"
        # Ensure 12 positions are selected
        if len(line_positions) != 12:
            QtWidgets.QMessageBox.warning(self, "Error", "You need to select exactly 12 line positions.")
            return

        sections = choose_sections(line_positions)
        baselined, baseline = self.calculate_baseline_and_power_yes(x, RefPower, sections) #baseline_correction(RefPower, x, sections)
        self.all_corrected_power.append(baselined)
        """Plot the baseline-corrected data."""
        canvas.plot_baseline_correction(x, RefPower, baseline, baselined, sections, case)

    def calculate_total_power(self):
        """Calculate the total power and standard deviation from all baseline-corrected data."""
        if not self.all_corrected_power:
            QtWidgets.QMessageBox.warning(self, "Error", "No baseline correction has been applied yet.")
            #self.label.setText("Total Power: missing data")
            return

        means = []
        stds = []

        for idx, corrected_data in enumerate(self.all_corrected_power):
            corrected_data = np.asarray(corrected_data).ravel()
            if len(corrected_data) == 0:
                continue
            mean_power = np.nanmean(corrected_data)
            std_power = np.nanstd(corrected_data)
            means.append(mean_power)
            stds.append(std_power)
            print(f" Power: {mean_power:.7f} mW ± {std_power:.7f} mW")
        #if len(means) != 3 or len(stds) != 3:
        #    QtWidgets.QMessageBox.warning(self, "Error", "Missing some baseline corrections")
        #    return

        total_power = np.mean(means)
        total_variance = np.sum(np.square(stds)) / len(means)
        total_std = np.sqrt(total_variance)

        self.plainTextEdit_3.setPlainText(f"{total_power:.2f}")
        self.plainTextEdit_4.setPlainText(f"{total_std:.2f}")
        self.I0_avg = total_power
        self.I0_err = total_std
        #self.label.setText(f"Total Power: {total_power:.7f} mW ± {total_std:.7f} mW")
        #print(f"Total Power: {total_power:.7f} mW ± {total_std:.7f} mW")
    ####################################################################################################################################

    def load_file(self, file_type):
        """Load a file (LED emission, spectral data, epsilons) based on the file type."""
        options = QtWidgets.QFileDialog.Options()

        # File dialog for selecting files
        file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self, 
                                                             f"Load {file_type} File", "",
                                                             "CSV, DAT Files (*.csv *dat);;DAT Files (*.dat);;All Files (*)", 
                                                             options=options)
        if not file_name:
            QtWidgets.QMessageBox.warning(self, "Error", f"No {file_type} file selected")
            return

        ##################################################
        ## Move Load to .utils
        ##################################################
        # Store the file path in the appropriate attribute based on the file type
        file_ext = os.path.splitext(file_name)[1].lower()
        if file_ext == '.csv':
            self.load_csv(file_name, file_type)
        elif file_ext == '.dat':
            self.load_dat(file_name, file_type)
        else:
            QtWidgets.QMessageBox.warning(self, "Error", f"Unknown file type: {file_type}")

        QtWidgets.QMessageBox.information(self, "Success", f"{file_type} file loaded successfully!")

    def load_dat(self, file_path, file_type):
        """Load the data file depending on its format (.csv or .dat)."""
        if file_type == "LED Emission":
            self.emission_wavelengths, self.emission_Intensity = Integration.Import_LEDemission("Spectragryph", file_path)
            self.filename_LED = file_path
            threshold_LED = 500     #!!! HARDCODED IN THE WRONG PLACE ###############################
            self.emission_Intensity_proc, self.LEDindex_first, self.LEDindex_last, self.wavelength_low, self.wavelength_high = \
                Integration.Processing_LEDemission(
                    self.emission_wavelengths, self.emission_Intensity, threshold_LED)
        elif file_type == "Epsilons A":
            self.epsilon_A_wavelengths, self.epsilon_A_values = Integration.Import_Epsilons("Spectragryph", file_path)
        elif file_type == "Epsilons B":
            self.epsilon_B_wavelengths, self.epsilon_B_values = Integration.Import_Epsilons("Spectragryph", file_path)
        elif file_type == "Spectral Data":
            self.SpectralData_Wavelengths, self.SpectralData_Abs, self.SpectralData_Index = Integration.Import_SpectralData("Spectragryph", file_path, 
                                                                                                             self.wavelength_low, self.wavelength_high, self.LEDw) #HARDCODED IN THE WRONG PLACE

    def load_csv(self, file_path, file_type):
        """Load the data file depending on its format (.csv or .dat)."""
        
        if file_type == "LED Emission": ###### A LOT OF PROBLEMS WITH '
            self.emission_wavelengths, self.emission_Intensity = Integration.Import_LEDemission("Not", file_path)
            self.filename_LED = file_path
            threshold_LED = 500     #!!! HARDCODED IN THE WRONG PLACE ###############################
            self.emission_Intensity_proc, self.LEDindex_first, self.LEDindex_last, self.wavelength_low, self.wavelength_high = \
                Integration.Processing_LEDemission(
                    self.emission_wavelengths, self.emission_Intensity, threshold_LED)
        elif file_type == "Epsilons A":
            self.epsilon_A_wavelengths, self.epsilon_A_values = Integration.Import_Epsilons("Not", file_path)
        elif file_type == "Epsilons B":
            self.epsilon_B_wavelengths, self.epsilon_B_values = Integration.Import_Epsilons("Not", file_path)
        elif file_type == "Spectral Data":
            self.SpectralData_Wavelengths, self.SpectralData_Abs, self.SpectralData_Index = Integration.Import_SpectralData("Not", file_path, 
                                                                                                                        self.wavelength_low, self.wavelength_high, self.LEDw) #HARDCODED IN THE WRONG PLACE


    ##DEFINE OUTSIDE OF CLASS    
    def plot_LEDprocessed(self,canvas):
            
        canvas.Plot_LEDemission_Processed(self.SpectralData_Abs, self.SpectralData_Wavelengths,
            self.emission_wavelengths, self.emission_Intensity,
            self.LEDindex_first, self.LEDindex_last, 
            self.emission_Intensity_proc)
            
    def process_LED(self):
        """Process and visualize the data based on the loaded files."""
        if self.emission_wavelengths is None or self.emission_Intensity is None or self.SpectralData_Abs is None:
            QtWidgets.QMessageBox.warning(self, "Error", "Please load LED emission file and Spectra.")
            return
        self.add_new_tab(self.plot_LEDprocessed, "LED Emission and Spectral Data")
        
        
    ##DEFINE OUTSIDE OF CLASS    
    def plot_epsilon(self):
        if self.emission_Intensity_proc is None:
            QtWidgets.QMessageBox.warning(self, "Error", "Please process the LED emission data first.")
            return
        self.epsilon_A_interp, self.epsilon_B_interp, self.emission_interp = Integration.Interpolate_Epsilons(self.SpectralData_Wavelengths,
                     self.epsilon_A_wavelengths, self.epsilon_A_values,
                     self.epsilon_B_wavelengths, self.epsilon_B_values,
                     self.emission_wavelengths, self.emission_Intensity_proc)
        #print(self.epsilon_A_interp)
        #print(self.epsilon_B_interp)
        #print(self.emission_interp)
        # Example: Plot the data using MplCanvas
        def plot_func(canvas):
            canvas.plot_Epsilons(self.SpectralData_Wavelengths,
                          self.epsilon_A_interp, self.epsilon_B_interp,
                          self.emission_interp)
        
        self.add_new_tab(plot_func, "Epsilons")

    def plot_auto_quant(self, canvas):
        """Extract quantum yields, calculate concentrations, and plot/save results."""
        #if not hasattr(self, 'fit_results') or not hasattr(self, 'N'):
        #    QtWidgets.QMessageBox.warning(self, "Error", "Please complete the quantum yield minimization first.")
        #    return

        # Create parameters needed for fitting
        initial_conc_A, initial_conc_B, total_absorbance, lambda_meters, normalized_emission = \
            Integration.CreateParameters(self.SpectralData_Abs, self.SpectralData_Wavelengths,
                                        self.epsilon_A_interp, self.emission_interp)

        log_file = self.filename_LED
        I0_list = self.GetPowerList(self.I0_avg, self.I0_err)   
        self.timestamps = self.GetTimestamps(log_file)
        
        N, fit_results = Integration.MinimizeQYs(I0_list, normalized_emission,
                                                lambda_meters, 
                                                initial_conc_A, initial_conc_B,
                                                self.timestamps, self.SpectralData_Abs,
                                                self.epsilon_A_interp, self.epsilon_B_interp,
                                                self.V)

        # Extract results from the fit
        self.QY_AB_opt, self.QY_BA_opt, self.error_QY_AB, self.error_QY_BA = Integration.ExtractResults(fit_results)

        # Calculate optimized concentrations
        self.conc_opt = Integration.CalculateConcentrations(lambda_meters,
                                                    initial_conc_A, initial_conc_B, 
                                                    self.timestamps,
                                                    self.QY_AB_opt, self.QY_BA_opt, 
                                                    self.epsilon_A_interp, self.epsilon_B_interp,
                                                    N, self.V)

        # Calculate total absorbance and residuals
        self.total_abs_fit, self.residuals = Integration.GetFittedAbs(fit_results, self.conc_opt,
                                                            self.epsilon_A_interp, self.epsilon_B_interp,
                                                            self.timestamps,
                                                            self.SpectralData_Wavelengths)
        # print('##########################',self.timestamps)
        
        # Plot and save the results
        self.add_new_tab(self.Plot_QY, "QY")
        QtWidgets.QMessageBox.information(self, "Success", "Results extracted, plotted, and saved!")

    def Plot_QY(self, canvas):
        canvas.PlotAndSave(self.LEDw,
                           self.timestamps,
                           self.conc_opt,
                           self.SpectralData_Abs,
                           self.SpectralData_Index,
                           self.total_abs_fit,
                           self.residuals,
                           self.QY_AB_opt, self.QY_BA_opt,
                           self.error_QY_AB, self.error_QY_BA)

    def GetPowerList(self, I0_avg, I0_err): ###??????????
        I0_list = [I0_avg, I0_avg+I0_err, I0_avg-I0_err]
        return I0_list

    def GetTimestamps(self, LogFile):
        LogFile = r'C:\Users\jorst136\Documents\Postdoc\Projects\Setup_QY\GUI-QY\SteGUI_v0.5\ExampleData\Azobenzene_340nm_log.csv'
        #LogFile = 'c:/Users/Hallina 9000/Desktop/SteGUI/ExampleData_autoQuant/Azobenzene_340nm_log.csv' 
        #CSV not DAT
        log = pd.read_csv(LogFile,
                        sep = ",", decimal = ".", skiprows = 1, header=None,)
        log_t=log[log[3] == 'Measure']
        # print("#################", log_t)
        t=log_t[2]
        timestamps=t.to_numpy()
        return timestamps


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = PowerProcessingApp()
    window.show()
    sys.exit(app.exec_())
