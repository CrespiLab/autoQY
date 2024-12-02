"""
@author: Alfredo and Jorn
=============================================================================
autoQuant

=============================================================================
The PowerProcessing part is used for processing of the data measured with the 
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

=============================================================================

TO DO:
[DONE] add feature: pick two or more input files; 
   then perform averaging and error calculation to yield final result
[] Output a Results file containing all the input values (incl. Exp. Param.)
[] Add feature to input starting concentrations

GUI features to add:
[DONE] Pop-up window to select file
[DONE] Show percentages at PSS
=============================================================================

"""
import sys, os
import numpy as np
import pandas as pd
# from scipy.optimize import curve_fit #change in baseline correction

from PyQt5 import QtWidgets, uic
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

import autoQuant.Integration as Integration
import autoQuant.ExpParams as ExpParams
import autoQuant.LoadedData as LoadedData
import autoQuant.SingleWavelength as SingleWavelength

import tools.load_data as LoadData
from tools.plotting import MplCanvas
import tools.baseline_power as BaselinePower
import tools.extractresults as ExtractResults
#from tools.style import apply_dark_theme

import tools.WindowPowerProcessing as WindowPowerProcessing

class PowerProcessingApp(QtWidgets.QMainWindow):
        ############################
        #Error Handling  ***********
        ############################
    def __init__(self):
        super(PowerProcessingApp, self).__init__()
        uic.loadUi('UIs/qt.ui', self)  # Load the UI file you provided

        ############################
        #Change style  *************
        ############################
        #self.setStyleSheet(get_stylesheet())

        ##!!! REMOVE UNNECESSARY VARIABLES HERE

        self.led_file = None
        self.spectral_data_path = './'
        self.spectral_data_file = 'output_AutoQuant'
        self.eps_file_a = None
        self.eps_file_b = None

        self.all_corrected_power = []

        self.filename_LED = None
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
        
        self.CalculationMethod = "Integration" # default Calculation Method

        ## Button connections
        self.loadDataButton.clicked.connect(self.OpenWindow_PowerProcessing_1)
        self.loadDataButton_2.clicked.connect(self.OpenWindow_PowerProcessing_2)
        self.loadDataButton_3.clicked.connect(self.OpenWindow_PowerProcessing_3)

        ##!!! FINISH FUNCTION        
        self.calculatePowerButton.clicked.connect(self.calculate_total_power) ## Calculates average power+error

        #######

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
        self.LoadLog.clicked.connect(lambda: self.load_file("Log Irr"))

        self.SaveResultsBtn.clicked.connect(self.Save_QY)

        # Set text in QPlainTextEdit using setPlainText
        self.plainTextEdit.setPlainText(str(ExpParams.V)) # volume
        self.plainTextEdit_2.setPlainText(str(ExpParams.k_BA)) # rate constant
        self.plainTextEdit_3.setPlainText(str(ExpParams.I0_avg)) # Power Manual
        self.plainTextEdit_4.setPlainText(str(ExpParams.I0_err)) # Power Manual
        self.plainTextEdit_5.setPlainText(str(ExpParams.LEDw)) # Integration Mode (default)

        self.plainTextEdit_6.setPlainText(str(ExpParams.I0_avg)) # PowerProcessing: Calculated Power
        self.plainTextEdit_8.setPlainText(str(ExpParams.I0_err)) # PowerProcessing: Error
        self.plainTextEdit_9.setPlainText(str(ExpParams.threshold)) # Threshold for LED Emission spectrum

        ## Connect the textChanged signal to the update functions
        self.plainTextEdit.textChanged.connect(self.update_V)
        self.plainTextEdit_2.textChanged.connect(self.update_k_BA)
        self.plainTextEdit_3.textChanged.connect(self.update_I0_avg)
        self.plainTextEdit_4.textChanged.connect(self.update_I0_err)
        self.plainTextEdit_5.textChanged.connect(self.update_LEDw_Integration) # Integration Mode
        self.plainTextEdit_7.textChanged.connect(self.update_LEDw_SingleWl) # SingleWavelength Mode

        self.plainTextEdit_epsilonR.textChanged.connect(self.update_epsilon_R) # Epsilon Reactant (SingleWavelength)
        self.plainTextEdit_epsilonP.textChanged.connect(self.update_epsilon_P) # Epsilon Product (SingleWavelength)

        self.plainTextEdit_9.textChanged.connect(self.update_threshold) # 

        ## Adding new buttons for custom plot functions from your scripts
        self.plotLEDButton.clicked.connect(self.process_LED)
        self.plotEpsilonButton.clicked.connect(self.plot_epsilon)
        self.plotQuantButton.clicked.connect(self.plot_auto_quant)

        self.radioButton_3.setEnabled(True) # Power Manual: default enabled
        self.loadDataButton.setEnabled(False) # default Manual Power input
        self.baselineCorrectionButton.setEnabled(False) # default Manual Power input
        self.calculatePowerButton.setEnabled(False) # default Manual Power input
        
        self.radioButton_2.setEnabled(True) # Irradiation Integration: default 
        self.LoadLED.setEnabled(True) # default Integration Mode
        
        #### INITIALISATION ####
        self.handle_radio_selection() 

    ## Update methods for the parameters
    def update_V(self):
        try:
            ExpParams.V = float(self.plainTextEdit.toPlainText())  # Convert the input to a float
            print(f"Updated V to {ExpParams.V}")
        except ValueError:
            pass  # Handle the case where the input is not a valid number

    def update_k_BA(self):
        try:
            ExpParams.k_BA = float(self.plainTextEdit_2.toPlainText())  # Convert the input to a float
            print(f"Updated k_BA to {ExpParams.k_BA}")
        except ValueError:
            pass

    def update_I0_avg(self):
        try:
            ExpParams.I0_avg = float(self.plainTextEdit_3.toPlainText())  # Convert the input to an integer
            print(f"Updated I0_avg to {ExpParams.I0_avg}")
        except ValueError:
            pass

    def update_I0_err(self):
        try:
            ExpParams.I0_err = float(self.plainTextEdit_4.toPlainText())  # Convert the input to an integer
            print(f"Updated I0_err to {ExpParams.I0_err}")
        except ValueError:
            pass

    def update_I0_avg_PP(self):
        try:
            ExpParams.I0_avg = float(self.plainTextEdit_6.toPlainText())  # Convert the input to an integer
            print(f"Updated I0_avg to {ExpParams.I0_avg}")
        except ValueError:
            pass

    def update_I0_err_PP(self):
        try:
            ExpParams.I0_err = float(self.plainTextEdit_8.toPlainText())  # Convert the input to an integer
            print(f"Updated I0_err to {ExpParams.I0_err}")
        except ValueError:
            pass

    def update_LEDw_SingleWl(self):
        """ For SingleWavelength Mode """
        try:
            ExpParams.LEDw = int(self.plainTextEdit_7.toPlainText())  # Convert the input to an integer
            print(f"Updated LEDw to {ExpParams.LEDw}")
        except ValueError:
            pass

    def update_epsilon_R(self):
        """ Epsilon of Reactant for SingleWavelength Mode """
        try:
            ExpParams.epsilon_R = float(self.plainTextEdit_epsilonR.toPlainText())  # Convert the input to a float
            print(f"Updated epsilon_R to {ExpParams.epsilon_R}")
        except ValueError:
            pass

    def update_epsilon_P(self):
        """ Epsilon of Product for SingleWavelength Mode """
        try:
            ExpParams.epsilon_P = float(self.plainTextEdit_epsilonP.toPlainText())  # Convert the input to a float
            print(f"Updated epsilon_P to {ExpParams.epsilon_P}")
        except ValueError:
            pass

    def update_LEDw_Integration(self):
        """ For Integration Mode (default) """
        try:
            ExpParams.LEDw = int(self.plainTextEdit_5.toPlainText())  # Convert the input to an integer
            print(f"Updated LEDw to {ExpParams.LEDw}")
        except ValueError:
            pass

    def update_threshold(self):
        """ Update value for threshold used for LED emission spectrum """
        try:
            ExpParams.threshold = int(self.plainTextEdit_9.toPlainText())  # Convert the input to an integer
            print(f"Updated threshold to {ExpParams.threshold}")
        except ValueError:
            pass

    def handle_radio_selection(self):
        if self.radioButton_3.isChecked(): # Power Manual Input
            self.update_I0_avg # set I0_avg to current text
            self.update_I0_err # set I0_err to current text
            
            ## Enable TextLabel and disable Load button
            self.plainTextEdit_3.setEnabled(True)
            self.plainTextEdit_4.setEnabled(True)
            self.loadDataButton.setEnabled(False)
            self.baselineCorrectionButton.setEnabled(False)
            self.plainTextEdit_Power_1.setEnabled(False)
            self.plainTextEdit_PowerError_1.setEnabled(False)
            self.loadDataButton_2.setEnabled(False)
            self.baselineCorrectionButton_2.setEnabled(False)
            self.plainTextEdit_Power_2.setEnabled(False)
            self.plainTextEdit_PowerError_2.setEnabled(False)
            self.loadDataButton_3.setEnabled(False)
            self.baselineCorrectionButton_3.setEnabled(False)
            self.plainTextEdit_Power_3.setEnabled(False)
            self.plainTextEdit_PowerError_3.setEnabled(False)
            self.calculatePowerButton.setEnabled(False)
            self.plainTextEdit_6.setEnabled(False)
            self.plainTextEdit_8.setEnabled(False)
        
        if self.radioButton_4.isChecked(): # PowerProcessing
            self.update_I0_avg_PP # set I0_avg to current text in PowerProcessing field
            self.update_I0_err_PP # set I0_err to current text in PowerProcessing field
            
            ## Disable TextLabel and enable Load button
            self.plainTextEdit_3.setEnabled(False) # turn off Manual Input: Power
            self.plainTextEdit_4.setEnabled(False) # turn off Manual Input: Error
            self.loadDataButton.setEnabled(True)
            self.baselineCorrectionButton.setEnabled(True)
            self.plainTextEdit_Power_1.setEnabled(True)
            self.plainTextEdit_PowerError_1.setEnabled(True)
            self.loadDataButton_2.setEnabled(True)
            self.baselineCorrectionButton_2.setEnabled(True)
            self.plainTextEdit_Power_2.setEnabled(True)
            self.plainTextEdit_PowerError_2.setEnabled(True)
            self.loadDataButton_3.setEnabled(True)
            self.baselineCorrectionButton_3.setEnabled(True)
            self.plainTextEdit_Power_3.setEnabled(True)
            self.plainTextEdit_PowerError_3.setEnabled(True)
            self.calculatePowerButton.setEnabled(True)
            self.plainTextEdit_6.setEnabled(True)
            self.plainTextEdit_8.setEnabled(True)
        
        if self.radioButton.isChecked(): # LED SingleWavelength Mode
            self.update_LEDw_SingleWl # set variable to current text
            self.update_epsilon_R
            self.update_epsilon_P
            self.plainTextEdit_7.setEnabled(True) # SingleWavelength wavelength (nm)
            self.plainTextEdit_epsilonR.setEnabled(True) # SingleWavelength epsilon Reactant
            self.plainTextEdit_epsilonP.setEnabled(True) # SingleWavelength epsilon Product
            self.plainTextEdit_5.setEnabled(False) # Integration wavelength (nm)
            self.LoadLED.setEnabled(False)
            self.CalculationMethod = "SingleWavelength"

        if self.radioButton_2.isChecked(): # LED Integration Mode
            self.update_LEDw_Integration # set variable to current text
            self.plainTextEdit_7.setEnabled(False) # SingleWavelength wavelength (nm)
            self.plainTextEdit_epsilonR.setEnabled(False) # SingleWavelength epsilon Reactant
            self.plainTextEdit_epsilonP.setEnabled(False) # SingleWavelength epsilon Product
            self.plainTextEdit_5.setEnabled(True) # Integration wavelength (nm)
            self.LoadLED.setEnabled(True)
            self.CalculationMethod = "Integration"

    def OpenWindow_PowerProcessing_1(self):
        """Load the power data from a file and plot it in a new window."""
        LoadedData.count = 1
        # self.window_PP = WindowPowerProcessing.WindowPowerProcessing() # load Class that includes loadUi
        self.window_PP = WindowPowerProcessing.WindowPowerProcessing(parent=self) # load Class that includes loadUi
        
        self.window_PP.show()
        # self.plainTextEdit_Power_1.setPlainText(LoadedData.PowersAtCuvette[LoadedData.count]) # calculated power


    def OpenWindow_PowerProcessing_2(self):
        """Load the power data from a file and plot it in a new window."""
        LoadedData.count = 2
        self.window_PP = WindowPowerProcessing.WindowPowerProcessing(parent=self) # load Class that includes loadUi
        self.window_PP.show()

    def OpenWindow_PowerProcessing_3(self):
        """Load the power data from a file and plot it in a new window."""
        LoadedData.count = 3
        self.window_PP = WindowPowerProcessing.WindowPowerProcessing(parent=self) # load Class that includes loadUi
        self.window_PP.show()

    def add_new_tab(self, plot_func, title):  # , idx=None, *args
        """Create a new tab with a plot and a navigation toolbar."""
        try:
            tab = QtWidgets.QWidget()
            layout = QtWidgets.QVBoxLayout()  # Use QVBoxLayout to stack toolbar and canvas vertically

            ## Create the custom MplCanvas
            canvas = MplCanvas(self)  # idx not required in this implementation

            ## Create the navigation toolbar for the canvas
            toolbar = NavigationToolbar(canvas, self)

            ## Add the toolbar and canvas to the layout
            layout.addWidget(toolbar)  # Add the toolbar at the top
            layout.addWidget(canvas)   # Add the canvas (plot area) below the toolbar

            tab.setLayout(layout)

            ## Call the plotting function to populate the canvas
            plot_func(canvas) # idx functionality removed as it's unused in this version

            ## Add the new tab to the tab widget
            self.tabWidget.addTab(tab, title)

            ## Make tabs closable
            self.tabWidget.setTabsClosable(True)

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to add new tab: {e}")

    ############################################################################################################
    ############################ plot_sections ############################
    ############################################################################################################

    def calculate_total_power(self):
        """Calculate the total power and standard deviation from all baseline-corrected data."""
        ##!!! ADJUST THIS FUNCTION TO CALCULATE AVERAGE OF POWER DATA LOADED FROM 1-3 DATASETS
        
        if not self.all_corrected_power:
            QtWidgets.QMessageBox.warning(self, "Error", "No baseline correction has been applied yet.")
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

        #if len(means) != 3 or len(stds) != 3:
        #    QtWidgets.QMessageBox.warning(self, "Error", "Missing some baseline corrections")
        #    return

        total_power = np.mean(means)
        total_variance = np.sum(np.square(stds)) / len(means)
        total_std = np.sqrt(total_variance)

        self.plainTextEdit_6.setPlainText(f"{total_power:.2f}") # PowerProcessing: Power
        self.plainTextEdit_8.setPlainText(f"{total_std:.2f}") # PowerProcessing: Error
        ExpParams.I0_avg = total_power # set to calculated power
        ExpParams.I0_err = total_std # set to calculated error
        
        
        
        print(f"I0_avg: {ExpParams.I0_avg}\nI0_err: {ExpParams.I0_err}")
        
        ##!!! MAKE I_O LIST in ExpParams
        
    ####################################################################################################################################

    def load_file(self, file_type):
        """Load a file (LED emission, spectral data, epsilons) based on the file type."""
        try:
            options = QtWidgets.QFileDialog.Options()

            ## File dialog for selecting files
            file_name, _ = QtWidgets.QFileDialog.getOpenFileName(self, 
                                                                f"Load {file_type} File", "",
                                                                "CSV, DAT Files (*.csv *dat);;DAT Files (*.dat);;All Files (*)", 
                                                                options=options)
            if not file_name:
                QtWidgets.QMessageBox.warning(self, "Error", f"No {file_type} file selected")
                return

            ##################################################
            ##!!! Move Load to .utils
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
        except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to load {file_type} file: {e}")

    def load_dat(self, file_path, file_type):
        """Load the data file depending on its format (.csv or .dat)."""
        try:
            if file_type == "LED Emission":
                self.emission_wavelengths, self.emission_Intensity = Integration.Import_LEDemission("Spectragryph", file_path)
                self.filename_LED = file_path
                
            elif file_type == "Epsilons A":
                self.epsilon_A_wavelengths, self.epsilon_A_values = Integration.Import_Epsilons("Spectragryph", file_path)
            elif file_type == "Epsilons B":
                self.epsilon_B_wavelengths, self.epsilon_B_values = Integration.Import_Epsilons("Spectragryph", file_path)
            elif file_type == "Spectral Data":
                LoadedData.SpectralData_Full = \
                    Integration.Import_SpectralData("Spectragryph",file_path) #HARDCODED IN THE WRONG PLACE # STILL??
            ########################################
            ##!!! ADD in case of .dat format
            # elif file_type == "Log Irr":
            #      self.timestamps = self.GetTimestamps(file_path)
            else:
                    raise ValueError(f"Unknown file type: {file_type}")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to process .dat file for {file_type}: {e}")

    def load_csv(self, file_path, file_type):
        """Load the data file depending on its format (.csv or .dat)."""
        try:       
            if file_type == "LED Emission": ###### A LOT OF PROBLEMS WITH '
                self.emission_wavelengths, self.emission_Intensity = Integration.Import_LEDemission("Not", file_path)
                self.filename_LED = file_path
            elif file_type == "Epsilons A":
                self.epsilon_A_wavelengths, self.epsilon_A_values = Integration.Import_Epsilons("Not", file_path)
            elif file_type == "Epsilons B":
                self.epsilon_B_wavelengths, self.epsilon_B_values = Integration.Import_Epsilons("Not", file_path)
            elif file_type == "Spectral Data":
                LoadedData.SpectralData_Full = \
                    Integration.Import_SpectralData("Not", file_path)
            elif file_type == "Log Irr":
                self.timestamps = self.GetTimestamps(file_path)
            else:
                    raise ValueError(f"Unknown file type: {file_type}")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to process .csv file for {file_type}: {e}")

    ##!!! DEFINE OUTSIDE OF CLASS
    def plot_LEDprocessed(self,canvas):
            
        canvas.Plot_LEDemission_Processed(self.SpectralData_Abs, self.SpectralData_Wavelengths,
            self.emission_wavelengths, self.emission_Intensity,
            self.LEDindex_first, self.LEDindex_last, 
            self.emission_Intensity_proc)
            
    def process_LED(self):
        """Process and visualize the data based on the loaded files."""
        
        ## Integration mode
        if self.CalculationMethod == "Integration":
            if self.emission_wavelengths is None or self.emission_Intensity is None or LoadedData.SpectralData_Full is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load LED emission file and Spectra.")
                return
            
            if ExpParams.LEDw is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please set Wavelength (nm).")
                return
            
            ########################################
            threshold_LED = ExpParams.threshold
            self.emission_Intensity_proc, self.LEDindex_first, self.LEDindex_last, self.wavelength_low, self.wavelength_high = \
                Integration.Processing_LEDemission(
                    self.emission_wavelengths, self.emission_Intensity, threshold_LED)
            
            ########################################
            ## Process Spectral data: cut to part of spectrum according to LED emission band
            self.SpectralData_Wavelengths, self.SpectralData_Abs, self.SpectralData_Index = \
                Integration.Process_SpectralData(LoadedData.SpectralData_Full, self.wavelength_low, self.wavelength_high, ExpParams.LEDw)
        
        ##!!! ADJUST CODE TO PLOT SPECTRA AND JUST INDICATE PART OF SPECTRUM (with a box)
            self.add_new_tab(self.plot_LEDprocessed, "LED Emission and Spectral Data")
        
        ## SingleWavelength mode
        elif self.CalculationMethod == "SingleWavelength":
            if LoadedData.SpectralData_Full is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Spectra.")
                return


        ##!!! WORKING ON IT HERE NOT FINISHED YET AAAHHH
            self.SpectralData_Wavelengths, self.SpectralData_Abs, self.SpectralData_Index = \
                SingleWavelength.Process_SpectralData(LoadedData.SpectralData_Full, self.wavelength_low, self.wavelength_high, ExpParams.LEDw)

            ##!!! ADJUST CODE TO PLOT SPECTRA AND ONLY INDICATE EITHER
                ## VERTICAL LINE: SINGLE WAVELENGTH
            # self.add_new_tab(self.plot_LEDprocessed, "LED Emission and Spectral Data")
        
        else:
            QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the self.CalculationMethod variable")
        
        #############
        
    ##!!! DEFINE OUTSIDE OF CLASS    
    def plot_epsilon(self):
        """ Interpolate epsilons spectra and plot """
        if self.emission_Intensity_proc is None:
            QtWidgets.QMessageBox.warning(self, "Error", "Please process the LED emission data first.")
            return
        self.epsilon_A_interp, self.epsilon_B_interp, self.emission_interp = Integration.Interpolate_Epsilons(self.SpectralData_Wavelengths,
                     self.epsilon_A_wavelengths, self.epsilon_A_values,
                     self.epsilon_B_wavelengths, self.epsilon_B_values,
                     self.emission_wavelengths, self.emission_Intensity_proc)

        def plot_func(canvas):
            """ Plot the data using MplCanvas """
            canvas.plot_Epsilons(self.SpectralData_Wavelengths,
                          self.epsilon_A_interp, self.epsilon_B_interp,
                          self.emission_interp)
        
        self.add_new_tab(plot_func, "Epsilons")

    def plot_auto_quant(self, canvas):
        """
        Calculate quantum yields by numerically solving the differential equations.
        Then calculate the concentrations, and plot the results.
        """
        I0_list = self.GetPowerList(ExpParams.I0_avg, ExpParams.I0_err) # make a list of the three values for I0

        if self.CalculationMethod == "Integration":
            ## Create parameters needed for fitting
            initial_conc_A, initial_conc_B, spectraldata_meters, normalized_emission = \
                Integration.CreateParameters(self.SpectralData_Abs, self.SpectralData_Wavelengths,
                                            self.epsilon_A_interp, self.emission_interp)
    
            
            N, fit_results = Integration.MinimizeQYs(I0_list, normalized_emission,
                                                    spectraldata_meters, 
                                                    initial_conc_A, initial_conc_B,
                                                    self.timestamps, self.SpectralData_Abs,
                                                    self.epsilon_A_interp, self.epsilon_B_interp,
                                                    ExpParams.V)

        elif self.CalculationMethod == "SingleWavelength":
            ##!!! WORKING ON THIS

            print(f"plot_auto_quant SingleWavelength SpectralData_Abs:\n{self.SpectralData_Abs}")
            
            ## Create parameters needed for fitting
            initial_conc_A, initial_conc_B, self.SpectralData_AbsAtLEDw = \
                SingleWavelength.CreateParameters(self.SpectralData_Abs, ExpParams.LEDw, ExpParams.epsilon_R)
            
            
            N, fit_results = SingleWavelength.MinimizeQYs(I0_list, ExpParams.LEDw,
                                                          initial_conc_A, initial_conc_B,
                                                          self.timestamps, self.SpectralData_AbsAtLEDw,
                                                          ExpParams.epsilon_R, ExpParams.epsilon_P,
                                                          ExpParams.V)
            
        else:
            QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the self.CalculationMethod variable")
            


        ######################################################################
        ##!!! MOVE TO SEPARATE PY SCRIPT (TOOLS OR SOMETHING)
        ## Extract results from the fit
        
        self.QY_AB_opt, self.QY_BA_opt, self.error_QY_AB, self.error_QY_BA = ExtractResults.ExtractResults(fit_results)

        ## Calculate optimized concentrations
        self.conc_opt, self.PSS_Reactant, self.PSS_Product = Integration.CalculateConcentrations(spectraldata_meters,
                                                    initial_conc_A, initial_conc_B, 
                                                    self.timestamps,
                                                    self.QY_AB_opt, self.QY_BA_opt, 
                                                    self.epsilon_A_interp, self.epsilon_B_interp,
                                                    N, ExpParams.V)

        ## Calculate total absorbance and residuals
        self.total_abs_fit, self.residuals = Integration.GetFittedAbs(fit_results, self.conc_opt,
                                                            self.epsilon_A_interp, self.epsilon_B_interp,
                                                            self.timestamps,
                                                            self.SpectralData_Wavelengths)
        
        ######################################################################
        
        ## Update labels
        ##!!! MAKE THIS A SEPARATE FUNCTION
        self.textEdit_QY_RtoP.setText(f"{self.QY_AB_opt:.3f}") # optimised QY R to P
        self.textEdit_QY_PtoR.setText(f"{self.QY_BA_opt:.3f}") # optimised QY P to R
        
        self.textEdit_QYerror_RtoP.setText(f"{self.error_QY_AB:.3f}") # error R to P
        self.textEdit_QYerror_PtoR.setText(f"{self.error_QY_BA:.3f}") # error P to R
        
        self.textEdit_PSS_R.setText(f"{self.PSS_Reactant:.1f}") # %R at PSS
        self.textEdit_PSS_P.setText(f"{self.PSS_Product:.1f}") # %P at PSS
        
        ## Plot and save the results
        self.add_new_tab(self.Plot_QY, "QY")
        QtWidgets.QMessageBox.information(self, "Success", "Results extracted, plotted, and saved!")

    def Plot_QY(self, canvas):
        canvas.PlotResults(ExpParams.LEDw,
                           self.timestamps,
                           self.conc_opt,
                           self.SpectralData_Abs,
                           self.SpectralData_Index,
                           self.total_abs_fit,
                           self.residuals,
                           self.QY_AB_opt, self.QY_BA_opt,
                           self.error_QY_AB, self.error_QY_BA,
                           self.CalculationMethod)

    def Save_QY(self):
        """ Save results """
        options = QtWidgets.QFileDialog.Options()

        ## File dialog for selecting files
        savefilename, _ = QtWidgets.QFileDialog.getSaveFileName(self, 
                                                             "Choose name of savefile", "",
                                                             options=options)
        
        canvas = MplCanvas(self)
        
        canvas.PlotResults(ExpParams.LEDw,
                           self.timestamps,
                           self.conc_opt,
                           self.SpectralData_Abs,
                           self.SpectralData_Index,
                           self.total_abs_fit,
                           self.residuals,
                           self.QY_AB_opt, self.QY_BA_opt,
                           self.error_QY_AB, self.error_QY_BA,
                           self.CalculationMethod,
                           SaveResults = "Yes",
                           SaveFileName = savefilename)

        ##!!! ADD: output a "Results" file containing all the used parameters and the obtained results

        QtWidgets.QMessageBox.information(self, "Success", f"{savefilename} file saved successfully!")


    def GetPowerList(self, I0_avg, I0_err): ###??????????
        I0_list = [I0_avg, I0_avg+I0_err, I0_avg-I0_err]
        return I0_list

    def GetTimestamps(self, LogFile):
        ## CSV not DAT
        log = pd.read_csv(LogFile,
                        sep = ",", decimal = ".", skiprows = 1, header=None,)
        log_t=log[log[3] == 'Measure']
        t=log_t[2]
        timestamps=t.to_numpy()
        return timestamps

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = PowerProcessingApp()
    window.show()
    sys.exit(app.exec_())
