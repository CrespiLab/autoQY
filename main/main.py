"""
Authors: Alfredo and Jorn
Date started: Sep 16th 2024
Based on: script by Anouk Volker from our publication in Beilstein Journal of Organic Chemistry 2024

Programme with GUI to calculate the QYs from the absorption data and the 
emission spectrum of the LED
=============================================================================
autoQY
=============================================================================

"""
import sys, os
import copy
import numpy as np


from PyQt5 import QtWidgets
# from PyQt5 import uic
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

import QY.integration as Integration
import QY.fractions as Fractions
import data.experimental_parameters as ExpParams
import data.loaded_data as LoadedData
import data.results as Results
import data.calc_settings as CalcSettings
import data.datasets as Datasets
import user_config.defaults as Defaults

import tools.load_data as LoadData
from tools.plotting import MplCanvas
import tools.extractresults as ExtractResults
#from tools.style import apply_dark_theme

import tools.power_processing as PowerProcessing

from UIs.MainWindow import Ui_MainWindow

##!!! ADD statusBar messages: self.statusBar.showMessage(f"success, successful: {err}")    
    ## for successful loading
## make it a statusFIELD so that there is a log of messages visible
## for severe issues: use pop-up window still

def main():
    class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
            ############################
            #Error Handling  ***********
            ############################
        def __init__(self):
            super(MainWindow, self).__init__()
            # uic.loadUi('UIs/MainWindow-large.ui', self)  # Old way: load the UI file
            self.setupUi(self)
            self.SetDefaultSettings() ## set default settings from defaults.py (editable by user)

            ## Save a copy of original state            
            ## Automatically capture only user-defined variables (not built-ins)
            self.default_state = {
                ExpParams: {
                    name: copy.deepcopy(value)
                    for name, value in vars(ExpParams).items()
                    if not name.startswith('__') and not callable(value)
                },
                CalcSettings: {
                    name: copy.deepcopy(value)
                    for name, value in vars(CalcSettings).items()
                    if not name.startswith('__') and not callable(value)
                }
            }

            ############################
            #Change style  *************
            ############################
            #self.setStyleSheet(get_stylesheet())
    
            self.labels_power = {1: self.lineEdit_Power_1,
                                       2: self.lineEdit_Power_2,
                                       3: self.lineEdit_Power_3}
            self.labels_error = {1: self.lineEdit_PowerError_1,
                                       2: self.lineEdit_PowerError_2,
                                       3: self.lineEdit_PowerError_3}
            
            ## Default Calculation Method is integration: see ExpParams
    
            ## Button connections
            self.loadDataButton_1.clicked.connect(lambda: self.OpenWindow_PowerProcessing(1))
            self.loadDataButton_2.clicked.connect(lambda: self.OpenWindow_PowerProcessing(2))
            self.loadDataButton_3.clicked.connect(lambda: self.OpenWindow_PowerProcessing(3))
            
            self.ButtonClearLoadedData.clicked.connect(self.ClearLoadedData)
            
            self.DeletePowerDataButton_1.clicked.connect(lambda: self.DeletePowerData(1))
            self.DeletePowerDataButton_2.clicked.connect(lambda: self.DeletePowerData(2))
            self.DeletePowerDataButton_3.clicked.connect(lambda: self.DeletePowerData(3))
    
            self.calculatePowerButton.clicked.connect(self.calculate_total_power) ## Calculates average power+error
    
            #######
            self.radioButton_ODE_Concentration.toggled.connect(self.handle_radio_selection_ODE)
            self.radioButton_ODE_Emission.toggled.connect(self.handle_radio_selection_ODE)
            self.radioButton_PowerManual.toggled.connect(self.handle_radio_selection_Power)
            self.radioButton_PowerProcessing.toggled.connect(self.handle_radio_selection_Power)
            self.radioButton_Log_Default.toggled.connect(self.handle_radio_selection_Log) # Timestamps Default
            self.radioButton_Log_AHK.toggled.connect(self.handle_radio_selection_Log) # Timestamps AHK
            self.radioButton_blcorrLED_on.toggled.connect(self.handle_radio_selection_blcorrLED)
            self.radioButton_blcorrLED_off.toggled.connect(self.handle_radio_selection_blcorrLED)
    
            # Connecting buttons to their respective methods
            self.LoadLED.clicked.connect(lambda: self.load_file("LED Emission"))
            self.LoadEpsilons_Reactant.clicked.connect(lambda: self.load_file("Epsilons Reactant"))
            self.LoadEpsilons_Product.clicked.connect(lambda: self.load_file("Epsilons Product"))
            self.LoadSpectra.clicked.connect(lambda: self.load_file("Spectral Data"))
            self.LoadLog.clicked.connect(lambda: self.load_file("Log Irr"))
            self.SaveResultsBtn.clicked.connect(self.Save_QY)
    
            ## Connect the textChanged signal to the update functions ##
            self.lineEdit_Volume.textChanged.connect(self.update_V)
            self.lineEdit_k.textChanged.connect(self.update_k_BA)
            self.lineEdit_ManPower.textChanged.connect(self.update_I0_avg)
            self.lineEdit_ManPowerError.textChanged.connect(self.update_I0_err)
            self.lineEdit_LEDWavelength.textChanged.connect(self.update_LEDw_Integration) # Integration Mode
            self.lineEdit_Threshold.textChanged.connect(self.update_threshold) # 
            self.lineEdit_WavelengthRange_Min.textChanged.connect(self.update_wl_range) # 
            self.lineEdit_WavelengthRange_Max.textChanged.connect(self.update_wl_range) # 
    
            ## Buttons for processing and plot functions ##
            self.plotEpsilonButton.clicked.connect(self.plot_epsilon)
            self.PlotSpectraButton.clicked.connect(self.plot_spectra)
            self.PlotLEDEmissionButton.clicked.connect(self.plot_LEDfull)
            self.ProcessPlotDataButton_Concentrations.clicked.connect(self.process_LED) # Concentrations ODE method
            self.CalculateFractionsButton.clicked.connect(self.Calc_Fractions)
            self.PlotFractionsResidualsButton.clicked.connect(self.plot_fractions_residuals) ## new window: plot obtained fractions with residuals
            self.ProcessPlotDataButton_Emission.clicked.connect(self.process_LED) # Emissions ODE method
            self.CalcQYButton.clicked.connect(self.Calc_QY)
            
            #### INITIALISATION ####
            self.SetButtons()
            self.SetTextfields() # Set ExpParams in text fields
            self.handle_radio_selections()
        
        ######################################################################
        
        def SetTextfields(self):
            """ Display experimental parameters in text fields """
            self.lineEdit_Volume.setText(str(ExpParams.V)) # volume
            self.lineEdit_k.setText(str(ExpParams.k_BA)) # rate constant
            self.lineEdit_ManPower.setText(str(ExpParams.I0_avg)) # Power Manual
            self.lineEdit_ManPowerError.setText(str(ExpParams.I0_err)) # Power Manual
            self.lineEdit_LEDWavelength.setText(str(ExpParams.LEDw)) # Integration Mode (default)
            self.lineEdit_AvgPower.setText(str(ExpParams.I0_avg)) # PowerProcessing: Calculated Power
            self.lineEdit_AvgPowerError.setText(str(ExpParams.I0_err)) # PowerProcessing: Error
            self.lineEdit_Threshold.setText(str(CalcSettings.threshold)) # Threshold for LED Emission spectrum
            self.lineEdit_WavelengthRange_Min.setText(str(CalcSettings.wl_low)) # Low wavelength range for Concentrations ODE METHOD
            self.lineEdit_WavelengthRange_Max.setText(str(CalcSettings.wl_high)) # High wavelength range for Concentrations ODE METHOD
    
        def SetDefaultSettings(self):
            CalcSettings.format_timestamps = Defaults.format_timestamps
            CalcSettings.CalculationMethod = Defaults.CalculationMethod
            CalcSettings.ODEMethod = Defaults.ODEMethod
            CalcSettings.PowerMethod = Defaults.PowerMethod
            CalcSettings.threshold = Defaults.threshold
            CalcSettings.BaselineCorrection_LED = Defaults.BaselineCorrection_LED
            CalcSettings.xlim_min_ProcessedData = Defaults.xlim_min_ProcessedData
            CalcSettings.xlim_max_ProcessedData = Defaults.xlim_max_ProcessedData
            CalcSettings.ylim_min_ProcessedData = Defaults.ylim_min_ProcessedData
            CalcSettings.ylim_max_ProcessedData = Defaults.ylim_max_ProcessedData
            CalcSettings.wl_low = Defaults.wl_low
            CalcSettings.wl_high = Defaults.wl_high
            ExpParams.LEDw = Defaults.LEDw
    
        def SetButtons(self):
            if CalcSettings.format_timestamps == "AHK":
                self.radioButton_Log_AHK.setChecked(True)
            elif CalcSettings.format_timestamps == "Default":
                self.radioButton_Log_Default.setChecked(True)

            if CalcSettings.ODEMethod == "Concentrations":
                self.radioButton_ODE_Concentration.setChecked(True)
            elif CalcSettings.ODEMethod == "Emission":
                self.radioButton_ODE_Emission.setChecked(True)

            if CalcSettings.PowerMethod == "Manual":
                self.radioButton_PowerManual.setChecked(True)
            elif CalcSettings.PowerMethod == "PowerProcessing":
                self.radioButton_PowerProcessing.setChecked(True)

            if CalcSettings.BaselineCorrection_LED == "ON":
                self.radioButton_blcorrLED_on.setChecked(True) # Baseline Correction of LED on
            elif CalcSettings.BaselineCorrection_LED == "OFF":
                self.radioButton_blcorrLED_off.setChecked(True) # Baseline Correction of LED off
            
            self.CalcQYButton.setEnabled(False)
            self.SaveResultsBtn.setEnabled(False)

        def SetResultTextfields(self):
            self.lineEdit_QY_RtoP.setText(f"{Results.QY_AB_opt}") # optimised QY R to P
            self.lineEdit_QY_PtoR.setText(f"{Results.QY_BA_opt}") # optimised QY P to R
            
            self.lineEdit_QYerror_RtoP.setText(f"{Results.error_QY_AB}") # error R to P
            self.lineEdit_QYerror_PtoR.setText(f"{Results.error_QY_BA}") # error P to R
            
            self.lineEdit_PSS_R.setText(f"{Results.PSS_Reactant}") # %R at PSS
            self.lineEdit_PSS_P.setText(f"{Results.PSS_Product}") # %P at PSS

        ## Update methods for the parameters
        def update_V(self):
            try:
                ExpParams.V = float(self.lineEdit_Volume.text())  # Convert the input to a float
                print(f"Updated V to {ExpParams.V}")
            except ValueError:
                pass  # Handle the case where the input is not a valid number
    
        def update_k_BA(self):
            try:
                ExpParams.k_BA = float(self.lineEdit_k.text())  # Convert the input to a float
                print(f"Updated k_BA to {ExpParams.k_BA}")
            except ValueError:
                pass
    
        def update_I0_avg(self):
            try:
                ExpParams.I0_avg = float(self.lineEdit_ManPower.text())  # Convert the input to an integer
                print(f"Updated I0_avg to {ExpParams.I0_avg}")
            except ValueError:
                pass
    
        def update_I0_err(self):
            try:
                ExpParams.I0_err = float(self.lineEdit_ManPowerError.text())  # Convert the input to an integer
                print(f"Updated I0_err to {ExpParams.I0_err}")
            except ValueError:
                pass
    
        def update_I0_avg_PP(self):
            try:
                ExpParams.I0_avg = float(self.lineEdit_AvgPower.text())  # Convert the input to an integer
                print(f"Updated I0_avg to {ExpParams.I0_avg}")
            except ValueError:
                pass
    
        def update_I0_err_PP(self):
            try:
                ExpParams.I0_err = float(self.lineEdit_AvgPowerError.text())  # Convert the input to an integer
                print(f"Updated I0_err to {ExpParams.I0_err}")
            except ValueError:
                pass

        def update_LEDw_Integration(self):
            """ For Integration Mode (default) """
            try:
                ExpParams.LEDw = int(self.lineEdit_LEDWavelength.text())  # Convert the input to an integer
                print(f"Updated LEDw to {ExpParams.LEDw}")
            except ValueError:
                pass
    
        def update_threshold(self):
            """ Update value for threshold used for LED emission spectrum """
            try:
                CalcSettings.threshold = int(self.lineEdit_Threshold.text())  # Convert the input to an integer
                print(f"Updated threshold to {CalcSettings.threshold}")
            except ValueError:
                pass
    
        def update_wl_range(self):
            """Update values for wavelength range for processing data for ODE Concentrations method"""
            try:
                CalcSettings.wl_low = int(self.lineEdit_WavelengthRange_Min.text())
                CalcSettings.wl_high = int(self.lineEdit_WavelengthRange_Max.text())
                print(f"Updated wavelength range to {CalcSettings.wl_low}-{CalcSettings.wl_high}")
            except ValueError:
                pass
    
        def handle_radio_selections(self):
            self.handle_radio_selection_Power()
            self.handle_radio_selection_ODE()
            self.handle_radio_selection_Log()
            self.handle_radio_selection_blcorrLED()
            
        def handle_radio_selection_Power(self):            
            if self.radioButton_PowerManual.isChecked(): # Power Manual Input
                self.update_I0_avg # set I0_avg to current text
                self.update_I0_err # set I0_err to current text
                
                ## Enable TextLabel and disable Load button
                self.lineEdit_ManPower.setEnabled(True)
                self.lineEdit_ManPowerError.setEnabled(True)
                
                self.loadDataButton_1.setEnabled(False)
                self.loadDataButton_2.setEnabled(False)
                self.loadDataButton_3.setEnabled(False)
                self.lineEdit_Power_1.setEnabled(False)
                self.lineEdit_PowerError_1.setEnabled(False)
                self.lineEdit_Power_2.setEnabled(False)
                self.lineEdit_PowerError_2.setEnabled(False)
                self.lineEdit_Power_3.setEnabled(False)
                self.lineEdit_PowerError_3.setEnabled(False)
                self.DeletePowerDataButton_1.setEnabled(False)
                self.DeletePowerDataButton_2.setEnabled(False)
                self.DeletePowerDataButton_3.setEnabled(False)
                self.calculatePowerButton.setEnabled(False)
                self.lineEdit_AvgPower.setEnabled(False)
                self.lineEdit_AvgPowerError.setEnabled(False)
                CalcSettings.PowerMethod = "Manual"
            
            if self.radioButton_PowerProcessing.isChecked(): # PowerProcessing Module
                self.update_I0_avg_PP # set I0_avg to current text in PowerProcessing field
                self.update_I0_err_PP # set I0_err to current text in PowerProcessing field
                
                ## Disable TextLabel and enable Load button
                self.lineEdit_ManPower.setEnabled(False) # turn off Manual Input: Power
                self.lineEdit_ManPowerError.setEnabled(False) # turn off Manual Input: Error
                
                self.loadDataButton_1.setEnabled(True)
                self.loadDataButton_2.setEnabled(True)
                self.loadDataButton_3.setEnabled(True)
                self.lineEdit_Power_1.setEnabled(True)
                self.lineEdit_PowerError_1.setEnabled(True)
                self.lineEdit_Power_2.setEnabled(True)
                self.lineEdit_PowerError_2.setEnabled(True)
                self.lineEdit_Power_3.setEnabled(True)
                self.lineEdit_PowerError_3.setEnabled(True)
                self.DeletePowerDataButton_1.setEnabled(True)
                self.DeletePowerDataButton_2.setEnabled(True)
                self.DeletePowerDataButton_3.setEnabled(True)
                self.calculatePowerButton.setEnabled(True)
                self.lineEdit_AvgPower.setEnabled(True)
                self.lineEdit_AvgPowerError.setEnabled(True)
                CalcSettings.PowerMethod = "PowerProcessing"

        def handle_radio_selection_ODE(self):
            if self.radioButton_ODE_Concentration.isChecked(): # ODE Solving Method Concentrations
                self.ProcessPlotDataButton_Concentrations.setEnabled(True)
                self.CalculateFractionsButton.setEnabled(False) # 
                self.PlotFractionsResidualsButton.setEnabled(False) # 
                self.radioButton_blcorrLED_on.setEnabled(True)
                self.radioButton_blcorrLED_off.setEnabled(True)
                # self.label_LEDwavelength.setEnabled(True)
                
                self.groupBox_Conc_Process.setEnabled(True)
                self.groupBox_blcorrLED.setEnabled(True)
                self.groupBox_Conc_wlrange.setEnabled(True)
                self.lineEdit_WavelengthRange_Min.setEnabled(True)
                self.label_Conc_wlrange_dash.setEnabled(True)
                self.lineEdit_WavelengthRange_Max.setEnabled(True)
                self.label_Conc_wlrange_nm.setEnabled(True)
                                
                self.label_Threshold.setEnabled(False) # threshold label
                self.lineEdit_Threshold.setEnabled(False) # 
                self.ProcessPlotDataButton_Emission.setEnabled(False)
                
                CalcSettings.ODEMethod = "Concentrations"
    
            if self.radioButton_ODE_Emission.isChecked(): # ODE Solving Method Emission
                self.ProcessPlotDataButton_Concentrations.setEnabled(False)
                self.CalculateFractionsButton.setEnabled(False) # SingleWavelength wavelength (nm)
                self.PlotFractionsResidualsButton.setEnabled(False) # SingleWavelength epsilon Reactant
                self.radioButton_blcorrLED_on.setEnabled(False)
                self.radioButton_blcorrLED_off.setEnabled(False)
                self.groupBox_Conc_Process.setEnabled(False)
                self.groupBox_blcorrLED.setEnabled(False)
                self.groupBox_Conc_wlrange.setEnabled(False)
                self.lineEdit_WavelengthRange_Min.setEnabled(False)
                self.label_Conc_wlrange_dash.setEnabled(False)
                self.lineEdit_WavelengthRange_Max.setEnabled(False)
                self.label_Conc_wlrange_nm.setEnabled(False)
                self.lineEdit_Threshold.setEnabled(True) # Integration wavelength (nm)
                self.ProcessPlotDataButton_Emission.setEnabled(True)
                self.label_Threshold.setEnabled(True)
                CalcSettings.ODEMethod = "Emission"
            
        def handle_radio_selection_Log(self):
            if self.radioButton_Log_Default.isChecked(): # Timestamps: Default
                CalcSettings.format_timestamps = "Default"
    
            if self.radioButton_Log_AHK.isChecked(): # Timestamps: AHK format (Crespi group)
                CalcSettings.format_timestamps = "AHK"

        def handle_radio_selection_blcorrLED(self):
            if self.radioButton_blcorrLED_on.isChecked(): # Baseline Correction of LED
                CalcSettings.BaselineCorrection_LED = "ON"
            if self.radioButton_blcorrLED_off.isChecked(): # Baseline Correction of LED
                CalcSettings.BaselineCorrection_LED = "OFF"

        ############################################################################################################

        def ClearLoadedData(self):
            """
            Clear all the loaded data and reset radio buttons.

            Returns
            -------
            None.

            """

            for module, vars_dict in self.default_state.items():
                for name, value in vars_dict.items():
                    setattr(module, name, copy.deepcopy(value))
            print("Variables reset!")

            self.SetTextfields()
            self.SetButtons()
            
        def DeletePowerData(self, count):
            """"""

            if count not in LoadedData.PowersAtCuvette:
                QtWidgets.QMessageBox.warning(self, "Error", "Data does not exist.")
                return
    
            LoadedData.PowersAtCuvette.pop(count) # remove result from dictionary
            LoadedData.ErrorsAtCuvette.pop(count) # remove result from dictionary
    
            # print(f"main-DeletePowerData_1===LoadedData.PowersAtCuvette:{LoadedData.PowersAtCuvette}")
            
            self.labels_power[count].setText("")
            self.labels_error[count].setText("")
    
        ############################################################################################################
        ############################ POWER ############################
        ############################################################################################################

        def OpenWindow_PowerProcessing(self, count):
            """Load the power data from a file and plot it in a new window."""
            LoadedData.count = count
            self.window_PP = PowerProcessing.WindowPowerProcessing(parent=self) # load Class that includes loadUi
            self.window_PP.show()
    
        def calculate_total_power(self):
            """Calculate the total power and standard deviation from all baseline-corrected data."""
            if not LoadedData.PowersAtCuvette:
                QtWidgets.QMessageBox.warning(self, "Error", "No power data has been processed yet.")
                return
    
            averages_alldatasets = list(LoadedData.PowersAtCuvette.values())
            errors_alldatasets = list(LoadedData.ErrorsAtCuvette.values())
    
            final_averaged_power = np.mean(averages_alldatasets)
            final_variance = np.sum(np.square(errors_alldatasets)) / len(averages_alldatasets)
            final_averaged_std = np.sqrt(final_variance)
            
            ExpParams.I0_avg = final_averaged_power # set to calculated power
            ExpParams.I0_err = final_averaged_std # set to calculated error
            
            self.lineEdit_AvgPower.setText(f"{ExpParams.I0_avg:.2f}") # set PowerProcessing: Power, final averaged
            self.lineEdit_AvgPowerError.setText(f"{ExpParams.I0_err:.2f}") # set PowerProcessing: Error, final averaged
            
            print(f"I0_avg: {ExpParams.I0_avg}\nI0_err: {ExpParams.I0_err}")
            
            self.Save_PowerResults()
    
        def Save_PowerResults(self):
            """ Save Power Processing results """
            
            savefile = Results.savefilename_power+".txt"
    
            allpowers = LoadedData.PowersAtCuvette
            allerrors = LoadedData.ErrorsAtCuvette
            avgdpowererror = {'Averaged Power (mW)': ExpParams.I0_avg,
                      'Averaged Error (mW)' : ExpParams.I0_err,
                      }
    
            try:
                os.remove(savefile)
            except:
                pass
            # except OSError as e:
            #     print(f"An error occurred for {e.filename} - {e.strerror}")
    
            try:
                file=savefile
                with open (file,'a') as file:
                    for i in allpowers:
                        file.write("Power "+str(i)+": "+str(allpowers[i])+'\n')
                    for i in allerrors:
                        file.write("Error "+str(i)+": "+str(allerrors[i])+'\n')
                    for i in avgdpowererror:
                        file.write(i+": "+str(avgdpowererror[i])+'\n')
            except IOError as e:
                print(f"An error occurred: {e}")
        
        ####################################################################################################################################
        ############################ LOAD DATA ############################
        ####################################################################################################################################
    
        def load_file(self, file_type):
            """Load a file (LED emission, spectral data, epsilons) based on the file type."""
            ##!!! Move to tools.load_data
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
                # Store the file path in the appropriate attribute based on the file type
                file_ext = os.path.splitext(file_name)[1].lower()
                if file_ext == '.csv':
                    self.load_csv(file_name, file_type)
                elif file_ext == '.dat':
                    self.load_dat(file_name, file_type)
                else:
                    QtWidgets.QMessageBox.warning(self, "Error", f"Unknown file type: {file_type}")
                
                ##!!! THIS MESSAGE STILL APPEARS EVEN IF A FILE IS UNSUCCESSFULLY LOADED...
                ## CHANGE SO THAT ABOVE IF/ELIF-STATEMENT NEEDS TO BE SUCCESSFUL
                QtWidgets.QMessageBox.information(self, "Success", f"{file_type} file loaded successfully!")
            except Exception as e:
                    QtWidgets.QMessageBox.critical(self, "Error", f"Failed to load {file_type} file: {e}")
    
        def load_dat(self, file_path, file_type):
            """Load the data file depending on its format (.csv or .dat)."""
            ##!!! Move to tools.load_data
            try:
                if file_type == "LED Emission":
                    LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity = LoadData.Import_LEDemission("Spectragryph", file_path)
                    LoadedData.filename_LED = file_path
                    
                elif file_type == "Epsilons Reactant":
                    LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_R_values = LoadData.Import_Epsilons("Spectragryph", file_path)
                elif file_type == "Epsilons Product":
                    LoadedData.epsilons_P_wavelengths, LoadedData.epsilons_P_values = LoadData.Import_Epsilons("Spectragryph", file_path)
                    print
                elif file_type == "Spectral Data":
                    LoadedData.SpectralData_Full, LoadedData.SpectralData_Wavelengths, LoadedData.SpectralData_Absorbance = \
                        LoadData.Import_SpectralData("Spectragryph",file_path) #HARDCODED IN THE WRONG PLACE # STILL??
                    LoadedData.filename_spectra = file_path
                    print(f"load_dat\n===LoadedData.filename_spectra: {LoadedData.filename_spectra}")
                        
                ########################################
                ##!!! ADD in case of .dat format
                # elif file_type == "Log Irr":
                #      LoadedData.timestamps = LoadData.GetTimestamps(file_path)
                else:
                        raise ValueError(f"Unknown file type: {file_type}")
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to import .dat file for {file_type}: {e}")
    
        def load_csv(self, file_path, file_type):
            """Load the data file depending on its format (.csv or .dat)."""
            try:       
                if file_type == "LED Emission": ###### A LOT OF PROBLEMS WITH '
                    LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity = LoadData.Import_LEDemission("Not", file_path)
                    LoadedData.filename_LED = file_path
                elif file_type == "Epsilons Reactant":
                    LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_R_values = LoadData.Import_Epsilons("Not", file_path)
                elif file_type == "Epsilons Product":
                    LoadedData.epsilons_P_wavelengths, LoadedData.epsilons_P_values = LoadData.Import_Epsilons("Not", file_path)
                elif file_type == "Spectral Data":
                    LoadedData.SpectralData_Full, LoadedData.SpectralData_Wavelengths,
                    LoadedData.SpectralData_Absorbance = \
                        LoadData.Import_SpectralData("Not", file_path)
                elif file_type == "Log Irr":
                    LoadedData.timestamps = LoadData.GetTimestamps(file_path)
                else:
                        raise ValueError(f"Unknown file type: {file_type}")
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to import .csv file for {file_type}: {e}")
    
        ###=========================================================================###
        ################################### PLOT ######################################
        ###=========================================================================###

        def add_new_tab(self, plot_func, title):
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
                self.tabWidget.tabCloseRequested.connect(self.tabWidget.removeTab)
    
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to add new tab: {e}")

        def plot_epsilon(self):
            """ Plot epsilons spectra (before interpolation) """
            if LoadedData.epsilons_R_wavelengths is None or LoadedData.epsilons_P_wavelengths is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load the epilons data first.")
                return
            
            def plot_func(canvas):
                """ Plot the data using MplCanvas """
                canvas.plot_EpsilonsOnly(LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_R_values,
                              LoadedData.epsilons_P_wavelengths, LoadedData.epsilons_P_values)
            
            self.add_new_tab(plot_func, "Epsilons")
    
        def plot_spectra(self):
            """ Plot spectra recorded during irradiation """
            if LoadedData.SpectralData_Full is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please first load Measurement Spectra.")
                return
    
            def plot_func(canvas):
                """ Plot the data using MplCanvas """
                canvas.PlotData_Full(LoadedData.SpectralData_Wavelengths,
                                     LoadedData.SpectralData_Absorbance)
            
            self.add_new_tab(plot_func, "Measurement Spectra")
                
        def plot_LEDfull(self):
            """ Plot LED emission spectrum (full) """
            if LoadedData.LEDemission_wavelengths is None or LoadedData.LEDemission_intensity is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load LED emission file.")
                return
            
            if ExpParams.LEDw == 0:
                QtWidgets.QMessageBox.warning(self, "Error", "Please set nominal wavelength (nm).")
                return
    
            def plot_func(canvas):
                """ Plot the data using MplCanvas """
                canvas.plot_LEDemission_full(LoadedData.LEDemission_wavelengths,LoadedData.LEDemission_intensity)
            
            self.add_new_tab(plot_func, "LED emission")

        def plot_fractions(self):
            """ Plot fractions obtained from fitting epsilons spectra to all spectra during irradiation"""
            if not Results.fractions_R or not Results.fractions_P:
                QtWidgets.QMessageBox.warning(self, "Error", "Something went wrong with the calculation of the fractions.")
                return
    
            def plot_func(canvas):
                """ Plot the data using MplCanvas """
                canvas.PlotFractions(LoadedData.SpectralDataCut_Wavelengths,
                                     LoadedData.SpectralDataCut_Abs,
                                     Results.reconstructed_spectra_fractions,
                                     Results.fractions_R,
                                     Results.fractions_P)
            
            self.add_new_tab(plot_func, "Fractions")


        def plot_fractions_residuals(self):
        ##!!! ADD function that starts a pop-up window 
            ## to show the residuals for each spectrum    
            '''
            Pop-up window that shows residuals for each fitted spectra
                        

            '''

        ################################################################################
        ################################################################################
        ################################################################################

        def process_LED(self):
            """
            Process and visualize the data based on the loaded files.
            - if ODEMethod = Emission: 
                Process spectra using the LED emission spectra for the wavelength limits.
            - if ODEMethod = Concentrations: 
                use the limits of the epsilons spectra
                
            Then:
            Process Spectral data: cut to part of spectrum according to limits set beforehand
                
            """
            self.CalculateFractionsButton.setEnabled(False)
            self.CalcQYButton.setEnabled(False)
            self.SaveResultsBtn.setEnabled(False)

            if LoadedData.LEDemission_wavelengths is None or LoadedData.LEDemission_intensity is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load LED emission file.")
                return

            if ExpParams.LEDw == 0:
                QtWidgets.QMessageBox.warning(self, "Error", "Please set nominal wavelength (nm).")
                return

            if LoadedData.SpectralData_Full is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Measurement Spectra.")
                return

            if LoadedData.epsilons_R_wavelengths is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Epsilons Reactant.")
                return

            if LoadedData.epsilons_P_wavelengths is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Epsilons Product.")
                return

            ''' Obtain smoothed and non-negatived LED emission data '''
            LoadedData.LEDemission_intensity_proc = \
                Integration.Process_LEDemission(LoadedData.LEDemission_wavelengths,
                                                LoadedData.LEDemission_intensity)
            
            ########################################
            if CalcSettings.ODEMethod == "Emission":
                ''' Find indices for wavelengths low and high end of LED emission data'''
                (Integration.wavelength_low, Integration.wavelength_high) = \
                    Integration.LEDemission_WavelengthLimits(
                        LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity_proc, CalcSettings.threshold)
            
            elif CalcSettings.ODEMethod == "Concentrations":
                ''' Find indices for wavelengths low and high end of epsilons data'''
                (Integration.wavelength_low, Integration.wavelength_high) = \
                    Integration.Epsilons_WavelengthLimits(
                        LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_P_wavelengths)
            
            else:
                QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the CalcSettings.ODEMethod variable")
            ########################################
            ''' 
            Process Spectral data: cut to part of spectrum according to 
            - LED emission band OR
            - epsilons data range 
            '''
            (LoadedData.SpectralDataCut_Wavelengths, 
             LoadedData.SpectralDataCut_Abs,
             LoadedData.SpectralDataCut_Index) = \
                Integration.Process_SpectralData(LoadedData.SpectralData_Full,
                                                 Integration.wavelength_low,
                                                 Integration.wavelength_high,
                                                 ExpParams.LEDw)

            (LoadedData.epsilons_R_interp, 
             LoadedData.epsilons_P_interp, 
             LoadedData.emission_interp) = \
                Integration.Interpolate_Epsilons(LoadedData.SpectralDataCut_Wavelengths,
                                                 LoadedData.epsilons_R_wavelengths,
                                                 LoadedData.epsilons_R_values,
                                                 LoadedData.epsilons_P_wavelengths,
                                                 LoadedData.epsilons_P_values,
                                                 LoadedData.LEDemission_wavelengths,
                                                 LoadedData.LEDemission_intensity_proc)
            ########################################
            def plot_func(canvas):
                """ 
                Plot the data using MplCanvas.
                The function PlotData_Cut plots the full spectra in grey,
                    and the cut spectra in colour
                """
                ##!!! probably can remove: LoadedData.LEDemission_wavelengths, 
                 
                canvas.PlotData_Cut(LoadedData.SpectralData_Absorbance, LoadedData.SpectralData_Wavelengths,
                                    LoadedData.SpectralDataCut_Abs, LoadedData.SpectralDataCut_Wavelengths,
                                    LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity,
                                    LoadedData.emission_interp) ## use interpolated (and cut) LED emission data
            
            if CalcSettings.ODEMethod == "Concentrations":
                self.CalculateFractionsButton.setEnabled(True)
            
            self.add_new_tab(plot_func, "Processed Spectra")
        #######################################


        def Calc_Fractions(self):
            '''
            Calculate Fractions
            Plot results
            Save results
            '''       
            self.CalcQYButton.setEnabled(False)
            self.SaveResultsBtn.setEnabled(False)
            
            if LoadedData.SpectralDataCut_Abs is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please perform processing of spectra.")
                return

            if LoadedData.epsilons_R_interp is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please perform processing of spectra.")
                return

            if LoadedData.epsilons_P_interp is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please perform processing of spectra.")
                return
            
            Results.fractions_R, Results.fractions_P, Results.reconstructed_spectra_fractions = \
                Fractions.CalculateFractions(LoadedData.SpectralDataCut_Abs,
                                             LoadedData.SpectralDataCut_Wavelengths,
                                             LoadedData.epsilons_R_interp,
                                             LoadedData.epsilons_P_interp)
            
            self.plot_fractions() ## plot retrieved fractions
            self.CalcQYButton.setEnabled(True)
            
            ### Save as .csv ###
            Fractions.Save_FractionsResults(Results.fractions_R, Results.fractions_P,
                                            LoadedData.filename_spectra)
            
        #######################################    

        def Calc_QY(self, canvas):
            """
            Calculate quantum yields by numerically solving the differential equations.
            Then calculate the concentrations, and plot the results.
            """
            
            self.SaveResultsBtn.setEnabled(False)

            if ExpParams.V == 0.0:
                QtWidgets.QMessageBox.warning(self, "Error", "Please set Volume.")
                return
            
            if ExpParams.I0_avg == 0.0:
                QtWidgets.QMessageBox.warning(self, "Error", "Please check Power.")
                return
            
            if LoadedData.LEDemission_wavelengths is None or LoadedData.LEDemission_intensity is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load LED emission file.")
                return
            
            if LoadedData.emission_interp is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please process LED emission spectrum.")
                return
    
            if LoadedData.timestamps is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load the Timestamps log file.")
                return

            if not Results.fractions_R or not Results.fractions_P:
                QtWidgets.QMessageBox.warning(self, "Error", "Please calculate the fractions first.")
                return

            ##!!! ADD CHECK: len(measurements) should equal len(timestamps)
            
            Datasets.I0_list = Datasets.ListPowers(ExpParams.I0_avg, ExpParams.I0_err) ## Make list of powers
            ######################################################################

            if CalcSettings.ODEMethod == "Emission":
                ## Create parameters needed for fitting
                (Datasets.initial_conc_R, Datasets.initial_conc_P, 
                 Datasets.wavelengths_meters, Datasets.normalized_emission) = \
                    Integration.CreateParameters(LoadedData.SpectralDataCut_Abs, 
                                                 LoadedData.SpectralDataCut_Wavelengths,
                                                LoadedData.epsilons_R_interp,
                                                LoadedData.emission_interp)
        
                
                Datasets.N, Datasets.fit_results = Integration.MinimizeQYs(Datasets.I0_list, 
                                                        Datasets.normalized_emission,
                                                        Datasets.wavelengths_meters, 
                                                        Datasets.initial_conc_R, Datasets.initial_conc_P,
                                                        LoadedData.timestamps, LoadedData.SpectralDataCut_Abs,
                                                        LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                        ExpParams.V)
                self.CalcQYButton.setEnabled(True)
                self.SaveResultsBtn.setEnabled(True)
                self.Extract_QY() # extract QY results and display
            ######################################################################    
            elif CalcSettings.ODEMethod == "Concentrations":
                print(f"Calc_QY ODEMethod:\n{CalcSettings.ODEMethod}")

                ## Create parameters needed for fitting
                (Datasets.initial_conc_R, Datasets.initial_conc_P, Datasets.concs_RP,
                 Datasets.wavelengths_meters, Datasets.normalized_emission) = \
                    Integration.CreateParameters_Conc(LoadedData.SpectralDataCut_Abs, 
                                                 LoadedData.SpectralDataCut_Wavelengths,
                                                 Results.fractions_R, Results.fractions_P,
                                                LoadedData.epsilons_R_interp,LoadedData.epsilons_P_interp,
                                                LoadedData.emission_interp)
                
                Datasets.N, Datasets.fit_results = Integration.MinimizeQYs_Conc(Datasets.I0_list, 
                                                        Datasets.normalized_emission,
                                                        Datasets.wavelengths_meters, 
                                                        Datasets.initial_conc_R, Datasets.initial_conc_P,
                                                        LoadedData.timestamps, Datasets.concs_RP,
                                                        LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                        ExpParams.V)
                self.CalcQYButton.setEnabled(True)
                self.SaveResultsBtn.setEnabled(True)
                self.Extract_QY() # extract QY results and display
            else:
                QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the CalcSettings.ODEMethod variable")
            ######################################################################

        def Extract_QY(self):
            '''
            Extract results from minimisation.

            Returns
            -------
            None.

            '''

            ## Extract results from the fit
            (Results.QY_AB_opt, Results.QY_BA_opt, Results.error_QY_AB, 
             Results.error_QY_BA) = ExtractResults.ExtractResults(Datasets.fit_results)

    
            ## Calculate optimized concentrations
            (Results.conc_opt, Results.PSS_Reactant, 
             Results.PSS_Product) = ExtractResults.CalculateConcentrations(Datasets.wavelengths_meters,
                                                Datasets.initial_conc_R, Datasets.initial_conc_P, 
                                                LoadedData.timestamps,
                                                Results.QY_AB_opt, Results.QY_BA_opt, 
                                                LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                Datasets.N, ExpParams.V)
    
                                                               
            ## Calculate total absorbance and residuals
            if CalcSettings.ODEMethod == "Emission":
                Results.total_abs_fit, Results.residuals = \
                    ExtractResults.GetFittedAbs(Datasets.fit_results,
                                                Results.conc_opt,
                                                LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                LoadedData.timestamps,
                                                LoadedData.SpectralDataCut_Wavelengths)
            else:
                pass
            ## Not needed for ODEMethod="Concentrations" => is done in the plotting.py function
            
            ######################################################################
            self.SetResultTextfields() ## Update textfields to show results
            self.add_new_tab(self.Plot_QY, "QY") ## Plot and save the results
            QtWidgets.QMessageBox.information(self, "Success", "Results extracted and plotted!")
    
        def Plot_QY(self, canvas):
            if CalcSettings.ODEMethod == "Emission":
                canvas.PlotResults(ExpParams.LEDw,
                                   LoadedData.timestamps,
                                   Results.conc_opt,
                                   LoadedData.SpectralDataCut_Abs,
                                   LoadedData.SpectralDataCut_Index,
                                   Results.total_abs_fit,
                                   Results.residuals,
                                   Results.QY_AB_opt, Results.QY_BA_opt,
                                   Results.error_QY_AB, Results.error_QY_BA,
                                   CalcSettings.ODEMethod)

            elif CalcSettings.ODEMethod == "Concentrations":
                canvas.PlotResults_Conc(ExpParams.LEDw,
                                   LoadedData.timestamps,
                                   Results.conc_opt,
                                   Datasets.concs_RP,
                                   LoadedData.SpectralDataCut_Wavelengths,
                                   LoadedData.epsilons_R_interp,
                                   LoadedData.epsilons_P_interp,
                                   Results.QY_AB_opt, Results.QY_BA_opt,
                                   Results.error_QY_AB, Results.error_QY_BA,
                                   CalcSettings.ODEMethod)
            
            else:
                QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the CalcSettings.ODEMethod variable")
        
        def Save_QY(self):
            """ Save results: plots """
            if Results.QY_AB_opt is None or Results.conc_opt is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please perform Calculate QY first.")
                return

            options = QtWidgets.QFileDialog.Options()
    
            ## File dialog for selecting files
            savefilename, _ = QtWidgets.QFileDialog.getSaveFileName(self, 
                                                                 "Choose name of savefile", "",
                                                                 options=options)
            
            path_split=savefilename.split('/')[0:-1] # leave only filepath (remove name)
            path='\\'.join(path_split) # re-join into string
    
            end_nameonly=savefilename.split('/')[-1] # only filename
            end=f"Results_{end_nameonly}_{CalcSettings.CalculationMethod}" # name with added info
    
            Results.savefilename = f"{path}\\{end}"
            
            canvas = MplCanvas(self)
            if not savefilename:
                return
            else:
                if CalcSettings.ODEMethod == "Emission":
                    canvas.PlotResults(ExpParams.LEDw,
                                       LoadedData.timestamps,
                                       Results.conc_opt,
                                       LoadedData.SpectralDataCut_Abs,
                                       LoadedData.SpectralDataCut_Index,
                                       Results.total_abs_fit,
                                       Results.residuals,
                                       Results.QY_AB_opt, Results.QY_BA_opt,
                                       Results.error_QY_AB, Results.error_QY_BA,
                                       CalcSettings.ODEMethod,
                                       SaveResults = "Yes",
                                       SaveFileName = Results.savefilename)
                
                elif CalcSettings.ODEMethod == "Concentrations":
                    canvas.PlotResults_Conc(ExpParams.LEDw,
                                       LoadedData.timestamps,
                                       Results.conc_opt,
                                       Datasets.concs_RP,
                                       LoadedData.SpectralDataCut_Wavelengths,
                                       LoadedData.epsilons_R_interp,
                                       LoadedData.epsilons_P_interp,
                                       Results.QY_AB_opt, Results.QY_BA_opt,
                                       Results.error_QY_AB, Results.error_QY_BA,
                                       CalcSettings.ODEMethod,
                                       SaveResults = "Yes",
                                       SaveFileName = Results.savefilename)
        
                self.Save_Results()
        
                QtWidgets.QMessageBox.information(self, "Success", f"{Results.savefilename} file saved successfully!")
    
        def Save_Results(self):
            """ Save results and all the parameters and settings used for the calculation """
            savefile = Results.savefilename+".txt"
    
            dict_results = {'PSS_Reactant (%)': Results.PSS_Reactant,
                     'PSS_Product (%)' : Results.PSS_Product,
                     'QY_AB_opt (%)' : Results.QY_AB_opt,
                     'QY_BA_opt (%)' : Results.QY_BA_opt,
                     'error_QY_AB (%)' : Results.error_QY_AB,
                     'error_QY_BA (%)' : Results.error_QY_BA}
    
            dict_expparams = {'Volume (ml)': ExpParams.V,
                              'k thermal back-reaction (s-1)': ExpParams.k_BA,
                              'Power average (mW)': ExpParams.I0_avg,
                              'Power error (mW)': ExpParams.I0_err,
                              'Wavelength of irradiation': ExpParams.LEDw}
            
            ##!!! ADD DATA DEPENDING ON METHOD
                ### if Emission: add Threshold
                ### if Concentrations: add ...
    
            dict_calcsettings = {'Calculation Method': CalcSettings.ODEMethod,
                                 'Threshold': CalcSettings.threshold}
    
            try:
                os.remove(savefile)
            except:
                pass
    
            try:
                file=savefile
                
                with open (file,'a') as file:
                    for i in dict_results:
                        file.write(i+": "+str(dict_results[i])+'\n')
                    # file.write('==== Experimental Parameters ====\n')
                    for i in dict_expparams:
                        file.write(i+": "+str(dict_expparams[i])+'\n')
                    for i in dict_calcsettings:
                        file.write(i+": "+str(dict_calcsettings[i])+'\n')
            except IOError as e:
                print(f"An error occurred: {e}")
    #####################################################################
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
