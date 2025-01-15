"""
@author: Alfredo and Jorn
=============================================================================
autoQuant
=============================================================================
TO DO:
[] Add feature to input starting concentrations
=============================================================================

"""
import sys, os
import numpy as np
import pandas as pd
# from scipy.optimize import curve_fit #change in baseline correction

from PyQt5 import QtWidgets, uic
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

import QY.Integration as Integration
import QY.ExpParams as ExpParams
import QY.Results as Results
import QY.LoadedData as LoadedData
import QY.SingleWavelength as SingleWavelength

import tools.load_data as LoadData
from tools.plotting import MplCanvas
import tools.extractresults as ExtractResults
#from tools.style import apply_dark_theme

import PowerProcessing.WindowPowerProcessing as WindowPowerProcessing

class autoQuant(QtWidgets.QMainWindow):
        ############################
        #Error Handling  ***********
        ############################
    def __init__(self):
        super(autoQuant, self).__init__()
        uic.loadUi('UIs/autoQuant.ui', self)  # Load the UI file you provided

        ############################
        #Change style  *************
        ############################
        #self.setStyleSheet(get_stylesheet())

        ##!!! REMOVE UNNECESSARY VARIABLES HERE

        self.labels_power = {1: self.plainTextEdit_Power_1,
                                   2: self.plainTextEdit_Power_2,
                                   3: self.plainTextEdit_Power_3}
        self.labels_error = {1: self.plainTextEdit_PowerError_1,
                                   2: self.plainTextEdit_PowerError_2,
                                   3: self.plainTextEdit_PowerError_3}
        
        ## Default Calculation Method is integration: see ExpParams

        ## Button connections
        self.loadDataButton_1.clicked.connect(lambda: self.OpenWindow_PowerProcessing(1))
        self.loadDataButton_2.clicked.connect(lambda: self.OpenWindow_PowerProcessing(2))
        self.loadDataButton_3.clicked.connect(lambda: self.OpenWindow_PowerProcessing(3))
        
        self.DeletePowerDataButton_1.clicked.connect(lambda: self.DeletePowerData(1))
        self.DeletePowerDataButton_2.clicked.connect(lambda: self.DeletePowerData(2))
        self.DeletePowerDataButton_3.clicked.connect(lambda: self.DeletePowerData(3))

        self.calculatePowerButton.clicked.connect(self.calculate_total_power) ## Calculates average power+error

        #######

        self.radioButton.toggled.connect(self.handle_radio_selection)
        self.radioButton_2.toggled.connect(self.handle_radio_selection)
        self.radioButton_3.toggled.connect(self.handle_radio_selection)
        self.radioButton_4.toggled.connect(self.handle_radio_selection)

        # Connecting buttons to their respective methods
        ##!!! implement the SingleWavelength option as well
        self.LoadLED.clicked.connect(lambda: self.load_file("LED Emission"))
        self.LoadEpsilons_Reactant.clicked.connect(lambda: self.load_file("Epsilons Reactant"))
        self.LoadEpsilons_Product.clicked.connect(lambda: self.load_file("Epsilons Product"))
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

        ## Connect the textChanged signal to the update functions ##
        self.plainTextEdit.textChanged.connect(self.update_V)
        self.plainTextEdit_2.textChanged.connect(self.update_k_BA)
        self.plainTextEdit_3.textChanged.connect(self.update_I0_avg)
        self.plainTextEdit_4.textChanged.connect(self.update_I0_err)
        self.plainTextEdit_5.textChanged.connect(self.update_LEDw_Integration) # Integration Mode
        self.plainTextEdit_7.textChanged.connect(self.update_LEDw_SingleWl) # SingleWavelength Mode
        self.plainTextEdit_epsilonR.textChanged.connect(self.update_epsilon_R) # Epsilon Reactant (SingleWavelength)
        self.plainTextEdit_epsilonP.textChanged.connect(self.update_epsilon_P) # Epsilon Product (SingleWavelength)
        self.plainTextEdit_9.textChanged.connect(self.update_threshold) # 

        ## Adding new buttons for custom plot functions from your scripts ##
        self.ProcessPlotDataButton.clicked.connect(self.process_LED)
        self.plotEpsilonButton.clicked.connect(self.plot_epsilon)
        self.PlotSpectraButton.clicked.connect(self.plot_spectra)
        self.PlotLEDEmissionButton.clicked.connect(self.plot_LEDfull)

        self.CalcQYButton.clicked.connect(self.Calc_QY)

        self.radioButton_3.setEnabled(True) # Power Manual: default enabled
        self.loadDataButton_1.setEnabled(False) # default Manual Power input
        self.loadDataButton_2.setEnabled(False) # default Manual Power input
        self.loadDataButton_3.setEnabled(False) # default Manual Power input
        self.DeletePowerDataButton_1.setEnabled(False)
        self.DeletePowerDataButton_2.setEnabled(False)
        self.DeletePowerDataButton_3.setEnabled(False)
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
            
            self.loadDataButton_1.setEnabled(False)
            self.loadDataButton_2.setEnabled(False)
            self.loadDataButton_3.setEnabled(False)
            self.plainTextEdit_Power_1.setEnabled(False)
            self.plainTextEdit_PowerError_1.setEnabled(False)
            self.plainTextEdit_Power_2.setEnabled(False)
            self.plainTextEdit_PowerError_2.setEnabled(False)
            self.plainTextEdit_Power_3.setEnabled(False)
            self.plainTextEdit_PowerError_3.setEnabled(False)
            self.DeletePowerDataButton_1.setEnabled(False)
            self.DeletePowerDataButton_2.setEnabled(False)
            self.DeletePowerDataButton_3.setEnabled(False)
            self.calculatePowerButton.setEnabled(False)
            self.plainTextEdit_6.setEnabled(False)
            self.plainTextEdit_8.setEnabled(False)
        
        if self.radioButton_4.isChecked(): # PowerProcessing
            self.update_I0_avg_PP # set I0_avg to current text in PowerProcessing field
            self.update_I0_err_PP # set I0_err to current text in PowerProcessing field
            
            ## Disable TextLabel and enable Load button
            self.plainTextEdit_3.setEnabled(False) # turn off Manual Input: Power
            self.plainTextEdit_4.setEnabled(False) # turn off Manual Input: Error
            
            self.loadDataButton_1.setEnabled(True)
            self.loadDataButton_2.setEnabled(True)
            self.loadDataButton_3.setEnabled(True)
            self.plainTextEdit_Power_1.setEnabled(True)
            self.plainTextEdit_PowerError_1.setEnabled(True)
            self.plainTextEdit_Power_2.setEnabled(True)
            self.plainTextEdit_PowerError_2.setEnabled(True)
            self.plainTextEdit_Power_3.setEnabled(True)
            self.plainTextEdit_PowerError_3.setEnabled(True)
            self.DeletePowerDataButton_1.setEnabled(True)
            self.DeletePowerDataButton_2.setEnabled(True)
            self.DeletePowerDataButton_3.setEnabled(True)
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
            ExpParams.CalculationMethod = "SingleWavelength"

        if self.radioButton_2.isChecked(): # LED Integration Mode
            self.update_LEDw_Integration # set variable to current text
            self.plainTextEdit_7.setEnabled(False) # SingleWavelength wavelength (nm)
            self.plainTextEdit_epsilonR.setEnabled(False) # SingleWavelength epsilon Reactant
            self.plainTextEdit_epsilonP.setEnabled(False) # SingleWavelength epsilon Product
            self.plainTextEdit_5.setEnabled(True) # Integration wavelength (nm)
            self.LoadLED.setEnabled(True)
            ExpParams.CalculationMethod = "Integration"

    def OpenWindow_PowerProcessing(self, count):
        """Load the power data from a file and plot it in a new window."""
        LoadedData.count = count
        self.window_PP = WindowPowerProcessing.WindowPowerProcessing(parent=self) # load Class that includes loadUi
        self.window_PP.show()

    def DeletePowerData(self, count):
        """"""
        print(f"main-DeletePowerData===count:{count} and type:{type(count)}")
        print(f"main-DeletePowerData_1===LoadedData.PowersAtCuvette:{LoadedData.PowersAtCuvette}")

        if count not in LoadedData.PowersAtCuvette:
            QtWidgets.QMessageBox.warning(self, "Error", "Data does not exist.")
            return

        LoadedData.PowersAtCuvette.pop(count) # remove result from dictionary
        LoadedData.ErrorsAtCuvette.pop(count) # remove result from dictionary

        print(f"main-DeletePowerData_1===LoadedData.PowersAtCuvette:{LoadedData.PowersAtCuvette}")
        
        self.labels_power[count].setPlainText("")
        self.labels_error[count].setPlainText("")

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

        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", f"Failed to add new tab: {e}")

    ############################################################################################################
    ############################ plot_sections ############################
    ############################################################################################################

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
        
        self.plainTextEdit_6.setPlainText(f"{ExpParams.I0_avg:.2f}") # set PowerProcessing: Power, final averaged
        self.plainTextEdit_8.setPlainText(f"{ExpParams.I0_err:.2f}") # set PowerProcessing: Error, final averaged
        
        print(f"I0_avg: {ExpParams.I0_avg}\nI0_err: {ExpParams.I0_err}")
        
        self.Save_PowerResults()

    def Save_PowerResults(self):
        """ Save results: """
        
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
            
            ##!!! THIS MESSAGE STILL APPEARS EVEN IF A FILE IS UNSUCCESSFULLY LOADED...
            QtWidgets.QMessageBox.information(self, "Success", f"{file_type} file loaded successfully!")
        except Exception as e:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to load {file_type} file: {e}")

    def load_dat(self, file_path, file_type):
        """Load the data file depending on its format (.csv or .dat)."""
        try:
            if file_type == "LED Emission":
                LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity = LoadData.Import_LEDemission("Spectragryph", file_path)
                LoadedData.filename_LED = file_path
                
            elif file_type == "Epsilons Reactant":
                LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_R_values = LoadData.Import_Epsilons("Spectragryph", file_path)
            elif file_type == "Epsilons Product":
                LoadedData.epsilons_P_wavelengths, LoadedData.epsilons_P_values = LoadData.Import_Epsilons("Spectragryph", file_path)
            elif file_type == "Spectral Data":
                LoadedData.SpectralData_Full, LoadedData.SpectralData_Wavelengths, LoadedData.SpectralData_Absorbance = \
                    LoadData.Import_SpectralData("Spectragryph",file_path) #HARDCODED IN THE WRONG PLACE # STILL??
                    
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
###=========================================================================###

    ##!!! DEFINE OUTSIDE OF CLASS    
    def plot_epsilon(self):
        """ Plot epsilons spectra (before interpolation) """
        if LoadedData.epsilons_R_wavelengths is None or LoadedData.epsilons_P_wavelengths is None:
            QtWidgets.QMessageBox.warning(self, "Error", "Please load the epilons data first.")
            return
        
        def plot_func(canvas):
            """ Plot the data using MplCanvas """
            canvas.plot_EpsilonsOnly(LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_R_values,
                          LoadedData.epsilons_P_wavelengths, LoadedData.epsilons_P_values)
        
        self.add_new_tab(plot_func, "Epsilons (before interpolation)")

    ##!!! DEFINE OUTSIDE OF CLASS    
    def plot_spectra(self):
        """ Plot spectra recorded during irradiation """
        if LoadedData.SpectralData_Full is None:
            QtWidgets.QMessageBox.warning(self, "Error", "Please first load Spectra during Irradiation.")
            return

        def plot_func(canvas):
            """ Plot the data using MplCanvas """
            canvas.plot_DataFull(LoadedData.SpectralData_Wavelengths,
                                 LoadedData.SpectralData_Absorbance)
        
        self.add_new_tab(plot_func, "Spectra during Irradiation")
            
    ##!!! DEFINE OUTSIDE OF CLASS    
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
        
        self.add_new_tab(plot_func, "LED emission (full)")


    def process_LED(self):
        """Process and visualize the data based on the loaded files."""
        
        ## Integration mode
        if ExpParams.CalculationMethod == "Integration":
            if LoadedData.LEDemission_wavelengths is None or LoadedData.LEDemission_intensity is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load LED emission file.")
                return

            if LoadedData.SpectralData_Full is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Spectra during Irradiation.")
                return

            if LoadedData.epsilons_R_wavelengths is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Epsilons Reactant.")
                return

            if LoadedData.epsilons_P_wavelengths is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Epsilons Product.")
                return
            
            if ExpParams.LEDw == 0:
                QtWidgets.QMessageBox.warning(self, "Error", "Please set nominal wavelength (nm).")
                return
            
            ########################################
            threshold_LED = ExpParams.threshold
            LoadedData.LEDemission_intensity_proc, Integration.LEDindex_first, Integration.LEDindex_last, Integration.wavelength_low, Integration.wavelength_high = \
                Integration.Processing_LEDemission(
                    LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity, threshold_LED)
            
            ########################################
            ## Process Spectral data: cut to part of spectrum according to LED emission band
            LoadedData.SpectralDataCut_Wavelengths, LoadedData.SpectralDataCut_Abs, LoadedData.SpectralDataCut_Index = \
                Integration.Process_SpectralData(LoadedData.SpectralData_Full, Integration.wavelength_low, Integration.wavelength_high, ExpParams.LEDw)
        
        ##!!! ADJUST CODE TO PLOT SPECTRA AND JUST INDICATE PART OF SPECTRUM (with a box)
            self.add_new_tab(self.PlotData_Cut, "LED Emission and Spectral Data")
        
        
            ##!!! MOVED HERE
            LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp, LoadedData.emission_interp = Integration.Interpolate_Epsilons(LoadedData.SpectralDataCut_Wavelengths,
                         LoadedData.epsilons_R_wavelengths, LoadedData.epsilons_R_values,
                         LoadedData.epsilons_P_wavelengths, LoadedData.epsilons_P_values,
                         LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity_proc)
        
        ## SingleWavelength mode
        elif ExpParams.CalculationMethod == "SingleWavelength":
            if LoadedData.SpectralData_Full is None:
                QtWidgets.QMessageBox.warning(self, "Error", "Please load Spectra.")
                return


        ##!!! WORKING ON IT HERE NOT FINISHED YET AAAHHH
            LoadedData.SpectralDataCut_Wavelengths, LoadedData.SpectralDataCut_Abs, LoadedData.SpectralDataCut_Index = \
                SingleWavelength.Process_SpectralData(LoadedData.SpectralData_Full, Integration.wavelength_low, Integration.wavelength_high, ExpParams.LEDw)

            ##!!! ADJUST CODE TO PLOT SPECTRA AND ONLY INDICATE EITHER
                ## VERTICAL LINE: SINGLE WAVELENGTH
            # self.add_new_tab(self.PlotData_Cut, "LED Emission and Spectral Data")
        
        else:
            QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the self.CalculationMethod variable")
        
        
        
        #############
        
        ##!!! DEFINE OUTSIDE OF CLASS

    def PlotData_Cut(self,canvas):
            
        canvas.PlotData_Cut(LoadedData.SpectralDataCut_Abs, LoadedData.SpectralDataCut_Wavelengths,
            LoadedData.LEDemission_wavelengths, LoadedData.LEDemission_intensity,
            Integration.LEDindex_first, Integration.LEDindex_last, 
            LoadedData.LEDemission_intensity_proc)
        

    def Calc_QY(self, canvas):
        """
        Calculate quantum yields by numerically solving the differential equations.
        Then calculate the concentrations, and plot the results.
        """
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
        
        ## Make list of powers
        I0_avg = ExpParams.I0_avg
        I0_err = ExpParams.I0_err
        I0_list = [I0_avg, I0_avg+I0_err, I0_avg-I0_err]

        if ExpParams.CalculationMethod == "Integration":
            ## Create parameters needed for fitting
            initial_conc_A, initial_conc_B, spectraldata_meters, normalized_emission = \
                Integration.CreateParameters(LoadedData.SpectralDataCut_Abs, LoadedData.SpectralDataCut_Wavelengths,
                                            LoadedData.epsilons_R_interp, LoadedData.emission_interp)
    
            
            N, fit_results = Integration.MinimizeQYs(I0_list, normalized_emission,
                                                    spectraldata_meters, 
                                                    initial_conc_A, initial_conc_B,
                                                    LoadedData.timestamps, LoadedData.SpectralDataCut_Abs,
                                                    LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                    ExpParams.V)

        elif ExpParams.CalculationMethod == "SingleWavelength":
            ##!!! WORKING ON THIS

            print(f"Calc_QY SingleWavelength SpectralData_Abs:\n{LoadedData.SpectralDataCut_Abs}")
            
            ## Create parameters needed for fitting
            initial_conc_A, initial_conc_B, LoadedData.SpectralData_AbsAtLEDw = \
                SingleWavelength.CreateParameters(LoadedData.SpectralData_Absorbance, ExpParams.LEDw, ExpParams.epsilon_R)
            
            N, fit_results = SingleWavelength.MinimizeQYs(I0_list, ExpParams.LEDw,
                                                          initial_conc_A, initial_conc_B,
                                                          LoadedData.timestamps, LoadedData.SpectralData_AbsAtLEDw,
                                                          ExpParams.epsilon_R, ExpParams.epsilon_P,
                                                          ExpParams.V)
            
        else:
            QtWidgets.QMessageBox.warning(self, "Error", "Something wrong with the self.CalculationMethod variable")
            


        ######################################################################

        ## Extract results from the fit
        Results.QY_AB_opt, Results.QY_BA_opt, Results.error_QY_AB, Results.error_QY_BA = ExtractResults.ExtractResults(fit_results)


        ##!!! MOVE TO SEPARATE PY SCRIPT (TOOLS OR SOMETHING)
        ## Calculate optimized concentrations
        Results.conc_opt, Results.PSS_Reactant, Results.PSS_Product = Integration.CalculateConcentrations(spectraldata_meters,
                                                    initial_conc_A, initial_conc_B, 
                                                    LoadedData.timestamps,
                                                    Results.QY_AB_opt, Results.QY_BA_opt, 
                                                    LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                    N, ExpParams.V)

        ## Calculate total absorbance and residuals
        Results.total_abs_fit, Results.residuals = Integration.GetFittedAbs(fit_results, Results.conc_opt,
                                                            LoadedData.epsilons_R_interp, LoadedData.epsilons_P_interp,
                                                            LoadedData.timestamps,
                                                            LoadedData.SpectralDataCut_Wavelengths)
        
        ######################################################################
        
        ## Update labels
        ##!!! MAKE THIS A SEPARATE FUNCTION
        self.textEdit_QY_RtoP.setText(f"{Results.QY_AB_opt:.3f}") # optimised QY R to P
        self.textEdit_QY_PtoR.setText(f"{Results.QY_BA_opt:.3f}") # optimised QY P to R
        
        self.textEdit_QYerror_RtoP.setText(f"{Results.error_QY_AB:.3f}") # error R to P
        self.textEdit_QYerror_PtoR.setText(f"{Results.error_QY_BA:.3f}") # error P to R
        
        self.textEdit_PSS_R.setText(f"{Results.PSS_Reactant:.1f}") # %R at PSS
        self.textEdit_PSS_P.setText(f"{Results.PSS_Product:.1f}") # %P at PSS
        
        ## Plot and save the results
        self.add_new_tab(self.Plot_QY, "QY")
        QtWidgets.QMessageBox.information(self, "Success", "Results extracted and plotted!")

    def Plot_QY(self, canvas):
        canvas.PlotResults(ExpParams.LEDw,
                           LoadedData.timestamps,
                           Results.conc_opt,
                           LoadedData.SpectralDataCut_Abs,
                           LoadedData.SpectralDataCut_Index,
                           Results.total_abs_fit,
                           Results.residuals,
                           Results.QY_AB_opt, Results.QY_BA_opt,
                           Results.error_QY_AB, Results.error_QY_BA,
                           ExpParams.CalculationMethod)

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
        end=f"Results_{end_nameonly}_{ExpParams.CalculationMethod}" # name with added info

        Results.savefilename = f"{path}\\{end}"
        
        canvas = MplCanvas(self)
        if Results.savefilename is None or '':
            return
        canvas.PlotResults(ExpParams.LEDw,
                           LoadedData.timestamps,
                           Results.conc_opt,
                           LoadedData.SpectralDataCut_Abs,
                           LoadedData.SpectralDataCut_Index,
                           Results.total_abs_fit,
                           Results.residuals,
                           Results.QY_AB_opt, Results.QY_BA_opt,
                           Results.error_QY_AB, Results.error_QY_BA,
                           ExpParams.CalculationMethod,
                           SaveResults = "Yes",
                           SaveFileName = Results.savefilename)


        self.Save_Results()

        QtWidgets.QMessageBox.information(self, "Success", f"{Results.savefilename} file saved successfully!")

    def Save_Results(self):
        """ Save results: """
        savefile = Results.savefilename+".txt"

        dict_results = {'PSS_Reactant (%)': Results.PSS_Reactant,
                 'PSS_Product (%)' : Results.PSS_Product,
                 'QY_AB_opt' : Results.QY_AB_opt,
                 'QY_BA_opt' : Results.QY_BA_opt,
                 'error_QY_AB' : Results.error_QY_AB,
                 'error_QY_BA' : Results.error_QY_BA}

        dict_expparams = {'Volume (ml)': ExpParams.V,
                          'k thermal back-reaction (s-1)': ExpParams.k_BA,
                          'Power average (mW)': ExpParams.I0_avg,
                          'Power error (mW)': ExpParams.I0_err,
                          'Calculation Method': ExpParams.CalculationMethod,
                          'Wavelength of irradiation': ExpParams.LEDw,
                          'Threshold': ExpParams.threshold}

        try:
            os.remove(savefile)
        except:
            pass
        # except OSError as e:
        #     print(f"An error occurred for {e.filename} - {e.strerror}")

        try:
            file=savefile
            
            with open (file,'a') as file:
                for i in dict_results:
                    file.write(i+": "+str(dict_results[i])+'\n')
                # file.write('==== Experimental Parameters ====\n')
                for i in dict_expparams:
                    file.write(i+": "+str(dict_expparams[i])+'\n')
        except IOError as e:
            print(f"An error occurred: {e}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = autoQuant()
    window.show()
    sys.exit(app.exec_())
