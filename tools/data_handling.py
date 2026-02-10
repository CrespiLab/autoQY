# -*- coding: utf-8 -*-
import os
import csv
import data.experimental_parameters as ExpParams
import data.loaded_data as LoadedData
import data.results as Results
import data.calc_settings as CalcSettings
import data.datasets as Datasets
from tools.plotting import MplCanvas

class FileHandler:
    """ Class for handling files (e.g., saving and loading)"""
    def __init__(self, filename, filetype, parent):
        self.parent = parent
        super(FileHandler, self).__init__()
        self.i = '' ## iteration of unique filename
        
        try:
            if filetype == "save":
                self.filename_user = filename ## filename given by user
                self.filename_mod = self.modify_filename(self.filename_user, "Results") ## filename modified for Results
                self.filename = self._get_unique_filename(self.filename_mod) ## unique Results filename
                self.obtain_filename_data() ## define filenames of data used in calculation
            elif filetype == "load":
                self.filename = filename ## load file with known filename
        except:
            print("opening file was unsuccessful (unknown error)")

    def modify_filename(self, filename, a, b = ''):
        """ 
        Modify filename:
            - add {a} to start of name
            - add {b} to end of name (if given)
        """
        name, ext = os.path.splitext(filename)
        
        path_split=name.split('/')[0:-1] # leave only filepath (remove name)
        path='\\'.join(path_split) # re-join into string
        
        end_nameonly=name.split('/')[-1] # only filename
        if b != '':
            end = f"{a}_{end_nameonly}_{b}" # name with added info
        else:
            end = f"{a}_{end_nameonly}" # name with added info
        
        name_mod = f"{path}\\{end}"
        return name_mod

    def _get_unique_filename(self, base_filename):
        """
        If 'filename.txt' exists, make 'filename_2.txt', etc.
        """
        ext='.txt'
        if not os.path.exists(f"{base_filename}{ext}"): ## search for filename.txt
            return base_filename

        i = 2
        while os.path.exists(f"{base_filename}_{i}"):
            i += 1
        self.i = i
        return f"{base_filename}_{i}"
    
    def obtain_filename_data(self):
        self.filename_spectra = LoadedData.filename_spectra ## full filepath of measurement data
        self.filename_epsilons_R = LoadedData.epsilons_R
        self.filename_epsilons_P = LoadedData.epsilons_P
        self.filename_LED = LoadedData.filename_LED
        self.filename_log = LoadedData.filename_log
    
    def build_dicts_results(self):
        """ Build the dictionaries that will be saved into the Results textfile """
        
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

        dict_calcsettings = {'Calculation Method': CalcSettings.CalculationMethod,
                             'ODE Solving Method': CalcSettings.ODEMethod,
                             # 'Baseline Correction LED Spectrum': CalcSettings.BaselineCorrection_LED ##!!! add when general option
                             # 'Smoothing of LED Spectrum': CalcSettings.Smoothing_LED ##!!! add when general option
                             }

        if CalcSettings.ODEMethod == "Emission":
            dict_calcsettings['Threshold'] = CalcSettings.threshold

        elif CalcSettings.ODEMethod == "Concentrations":
            dict_calcsettings['Baseline Correction LED Spectrum'] = CalcSettings.BaselineCorrection_LED
            dict_calcsettings['Wavelength Range'] = f"{CalcSettings.wl_low}-{CalcSettings.wl_high}"

        dict_data = {'Measurement Data': self.filename_spectra,
                     'Epsilons Reactant': self.filename_epsilons_R,
                     'Epsilons Product': self.filename_epsilons_P,
                     'LED Emission': self.filename_LED,
                     'Log (timestamps)': self.filename_log
                     }

        return dict_results, dict_expparams, dict_calcsettings, dict_data

    def write_to_textfile_results(self):
        """ Write Results, Experimental Parameters, and Calculation Settings to Results textfile """
        file = f"{self.filename}.txt"
        dict_results, dict_expparams, dict_calcsettings, dict_data = self.build_dicts_results()
         
        with open (file,'a') as file:
            file.write('Results Obtained and Parameters Used'+'\n')
            file.write('\n')
            
            for i in dict_results:
                file.write(i+": "+str(dict_results[i])+'\n')
            file.write('\n')
            
            for i in dict_expparams:
                file.write(i+": "+str(dict_expparams[i])+'\n')
            file.write('\n')
            
            for i in dict_calcsettings:
                file.write(i+": "+str(dict_calcsettings[i])+'\n')
            file.write('\n')
            
            for i in dict_data:
                file.write(i+": "+str(dict_data[i])+'\n')

    def save_plots_results(self):
        canvas = MplCanvas(self.parent)
        
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
                               SaveFileName = self.filename)
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
                               SaveFileName = self.filename)

        self.parent.message_console.append(f"Plots and textfile of results saved successfully as {self.filename}")

    def Save_FractionsResults(self):
        ''' 
        Save Fractions calculation results.
        Using the same name as that of the Measurement Data (spectra measured during irradiation)
        '''
        ### Output results ###
        # for i, (f1, f2) in enumerate(zip(f1_list, f2_list)):
        #     print(f"{i+1}   f₁ = {f1:.4f}   f₂ = {f2:.4f}")
        ##!!! OUTPUT RESULTS IN SOME WAY IN INFO SCREEN

        f1_list, f2_list = Results.fractions_R, Results.fractions_P

        filename_fractions = self.modify_filename(self.filename_user, "Fractions", self.i)
        output_file = f"{filename_fractions}.csv"

        	### Write to CSV ###
        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter=",")        
            writer.writerow(["Spectrum", "Fraction_Reactant", "Fraction_Product"])
            for i, (f1, f2) in enumerate(zip(f1_list, f2_list)):
                writer.writerow([i+1, f1, f2])
    
        self.parent.message_console.append(f"Fractions saved to: {output_file}")
        return
