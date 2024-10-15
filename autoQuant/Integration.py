"""
@author: Jorn
Date started: Sep 16th 2024
Based on: script by Anouk Volker from our publication

Programme with GUI to calculate the QYs from the absorption data and the 
emission spectrum of the LED

A <=> B
Starting from species A (100%)

##########
##!!!
TO DO:
[] Create functions for all the code
[] Add descriptions to functions

[] include option to set starting percentage (see below; StartPercentage_A as optional parameter)

##########
GUI TO DO:
[] include a user-selection option: data is with or without Wavenumbers column; for all datasets
[] add note/message in GUI: use raw (nanometer) data of LED emission

"""
import os
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid,odeint
from lmfit import minimize, Parameters
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec

#import Constants
#import ExpParam

k_BA = 7.240e-7  ####################N.B.: COMMENTED ExpParam
h = 6.626e-34           # Planck's constant in J·s
c = 299792458           # Speed of light in m/s
Avogadro = 6.022e23     # Avogadro's number

#%% DATA AND PARAMETERS
# ## FOR AUTOQUANT.PY

# script_dir = os.path.dirname(os.path.abspath(__file__)) ## Directory where the script is located
# datafolder = script_dir+r"\ExampleData"
# ###################################################
# ############# FILE FORMAT #############
# # FileFormat = "Spectragryph"
# FileFormat = "Not"
# ##############################

# ############# EPSILONS #############
# if FileFormat == "Spectragryph":
#     epsilonsfile_A = datafolder+"\\"+"trans_azobenzene_epsilons.dat"
#     epsilonsfile_B = datafolder+"\\"+"cis_azobenzene_epsilons.dat"
# elif FileFormat == "Not":
#     epsilonsfile_A = datafolder+"\\"+"trans_azobenzene_epsilons.csv"
#     epsilonsfile_B = datafolder+"\\"+"cis_azobenzene_epsilons.csv"
# ###################################################

# ############################################
# ############# IRRADIATION DATA #############
# datafilename = "Azobenzene_340nm_processed"

# data_file = datafolder+"\\"+datafilename+".dat" 
# log_file = datafolder+"\\"+"Azobenzene_340nm_log.csv"
# ###################################################

# ########################################
# ############# LED EMISSION #############
# #### Input LED emission data in nanometers: just the raw data
# ##=> it gets normalised during fitting!


# LED_file = datafolder+"\\"+"LED340nm_normalisedarea.csv"
# NominalWavelengthLED = 340    # Wavelength of interest
# ###################################################

# #######################
# ######## POWER ########
# I0_avg = 743             # Photon flux in microWatt
# I0_err = 4                # Error on photon fplux in microWatt
# I0_list=[I0_avg, I0_avg+I0_err, I0_avg-I0_err]
# ###################################################

# ############# EXPERIMENTAL PARAMETERS #############
# V = 3.0                 # Volume in ml
# k_BA = 7.240e-7         # Thermal back reaction rate s-1
# ###################################################

# ###################################################
# ################## CONSTANTS ######################
# h = 6.626e-34           # Planck's constant in J·s
# c = 299792458           # Speed of light in m/s
# Avogadro = 6.022e23     # Avogadro's number
# ###################################################

#%% CREATE DATASETS

##################################################
################### TIMESTAMPS ###################
##################################################
## FOR AUTOQUANT.PY

# #### Import log file of irradiation data
# def GetTimestamps(LogFile):
#     log = pd.read_csv(LogFile,
#                       sep = ",", decimal = ".", skiprows = 1, header=None,)
#     log_t=log[log[3] == 'Measure']
#     t=log_t[2]
#     timestamps=t.to_numpy()
#     return timestamps
# timestamps = GetTimestamps(log_file)

#%% LED 

############# ANALYSIS PARAMETERS #################
# threshold_LED = 500     # Threshold for part of spectrum where there is LED emission 
##!!! make this a user-defined option

##################################################
###################### LED #######################
##################################################
def Import_LEDemission(FileFormat, file_LEDemission_raw):
    ## Import (raw) emission LED_file
    if FileFormat == "Spectragryph":
        emission_data = pd.read_csv(file_LEDemission_raw, sep = '\t', usecols = lambda x: x not in ["Wavenumbers [1/cm]"]) # in nanometers
    elif FileFormat == "Not":
        emission_data = pd.read_csv(file_LEDemission_raw, delimiter=',') # in meters #MISSING DELIMITER ALFREDO
    
    #ALFREDO: I think this is not necessary
    # Clean up the column names by stripping leading/trailing spaces and invisible characters
    emission_data.columns = emission_data.columns.str.strip()

    emission_wavelengths = emission_data['Wavelength [nm]'].values
    emission_Intensity = emission_data['Intensity'].values ## not normalised
    return emission_wavelengths, emission_Intensity
    ####################################

####################################
def Processing_LEDemission(wavelengths_LED, intensity_LED, threshold):
    
    emission_Intensity_proc = savgol_filter(intensity_LED, 12,3) ## Smoothing
    emission_Intensity_proc[emission_Intensity_proc[:]<0] = 0 ## removal of zeroes 

    #### Cutting spectrum to match emission of LED 
    above_threshold_indices = np.where(emission_Intensity_proc > threshold)[0] 
    
    #### Find the first and last indices
    first_index = above_threshold_indices[0]
    last_index = above_threshold_indices[-1]
    
    low = wavelengths_LED[first_index]
    high = wavelengths_LED[last_index]
    #print(f"wavelength_low: {low}\nwavelength_high: {high}")
    
    return emission_Intensity_proc, first_index, last_index, low, high
    ####################################

##################################################
################### ABSORBANCE ###################
##################################################

def Import_SpectralData(FileFormat, file, wl_low, wl_high, wavelength_of_interest):
    ## Import Absorbance data_file (.dat from Spectragryph, so with Wavenumbers column)
    if FileFormat == "Spectragryph":
        data_pd = pd.read_csv(file, sep='\t', usecols=lambda x: x not in ["Wavenumbers [1/cm]"])  # in nanometers
    elif FileFormat == "Not":
        data_pd = pd.read_csv(file, delimiter=',', header=0)  # in nanometers ####ALFREDO: FIXED?

    # Clean up the first column by splitting space-separated values and selecting the first one
    def clean_column(value):
        if isinstance(value, str):
            parts = value.split()  # Split the string by spaces
            return float(parts[0])  # Take the first part and convert to float
        return value  # If it's already a number, just return it

    data_pd.iloc[:, 0] = data_pd.iloc[:, 0].apply(clean_column)
    data = data_pd.to_numpy()  # Convert DataFrame to Numpy array

    # Perform operations on the cleaned numeric data
    #(data)
    low = np.argmin(np.abs(data[:, 0].astype(float) - wl_low))
    high = np.argmin(np.abs(data[:, 0].astype(float) - wl_high))

    wavelengths_lowhigh = data[low:high, 0]
    #print(data[low:high, :].shape)
    absorbance_lowhigh = data[low:high, 1:]
    #print('shape : ', data[low:high, 1:].shape)

    index = np.abs(wavelengths_lowhigh - wavelength_of_interest).argmin()

    #print(wavelengths_lowhigh)
    #print(absorbance_lowhigh)

    return wavelengths_lowhigh, absorbance_lowhigh, index

##################################################
###################### PLOT ######################
##################################################

def Plot_LEDemission_Processed_(absorbance, wavelengths,
                               em_wl, em_int,
                               index_first,index_last,
                               em_int_proc):
    fig, axdata = plt.subplots(2,1,figsize=(6,6),
                                dpi=600,constrained_layout=True)
    
    for i in range(0,absorbance.shape[1]):
        axdata[0].plot(wavelengths, absorbance.T[i])
    
    for i in axdata:
        i.set_xlabel("Wavelength (nm)")
        i.set_xlim(220,650)
    
    axdata[0].set_title("Spectral Data")
    axdata[0].set_ylabel("Absorbance")
    axdata[0].set_ylim(-0.05,2.5)
    
    axdata[1].plot(em_wl,em_int,
                    label="Untreated (in this .py file at least)")
    axdata[1].plot(em_wl[index_first:index_last],
                    em_int_proc[index_first:index_last],
                    label="Smoothed, removed zeroes\nand applied threshold")
    axdata[1].legend(fontsize=8)
    axdata[1].set_title("LED emission")
    axdata[1].set_ylabel("Intensity")
    plt.show()
    ####################################

####################################
## FOR AUTOQUANT.PY 
# emission_wavelengths, emission_Intensity = Import_LEDemission(LED_file)
# emission_Intensity_proc, LEDindex_first, LEDindex_last, wavelength_low, wavelength_high = Processing_LEDemission(emission_wavelengths, emission_Intensity, threshold_LED)
# SpectralData_Wavelengths, SpectralData_Abs, SpectralData_Index = Import_SpectralData(data_file, wavelength_low, wavelength_high, NominalWavelengthLED)
# Plot_LEDemission_Processed()
####################################
#%% EPSILONS LOAD
""" 
Epsilons datafiles are in a certain format
TO DO for GUI:
    - select delimiter: ',' or '\t'
    - select; ignore Wavenumbers column or not
    
"""

##################N.B.: COMMENTED TO MERGE WITH MAIN.PY!
def Import_Epsilons(FileFormat, 
                    A): #NO B - ALFREDO
    if FileFormat == "Spectragryph":
        epsilon_A_data = pd.read_csv(A, delimiter='\t',
                                     usecols = lambda x: x not in ["Wavenumbers [1/cm]"])
    #    epsilon_B_data= pd.read_csv(B, delimiter='\t', 
    #                                usecols = lambda x: x not in ["Wavenumbers [1/cm]"])
    elif FileFormat == "Not":
        epsilon_A_data = pd.read_csv(A, delimiter=',')
    #    epsilon_B_data = pd.read_csv(B, delimiter=',')
        
    epsilon_A_data.columns = ['Wavelengths', 'Epsilons'] ## rename columns
    epsilon_A_wavelengths = epsilon_A_data['Wavelengths'].values
    epsilon_A_values = epsilon_A_data['Epsilons'].values
    
    #epsilon_B_data.columns = ['Wavelengths', 'Epsilons'] ## rename columns
    #epsilon_B_wavelengths = epsilon_B_data['Wavelengths'].values
    #epsilon_B_values = epsilon_B_data['Epsilons'].values
    return epsilon_A_wavelengths, epsilon_A_values#, epsilon_B_wavelengths, epsilon_B_values
##!!! reduce the number of returned objects

def Interpolate_Epsilons(Abs_wavelengths,
                         eA_wavelengths, eA_values,
                         eB_wavelengths, eB_values,
                         LED_wavelengths, LED_intensity_proc):
    ## Interpolate ε_A and ε_B at each wavelength so that all data sets use the same 
    ## wavelengths data
    epsilon_A_interp = np.interp(Abs_wavelengths, eA_wavelengths, 
                                 eA_values)
    epsilon_B_interp = np.interp(Abs_wavelengths, eB_wavelengths, 
                                 eB_values)
    emission_interp = np.interp(Abs_wavelengths, LED_wavelengths,
                                LED_intensity_proc)
    return epsilon_A_interp, epsilon_B_interp, emission_interp

##################################################
###################### PLOT ######################
##################################################

def Plot_Epsilons(wavelengths, 
                  e_A_inter, e_B_inter,
                  emission_inter):
    fig, axdata_interp = plt.subplots(2,1,figsize=(6,6),
                                dpi=600,constrained_layout=True)
    
    axdata_interp[0].plot(wavelengths, e_A_inter, label="A")
    axdata_interp[0].plot(wavelengths, e_B_inter, label="B")
    axdata_interp[0].legend(fontsize=12)
    
    for i in axdata_interp:
        i.set_xlabel("Wavelength (nm)")
        i.set_xlim(220,650)
    
    axdata_interp[0].set_title("Epsilons, interpolated")
    axdata_interp[0].set_ylabel(r"$\epsilon$ (M$^{-1}$ cm$^{-1}$)")
    # axdata_interp[0].set_ylim(-0.05,2.5)
    
    axdata_interp[1].set_title("LED emission, interpolated")
    axdata_interp[1].set_ylabel("Intensity")
    axdata_interp[1].plot(wavelengths,emission_inter)
    # axdata_interp[1].legend(fontsize=8)
    
    plt.show()

##################################################
## FOR AUTOQUANT.PY

# epsilon_A_wavelengths, epsilon_A_values, epsilon_B_wavelengths, epsilon_B_values = Import_Epsilons(epsilonsfile_A, epsilonsfile_B)

# epsilon_A_interp, epsilon_B_interp, emission_interp = Interpolate_Epsilons(SpectralData_Wavelengths,
#                      epsilon_A_wavelengths, epsilon_A_values,
#                      epsilon_B_wavelengths, epsilon_B_values,
#                      emission_wavelengths, emission_Intensity_proc)

# Plot_Epsilons()

##################################################
#%% SOLVE RATE EQUATIONS AND PLOT RESULTS
# StartPercentage_A = float(input("Starting Percentage of A: "))
# StartPercentage_A = float(100)

def CreateParameters(absorbance_values, wavelengths_data,
                     e_A_inter, emission_inter):
    #########################################
    ##### trying out: input starting percentages ######
    #########################################
    StartPercentage_A = float(100) ##!!! turn into optional parameter

    initial_conc_A_100 = trapezoid(absorbance_values[:,0],
                                   x=wavelengths_data) / trapezoid(e_A_inter,
                                                                   x=wavelengths_data)
    initial_conc_A = initial_conc_A_100/100*StartPercentage_A
    initial_conc_B_0 = 0
    initial_conc_B = initial_conc_B_0/100*(100-StartPercentage_A)
#########################################
#########################################
    # epsilon_A = epsilon_A_interp            # array of epsilon_A values
    # epsilon_B = epsilon_B_interp            # array of epsilon_B values
    total_absorbance = absorbance_values.T  # array of total_absorbance values
    # time = timestamps

    lambda_meters = wavelengths_data * 1e-9  ## Convert to meters

    ## Normalize the emission spectrum to ensure the area under the curve is 1
    normalized_emission = emission_inter / trapezoid(emission_inter, lambda_meters)

    return initial_conc_A, initial_conc_B, total_absorbance, lambda_meters, normalized_emission

def rate_equations(concentrations, time, 
                   lambda_meters,
                   QY_AB, QY_BA, 
                   epsilon_A_lambda, epsilon_B_lambda, 
                   N_lambda, V):
    """
    Rate equations for the isomerization, as described in 
    10.1038/srep41145 (Stranius & Börjesson)
    and
    H. Volfova, Q. Hu, E. Riedle, EPA Newsl. 2019.

    Parameters
    ----------
    concentrations : tuple
        Concentrations of A and B
    time : numpy array
        Time 
    QY_AB : float
        Quantum Yield of A --> B
    QY_BA : float
        Quantum Yield of B --> A

    ##!!! ELABORATE

    Returns
    -------
    [dAdt, dBdt] : list
        List with the functions for the changes in concentration over time

    """
    A, B = concentrations
    total_absorbance = A * epsilon_A_lambda + B * epsilon_B_lambda

    dAdt = QY_AB * trapezoid(-((A * epsilon_A_lambda / total_absorbance) \
            * (1 - 10 ** (-total_absorbance)) * N_lambda) / V,lambda_meters)+ \
            QY_BA * trapezoid(((B * epsilon_B_lambda / total_absorbance) \
            * (1 - 10 ** (-total_absorbance)) * N_lambda) / V,lambda_meters) \
            + k_BA * B

    dBdt = -dAdt

    return [dAdt, dBdt]
#########################################
# rate_equations()
#########################################

def objective_function(params, lambda_meters, 
                       init_conc_A, init_conc_B,
                       time, experimental_data, 
                       epsilon_A_lambda, epsilon_B_lambda, 
                       N_lambda, V):
    QY_AB = params['QY_AB'].value
    QY_BA = params['QY_BA'].value

    # print(f"lambda_meters: {lambda_meters}")    
    # print(f"initial concentrations: {init_conc_A}, {init_conc_B}")
    # print(f"time: {time}")
    # print(f"volume: {V}")

    concentrations_ode = odeint(rate_equations, [init_conc_A, init_conc_B], time,
                               args=(lambda_meters, QY_AB, QY_BA, epsilon_A_lambda, 
                                     epsilon_B_lambda, N_lambda, V))
    
    # Matrix multiplication using dot product
    total_absorbance_ode=concentrations_ode.dot(np.vstack([epsilon_A_lambda,
                                                            epsilon_B_lambda]))

    return total_absorbance_ode - experimental_data.T

def MinimizeQYs(I0_list,
                emission_norm, lambda_meters, 
                init_conc_A, init_conc_B,
                timestamps, absorbance_values,
                e_A_inter, e_B_inter,
                V):
    
    # print(f"lambda_meters: {lambda_meters}")
    # print(f"initial concentrations: {init_conc_A}, {init_conc_B}")
    # print(f"volume: {V}")
    # print(f"k_BA: {k_BA}")
    # print(f"timestamps: {timestamps}")
    
    
    ### loop over the powers in the list 
    fit_results=[]
    for total_power in I0_list:
        power_at_each_wavelength = emission_norm * total_power * 1e-6 ## Photon flux in mmol/s per wavelength
        N = power_at_each_wavelength / (h * c / lambda_meters) / Avogadro * 1000  
        
        ## Add wavelength-specific parameters to the lmfit Parameters object
        params = Parameters()
        params.add('QY_AB', value=0.5, min=0, max=1) ## with boundaries
        params.add('QY_BA', value=0.5, min=0, max=1) ## with boundaries
        
        ## Minimize the objective function using lmfit   
        result_lmfit = minimize(objective_function, 
                                params, args=(lambda_meters,
                                              init_conc_A, init_conc_B,
                                              timestamps, absorbance_values, 
                                              e_A_inter, e_B_inter, 
                                              N, V))
        # print(f"Results for I0={total_power} mW: ")
        # print(result_lmfit.params)
        # print("\n")
        fit_results.append(result_lmfit)
    return N, fit_results

#########################################################
## FOR AUTOQUANT.PY
# initial_conc_A, initial_conc_B, total_absorbance, lambda_meters, normalized_emission = CreateParameters(SpectralData_Abs,
#                                                                                          SpectralData_Wavelengths)

# N, fit_results = MinimizeQYs(SpectralData_Abs)

#########################################################

#%% Extract results, plot, and save
def ExtractResults(fit_results):
    result_lmfit=fit_results[0] ## results using I0_avg
    ## Extract the optimized parameters and their standard deviations
    QY_AB_opt = result_lmfit.params['QY_AB'].value ## value using I0_avg
    std_dev_fit_QY_AB = result_lmfit.params['QY_AB'].stderr ## error in the fit using I0_avg
    
    QY_BA_opt = result_lmfit.params['QY_BA'].value ## value using I0_avg
    std_dev_fit_QY_BA = result_lmfit.params['QY_BA'].stderr ## error in the fit using I0_avg
        
    QY_AB_opt_min=fit_results[1].params['QY_AB'].value \
                        - fit_results[1].params['QY_AB'].stderr ## results using I0_avg+I0_err
    QY_AB_opt_max=fit_results[2].params['QY_AB'].value \
                        - fit_results[2].params['QY_AB'].stderr ## results using I0_avg-I0_err
    
    QY_BA_opt_min=fit_results[1].params['QY_BA'].value \
                        - fit_results[1].params['QY_BA'].stderr ## results using I0_avg+I0_err
    QY_BA_opt_max=fit_results[2].params['QY_BA'].value \
                        - fit_results[2].params['QY_BA'].stderr ## results using I0_avg-I0_err
    
    ## Calculate error for QY_AB and QY_BA based on the error in the power 
    ## So: the error in the fit is not considered, since it is smaller
    error_QY_AB = max([QY_AB_opt-QY_AB_opt_min, QY_AB_opt_max-QY_AB_opt])
    error_QY_BA = max([QY_BA_opt-QY_BA_opt_min, QY_BA_opt_max-QY_BA_opt])
    
    ## Print the results
    print(f"Optimized QY_AB: {QY_AB_opt:.5f}" )
    print(f"Standard deviation of QY_AB in the fit using I0_avg: {std_dev_fit_QY_AB:.5f}")
    print(f"Error for QY_AB: {error_QY_AB:.5f} ")
    
    print(f"Optimized QY_BA: {QY_BA_opt:.5f}" )
    print(f"Standard deviation of QY_BA in the fit using I0_avg: {std_dev_fit_QY_BA:.5f}")
    print(f"Error for QY_BA: {error_QY_BA:.5f} ")
    
    return QY_AB_opt, QY_BA_opt, error_QY_AB, error_QY_BA

##################################################
def CalculateConcentrations(lambda_meters, 
                            init_conc_A, init_conc_B,
                            timestamps,
                            QY_AB_opt, QY_BA_opt,
                            e_A_inter, e_B_inter,
                            N, V):
    ## Integrate the rate equations with the optimized parameters
    
    conc_opt = odeint(rate_equations, [init_conc_A, init_conc_B], timestamps,
                               args=(lambda_meters, QY_AB_opt, QY_BA_opt, 
                                     e_A_inter, e_B_inter,
                                     N, V))
    
    
    PSS_A = conc_opt[-1,0]/(init_conc_A+init_conc_B)*100
    PSS_B = conc_opt[-1,1]/(init_conc_A+init_conc_B)*100
    
    print(f"At the PSS, {PSS_A:.2f} % of A and {PSS_B:.2f} % of B")
    return conc_opt
##################################################
def GetFittedAbs(fit_results, conc_opt,
                 e_A_inter, e_B_inter,
                 timestamps,
                 wavelengths_data):
    ##### Get the fitted absorbance
    result_lmfit=fit_results[0] ## results using I0_avg
    total_abs_fit = conc_opt.dot(np.vstack([e_A_inter, e_B_inter]))
    residuals = result_lmfit.residual.reshape((len(timestamps), 
                                              len(wavelengths_data))).T
    return total_abs_fit, residuals
##################################################
def PlotAndSave(LEDwavelength, datafolder, filename,
                timestamps, conc_opt,
                absorbance, index, 
                total_abs_fit, residuals,
                QY_AB_opt, QY_BA_opt, 
                error_QY_AB, error_QY_BA):
    # ToSave = input("Save plots (Yes or No)? ")
    ToSave = "No" ## for testing
    ##################################################
    
    ######################
    fig = plt.figure(figsize=(8, 4),dpi=600,constrained_layout=True)
    gs = gridspec.GridSpec(4, 2, figure=fig)
    fig.suptitle(f'{LEDwavelength} nm: Integration Method\n \
                  {filename}')
    
    axresults_conc = fig.add_subplot(gs[0:3, 0])
    axresults_Abs = fig.add_subplot(gs[0:3, 1])
    axresults_res = fig.add_subplot(gs[3, 1])
    axresults_notes = fig.add_subplot(gs[3, 0])
    
    axresults_conc.set_title("Concentrations over time")
    
    axresults_conc.plot(timestamps, conc_opt[:,0], label="Species A")
    axresults_conc.plot(timestamps, conc_opt[:,1], label = "Species B")
    
    axresults_conc.set_ylabel("Concentration (M)")
    axresults_conc.yaxis.set_minor_locator(AutoMinorLocator(2))
    axresults_conc.set_xlabel("Time (s)")
    
    axresults_conc.legend()
    
    ##################################################
    ################ NOTES ################
    axresults_notes.axis("off")
    # axresults_notes.text(0, 0, f"Starting Percentage A: {StartPercentage_A} %")
    ##################################################
    #### Plot the experimental data and the fitted total absorbance curve
    axresults_Abs.plot(timestamps, absorbance[index,:], linestyle='-', 
                color="#DCEEFF", linewidth=5, label='Experimental Data')
    axresults_Abs.plot(timestamps, total_abs_fit[:,index], linestyle='--', color="#9ACCFF", 
                label=f"Fitted Total Absorbance Curve\nQY_AB: {QY_AB_opt:.3f}"+u"\u00B1"+f"{error_QY_AB:.3f}\
                \nQY_BA: {QY_BA_opt:.3f}"+u"\u00B1"+f"{error_QY_BA:.3f}")
    axresults_Abs.set_title('Absorbance over time')
    axresults_Abs.set_xlabel('Time (s)')
    axresults_Abs.set_ylabel('Absorbance')
    axresults_Abs.legend()
    
    #### Plot the residuals over time
    axresults_res.plot(timestamps, residuals[index,:], color="#FF952A", label='Residuals')
    axresults_res.set_xlabel('Time (s)')
    axresults_res.set_ylabel('Residual Abs')
    
    #### SAVE PLOTS ####
    if ToSave == "Yes":
        savefilename = datafolder+"\\"+"QY-Integrated_"+str(LEDwavelength)+"nm_"+filename
        plt.savefig(savefilename+".svg",bbox_inches="tight")
        plt.savefig(savefilename+".png",bbox_inches="tight")
        print("Saved plots")
    elif ToSave == "No":
        print("Plots not saved")
    else:
        print("Wrong ToSave input")
    ####################
    plt.show()
#################################################
## FOR AUTOQUANT.PY
# QY_AB_opt, QY_BA_opt, error_QY_AB, error_QY_BA = ExtractResults()
# conc_opt = CalculateConcentrations()
# total_abs_fit, residuals = GetFittedAbs()
# PlotAndSave()
#################################################
#%%

def CreateParameters(absorbance_values, wavelengths_data,
                     e_A_inter, emission_inter):
    #########################################
    ##### trying out: input starting percentages ######
    #########################################
    StartPercentage_A = float(100) ##!!! turn into optional parameter

    initial_conc_A_100 = trapezoid(absorbance_values[:,0],
                                   x=wavelengths_data) / trapezoid(e_A_inter,
                                                                   x=wavelengths_data)
    initial_conc_A = initial_conc_A_100/100*StartPercentage_A
    initial_conc_B_0 = 0
    initial_conc_B = initial_conc_B_0/100*(100-StartPercentage_A)
#########################################
#########################################
    # epsilon_A = epsilon_A_interp            # array of epsilon_A values
    # epsilon_B = epsilon_B_interp            # array of epsilon_B values
    total_absorbance = absorbance_values.T  # array of total_absorbance values
    # time = timestamps

    lambda_meters = wavelengths_data * 1e-9  ## Convert to meters

    ## Normalize the emission spectrum to ensure the area under the curve is 1
    normalized_emission = emission_inter / trapezoid(emission_inter, lambda_meters)

    return initial_conc_A, initial_conc_B, total_absorbance, lambda_meters, normalized_emission

def rate_equations(concentrations, time, 
                   lambda_meters,
                   QY_AB, QY_BA, 
                   epsilon_A_lambda, epsilon_B_lambda, 
                   N_lambda, V):
    """
    Rate equations for the isomerization, as described in 
    10.1038/srep41145 (Stranius & Börjesson)
    and
    H. Volfova, Q. Hu, E. Riedle, EPA Newsl. 2019.

    Parameters
    ----------
    concentrations : tuple
        Concentrations of A and B
    time : numpy array
        Time 
    QY_AB : float
        Quantum Yield of A --> B
    QY_BA : float
        Quantum Yield of B --> A

    ##!!! ELABORATE

    Returns
    -------
    [dAdt, dBdt] : list
        List with the functions for the changes in concentration over time

    """
    A, B = concentrations
    total_absorbance = A * epsilon_A_lambda + B * epsilon_B_lambda

    dAdt = QY_AB * trapezoid(-((A * epsilon_A_lambda / total_absorbance) \
            * (1 - 10 ** (-total_absorbance)) * N_lambda) / V,lambda_meters)+ \
            QY_BA * trapezoid(((B * epsilon_B_lambda / total_absorbance) \
            * (1 - 10 ** (-total_absorbance)) * N_lambda) / V,lambda_meters) \
            + k_BA * B

    dBdt = -dAdt

    return [dAdt, dBdt]
#########################################
# rate_equations()
#########################################

def objective_function(params, lambda_meters, 
                       init_conc_A, init_conc_B,
                       time, experimental_data, 
                       epsilon_A_lambda, epsilon_B_lambda, 
                       N_lambda, V):
    QY_AB = params['QY_AB'].value
    QY_BA = params['QY_BA'].value

    #print(f"lambda_meters: {lambda_meters}")    
    #print(f"initial concentrations: {init_conc_A}, {init_conc_B}")
    #print(f"time: {time}")
    #print(f"volume: {V}")

    concentrations_ode = odeint(rate_equations, [init_conc_A, init_conc_B], time,
                               args=(lambda_meters, QY_AB, QY_BA, epsilon_A_lambda, 
                                     epsilon_B_lambda, N_lambda, V))
    
    # Matrix multiplication using dot product
    total_absorbance_ode=concentrations_ode.dot(np.vstack([epsilon_A_lambda,
                                                            epsilon_B_lambda]))

    return total_absorbance_ode - experimental_data.T

def MinimizeQYs(I0_list,
                emission_norm, lambda_meters, 
                init_conc_A, init_conc_B,
                timestamps, absorbance_values,
                e_A_inter, e_B_inter,
                V):
    
    # print(f"lambda_meters: {lambda_meters}")
    # print(f"initial concentrations: {init_conc_A}, {init_conc_B}")
    # print(f"volume: {V}")
    # print(f"k_BA: {k_BA}")
    # print(f"timestamps: {timestamps}")
    
    
    ### loop over the powers in the list 
    fit_results=[]
    for total_power in I0_list:
        power_at_each_wavelength = emission_norm * total_power * 1e-6 ## Photon flux in mmol/s per wavelength
        N = power_at_each_wavelength / (h * c / lambda_meters) / Avogadro * 1000  
        
        ## Add wavelength-specific parameters to the lmfit Parameters object
        params = Parameters()
        params.add('QY_AB', value=0.5, min=0, max=1) ## with boundaries
        params.add('QY_BA', value=0.5, min=0, max=1) ## with boundaries
        
        ## Minimize the objective function using lmfit   
        result_lmfit = minimize(objective_function, 
                                params, args=(lambda_meters,
                                              init_conc_A, init_conc_B,
                                              timestamps, absorbance_values, 
                                              e_A_inter, e_B_inter, 
                                              N, V))
        #print(f"Results for I0={total_power} mW: ")
        #print(result_lmfit.params)
        #print("\n")
        fit_results.append(result_lmfit)
    return N, fit_results

#########################################################
## FOR AUTOQUANT.PY
# initial_conc_A, initial_conc_B, total_absorbance, lambda_meters, normalized_emission = CreateParameters(SpectralData_Abs,
#                                                                                          SpectralData_Wavelengths)

# N, fit_results = MinimizeQYs(SpectralData_Abs)

#########################################################

#%% Extract results, plot, and save
def ExtractResults(fit_results):
    result_lmfit=fit_results[0] ## results using I0_avg
    ## Extract the optimized parameters and their standard deviations
    QY_AB_opt = result_lmfit.params['QY_AB'].value ## value using I0_avg
    std_dev_fit_QY_AB = result_lmfit.params['QY_AB'].stderr ## error in the fit using I0_avg
    
    QY_BA_opt = result_lmfit.params['QY_BA'].value ## value using I0_avg
    std_dev_fit_QY_BA = result_lmfit.params['QY_BA'].stderr ## error in the fit using I0_avg
        
    QY_AB_opt_min=fit_results[1].params['QY_AB'].value \
                        - fit_results[1].params['QY_AB'].stderr ## results using I0_avg+I0_err
    QY_AB_opt_max=fit_results[2].params['QY_AB'].value \
                        - fit_results[2].params['QY_AB'].stderr ## results using I0_avg-I0_err
    
    QY_BA_opt_min=fit_results[1].params['QY_BA'].value \
                        - fit_results[1].params['QY_BA'].stderr ## results using I0_avg+I0_err
    QY_BA_opt_max=fit_results[2].params['QY_BA'].value \
                        - fit_results[2].params['QY_BA'].stderr ## results using I0_avg-I0_err
    
    ## Calculate error for QY_AB and QY_BA based on the error in the power 
    ## So: the error in the fit is not considered, since it is smaller
    error_QY_AB = max([QY_AB_opt-QY_AB_opt_min, QY_AB_opt_max-QY_AB_opt])
    error_QY_BA = max([QY_BA_opt-QY_BA_opt_min, QY_BA_opt_max-QY_BA_opt])
    
    ## Print the results
    print(f"Optimized QY_AB: {QY_AB_opt:.5f}" )
    print(f"Standard deviation of QY_AB in the fit using I0_avg: {std_dev_fit_QY_AB:.5f}")
    print(f"Error for QY_AB: {error_QY_AB:.5f} ")
    
    print(f"Optimized QY_BA: {QY_BA_opt:.5f}" )
    print(f"Standard deviation of QY_BA in the fit using I0_avg: {std_dev_fit_QY_BA:.5f}")
    print(f"Error for QY_BA: {error_QY_BA:.5f} ")
    
    return QY_AB_opt, QY_BA_opt, error_QY_AB, error_QY_BA

##################################################
def CalculateConcentrations(lambda_meters, 
                            init_conc_A, init_conc_B,
                            timestamps,
                            QY_AB_opt, QY_BA_opt,
                            e_A_inter, e_B_inter,
                            N, V):
    ## Integrate the rate equations with the optimized parameters
    
    conc_opt = odeint(rate_equations, [init_conc_A, init_conc_B], timestamps,
                               args=(lambda_meters, QY_AB_opt, QY_BA_opt, 
                                     e_A_inter, e_B_inter,
                                     N, V))
    
    
    PSS_A = conc_opt[-1,0]/(init_conc_A+init_conc_B)*100
    PSS_B = conc_opt[-1,1]/(init_conc_A+init_conc_B)*100
    
    #print(f"At the PSS, {PSS_A:.2f} % of A and {PSS_B:.2f} % of B")    
    return conc_opt
##################################################
def GetFittedAbs(fit_results, conc_opt,
                 e_A_inter, e_B_inter,
                 timestamps,
                 wavelengths_data):
    ##### Get the fitted absorbance
    result_lmfit=fit_results[0] ## results using I0_avg
    total_abs_fit = conc_opt.dot(np.vstack([e_A_inter, e_B_inter]))
    residuals = result_lmfit.residual.reshape((len(timestamps), 
                                              len(wavelengths_data))).T
    return total_abs_fit, residuals
##################################################
def PlotAndSave(LEDwavelength, datafolder, filename,
                timestamps, conc_opt,
                absorbance, index, 
                total_abs_fit, residuals,
                QY_AB_opt, QY_BA_opt, 
                error_QY_AB, error_QY_BA):
    # ToSave = input("Save plots (Yes or No)? ")
    ToSave = "No" ## for testing
    ##################################################
    
    ######################
    fig = plt.figure(figsize=(8, 4),dpi=600,constrained_layout=True)
    gs = gridspec.GridSpec(4, 2, figure=fig)
    fig.suptitle(f'{LEDwavelength} nm: Integration Method\n \
                  {filename}')
    
    axresults_conc = fig.add_subplot(gs[0:3, 0])
    axresults_Abs = fig.add_subplot(gs[0:3, 1])
    axresults_res = fig.add_subplot(gs[3, 1])
    axresults_notes = fig.add_subplot(gs[3, 0])
    
    axresults_conc.set_title("Concentrations over time")
    
    axresults_conc.plot(timestamps, conc_opt[:,0], label="Species A")
    axresults_conc.plot(timestamps, conc_opt[:,1], label = "Species B")
    
    axresults_conc.set_ylabel("Concentration (M)")
    axresults_conc.yaxis.set_minor_locator(AutoMinorLocator(2))
    axresults_conc.set_xlabel("Time (s)")
    
    axresults_conc.legend()
    
    ##################################################
    ################ NOTES ################
    axresults_notes.axis("off")
    # axresults_notes.text(0, 0, f"Starting Percentage A: {StartPercentage_A} %")
    ##################################################
    #### Plot the experimental data and the fitted total absorbance curve
    axresults_Abs.plot(timestamps, absorbance[index,:], linestyle='-', 
                color="#DCEEFF", linewidth=5, label='Experimental Data')
    axresults_Abs.plot(timestamps, total_abs_fit[:,index], linestyle='--', color="#9ACCFF", 
                label=f"Fitted Total Absorbance Curve\nQY_AB: {QY_AB_opt:.3f}"+u"\u00B1"+f"{error_QY_AB:.3f}\
                \nQY_BA: {QY_BA_opt:.3f}"+u"\u00B1"+f"{error_QY_BA:.3f}")
    axresults_Abs.set_title('Absorbance over time')
    axresults_Abs.set_xlabel('Time (s)')
    axresults_Abs.set_ylabel('Absorbance')
    axresults_Abs.legend()
    
    #### Plot the residuals over time
    axresults_res.plot(timestamps, residuals[index,:], color="#FF952A", label='Residuals')
    axresults_res.set_xlabel('Time (s)')
    axresults_res.set_ylabel('Residual Abs')
    
    #### SAVE PLOTS ####
    if ToSave == "Yes":
        savefilename = datafolder+"\\"+"QY-Integrated_"+str(LEDwavelength)+"nm_"+filename
        plt.savefig(savefilename+".svg",bbox_inches="tight")
        plt.savefig(savefilename+".png",bbox_inches="tight")
        print("Saved plots")
    elif ToSave == "No":
        print("Plots not saved")
    else:
        print("Wrong ToSave input")
    ####################
    plt.show()
#################################################
#%%
