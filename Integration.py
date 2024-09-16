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

[] include option to set starting percentage (see below)

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

#%% DATA AND PARAMETERS
##!!! FOR AUTOQUANT.PY

script_dir = os.path.dirname(os.path.abspath(__file__)) ## Directory where the script is located
datafolder = script_dir+r"\ExampleData"
###################################################
############# FILE FORMAT #############
# FileFormat = "Spectragryph"
FileFormat = "Not"
##############################

############# EPSILONS #############
if FileFormat == "Spectragryph":
    epsilonsfile_A = datafolder+"\\"+"trans_azobenzene_epsilons.dat"
    epsilonsfile_B = datafolder+"\\"+"cis_azobenzene_epsilons.dat"
elif FileFormat == "Not":
    epsilonsfile_A = datafolder+"\\"+"trans_azobenzene_epsilons.csv"
    epsilonsfile_B = datafolder+"\\"+"cis_azobenzene_epsilons.csv"
###################################################

############################################
############# IRRADIATION DATA #############
datafilename = "Azobenzene_340nm_processed"

data_file = datafolder+"\\"+datafilename+".dat" 
log_file = datafolder+"\\"+"Azobenzene_340nm_log.csv"
###################################################

########################################
############# LED EMISSION #############
#### Input LED emission data in nanometers: just the raw data
##=> it gets normalised during fitting!


LED_file = datafolder+"\\"+"LED340nm_normalisedarea.csv"
NominalWavelengthLED = 340    # Wavelength of interest
###################################################

#######################
######## POWER ########
I0_avg = 743             # Photon flux in microWatt
I0_err = 4                # Error on photon fplux in microWatt
I0_list=[I0_avg, I0_avg+I0_err, I0_avg-I0_err]
###################################################

############# EXPERIMENTAL PARAMETERS #############
V = 3.0                 # Volume in ml
k_BA = 7.240e-7         # Thermal back reaction rate s-1
###################################################

###################################################
################## CONSTANTS ######################
h = 6.626e-34           # Planck's constant in J·s
c = 299792458           # Speed of light in m/s
Avogadro = 6.022e23     # Avogadro's number
###################################################

#%% CREATE DATASETS

##################################################
################### TIMESTAMPS ###################
##################################################
##!!! FOR AUTOQUANT.PY

#### Import log file of irradiation data
def GetTimestamps(LogFile):
    log = pd.read_csv(LogFile,
                      sep = ",", decimal = ".", skiprows = 1, header=None,)
    log_t=log[log[3] == 'Measure']
    t=log_t[2]
    timestamps=t.to_numpy()
    return timestamps
timestamps = GetTimestamps(log_file)

#%% LED 

############# ANALYSIS PARAMETERS #################
threshold_LED = 500     # Threshold for part of spectrum where there is LED emission 
##!!! make this a user-defined option

##################################################
###################### LED #######################
##################################################
def Import_LEDemission(file_LEDemission_raw):
    ## Import (raw) emission LED_file
    if FileFormat == "Spectragryph":
        emission_data = pd.read_csv(file_LEDemission_raw, sep = '\t', usecols = lambda x: x not in ["Wavenumbers [1/cm]"]) # in nanometers
    elif FileFormat == "Not":
        emission_data = pd.read_csv(file_LEDemission_raw) # in meters
    
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
    print(f"wavelength_low: {low}\nwavelength_high: {high}")
    
    return emission_Intensity_proc, first_index, last_index, low, high
    ####################################

##################################################
################### ABSORBANCE ###################
##################################################

def Import_SpectralData(file, wl_low, wl_high, wavelength_of_interest):
    ## import Absorbance data_file (.dat from Spectragryph, so with Wavenumbers column)
    if FileFormat == "Spectragryph":
        data_pd = pd.read_csv(file, sep = '\t', usecols = lambda x: x not in ["Wavenumbers [1/cm]"]) # in nanometers
    elif FileFormat == "Not":
        data_pd = pd.read_csv(file, sep = ' ', header=None) # in nanometers

    data = data_pd.to_numpy() ## convert DataFrame to Numpy array

    low = np.argmin(np.abs(data[:,0] - wl_low))
    high = np.argmin(np.abs(data[:,0] - wl_high))
    
    wavelengths_lowhigh = data[low:high, 0]
    absorbance_lowhigh = data[low:high, 1:]
    index = np.abs(wavelengths_lowhigh - wavelength_of_interest).argmin()

    return wavelengths_lowhigh, absorbance_lowhigh, index
    ####################################

##################################################
###################### PLOT ######################
##################################################
##!!! FOR AUTOQUANT.PY

def Plot_LEDemission_Processed():
    fig, axdata = plt.subplots(2,1,figsize=(6,6),
                                dpi=600,constrained_layout=True)
    
    for i in range(0,SpectralData_Abs.shape[1]):
        axdata[0].plot(SpectralData_Wavelengths, SpectralData_Abs.T[i])
    
    for i in axdata:
        i.set_xlabel("Wavelength (nm)")
        i.set_xlim(220,650)
    
    axdata[0].set_title("Spectral Data")
    axdata[0].set_ylabel("Absorbance")
    axdata[0].set_ylim(-0.05,2.5)
    
    axdata[1].plot(emission_wavelengths,emission_Intensity,
                    label="Untreated (in this .py file at least)")
    axdata[1].plot(emission_wavelengths[LEDindex_first:LEDindex_last],
                    emission_Intensity_proc[LEDindex_first:LEDindex_last],
                    label="Smoothed, removed zeroes\nand applied threshold")
    axdata[1].legend(fontsize=8)
    axdata[1].set_title("LED emission")
    axdata[1].set_ylabel("Intensity")
    plt.show()
    ####################################

####################################
emission_wavelengths, emission_Intensity = Import_LEDemission(LED_file)
emission_Intensity_proc, LEDindex_first, LEDindex_last, wavelength_low, wavelength_high = Processing_LEDemission(emission_wavelengths, emission_Intensity, threshold_LED)
SpectralData_Wavelengths, SpectralData_Abs, SpectralData_Index = Import_SpectralData(data_file, wavelength_low, wavelength_high, NominalWavelengthLED)
Plot_LEDemission_Processed()
####################################
#%% EPSILONS LOAD
""" 
Epsilons datafiles are in a certain format
TO DO for GUI:
    - select delimiter: ',' or '\t'
    - select; ignore Wavenumbers column or not
    
"""

def Import_Epsilons(A, B):
    if FileFormat == "Spectragryph":
        epsilon_A_data = pd.read_csv(A, delimiter='\t',
                                     usecols = lambda x: x not in ["Wavenumbers [1/cm]"])
        epsilon_B_data= pd.read_csv(B, delimiter='\t', 
                                    usecols = lambda x: x not in ["Wavenumbers [1/cm]"])
    elif FileFormat == "Not":
        epsilon_A_data = pd.read_csv(A, delimiter=',')
        epsilon_B_data = pd.read_csv(B, delimiter=',')
        
    epsilon_A_data.columns = ['Wavelengths', 'Epsilons'] ## rename columns
    epsilon_A_wavelengths = epsilon_A_data['Wavelengths'].values
    epsilon_A_values = epsilon_A_data['Epsilons'].values
    
    epsilon_B_data.columns = ['Wavelengths', 'Epsilons'] ## rename columns
    epsilon_B_wavelengths = epsilon_B_data['Wavelengths'].values
    epsilon_B_values = epsilon_B_data['Epsilons'].values
    return epsilon_A_wavelengths, epsilon_A_values, epsilon_B_wavelengths, epsilon_B_values
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
##!!! FOR AUTOQUANT.PY

def Plot_Epsilons():
    fig, axdata_interp = plt.subplots(2,1,figsize=(6,6),
                                dpi=600,constrained_layout=True)
    
    axdata_interp[0].plot(SpectralData_Wavelengths, epsilon_A_interp, label="A")
    axdata_interp[0].plot(SpectralData_Wavelengths, epsilon_B_interp, label="B")
    axdata_interp[0].legend(fontsize=12)
    
    for i in axdata_interp:
        i.set_xlabel("Wavelength (nm)")
        i.set_xlim(220,650)
    
    axdata_interp[0].set_title("Epsilons, interpolated")
    axdata_interp[0].set_ylabel(r"$\epsilon$ (M$^{-1}$ cm$^{-1}$)")
    # axdata_interp[0].set_ylim(-0.05,2.5)
    
    axdata_interp[1].set_title("LED emission, interpolated")
    axdata_interp[1].set_ylabel("Intensity")
    axdata_interp[1].plot(SpectralData_Wavelengths,emission_interp)
    # axdata_interp[1].legend(fontsize=8)
    
    plt.show()

##################################################

epsilon_A_wavelengths, epsilon_A_values, epsilon_B_wavelengths, epsilon_B_values = Import_Epsilons(epsilonsfile_A, epsilonsfile_B)

epsilon_A_interp, epsilon_B_interp, emission_interp = Interpolate_Epsilons(SpectralData_Wavelengths,
                     epsilon_A_wavelengths, epsilon_A_values,
                     epsilon_B_wavelengths, epsilon_B_values,
                     emission_wavelengths, emission_Intensity_proc)

Plot_Epsilons()

##################################################
#%% SOLVE RATE EQUATIONS AND PLOT RESULTS

"""
TO DO:
[] add StartPercentage_A as an optional variable
    
"""
# StartPercentage_A = float(input("Starting Percentage of A: "))
StartPercentage_A = float(100)

def rate_equations(concentrations, time, QY_AB, QY_BA, epsilon_A_lambda, 
                   epsilon_B_lambda, N_lambda, V):
    A, B = concentrations
    total_absorbance = A * epsilon_A_lambda + B * epsilon_B_lambda

    dAdt = QY_AB * trapezoid(-((A * epsilon_A_lambda / total_absorbance) \
            * (1 - 10 ** (-total_absorbance)) * N_lambda) / V,lambda_meters)+ \
            QY_BA * trapezoid(((B * epsilon_B_lambda / total_absorbance) \
            * (1 - 10 ** (-total_absorbance)) * N_lambda) / V,lambda_meters) \
            + k_BA * B

    dBdt = -dAdt

    return [dAdt, dBdt]

def objective_function(params, time, experimental_data, epsilon_A_lambda, 
                       epsilon_B_lambda, N_lambda, V):
    QY_AB = params['QY_AB'].value
    QY_BA = params['QY_BA'].value

    concentrations_ode = odeint(rate_equations, [initial_conc_A, 
                                initial_conc_B], time,
                               args=(QY_AB, QY_BA, epsilon_A_lambda, 
                                     epsilon_B_lambda, N_lambda, V))
    
    # Matrix multiplication using dot product
    total_absorbance_ode=concentrations_ode.dot(np.vstack([epsilon_A_lambda,
                                                            epsilon_B_lambda]))

    return total_absorbance_ode - experimental_data.T

# loop over the powers in the list 
fit_results=[]
for total_power in I0_list:
    
    lambda_meters = wavelengths_data * 1e-9  # Convert to meters
    
    ## Normalize the emission spectrum to ensure the area under the curve is 1
    normalized_emission=emission_interp / trapezoid(emission_interp, lambda_meters)
    
    power_at_each_wavelength= normalized_emission * total_power * 1e-6 ## Photon flux in mmol/s per wavelength
    N = power_at_each_wavelength / (h * c / lambda_meters) / Avogadro * 1000  
    
    #########################################
    ##### trying out: input starting percentages ######
    #########################################
    
    initial_conc_A_100 = trapezoid(absorbance_values[:,0],x=wavelengths_data)\
        / trapezoid(epsilon_A_interp,x=wavelengths_data)
    initial_conc_A = initial_conc_A_100/100*StartPercentage_A
    
    initial_conc_B_0 = 0
    initial_conc_B = initial_conc_B_0/100*(100-StartPercentage_A)
    
    #########################################
    #########################################
    
    epsilon_A = epsilon_A_interp            # array of epsilon_A values
    epsilon_B = epsilon_B_interp            # array of epsilon_B values
    total_absorbance = absorbance_values.T  # array of total_absorbance values
    N = N                                   # array of N values
    time = timestamps
    
    ## Add wavelength-specific parameters to the lmfit Parameters object
    params = Parameters()
    params.add('QY_AB', value=0.5, min=0, max=1) ## with boundaries
    params.add('QY_BA', value=0.5, min=0, max=1) ## with boundaries
    
    ## Minimize the objective function using lmfit   
    result_lmfit = minimize(objective_function, 
                            params, args=(timestamps, absorbance_values, 
                                          epsilon_A_interp, epsilon_B_interp, 
                                          N, V))
    print(f"Results for I0={total_power} mW: ")
    print(result_lmfit.params)
    print("\n")
    fit_results.append(result_lmfit)

##!!! return fit_results from the function

#%% FOR AUTOQUANT.PY

##################################################
################# EXTRACT RESULTS#################
##################################################
## Extract the optimized parameters and their standard deviations
result_lmfit=fit_results[0]
## Extract the optimized parameters and their standard deviations
QY_AB_opt = result_lmfit.params['QY_AB'].value
std_dev_QY_AB = result_lmfit.params['QY_AB'].stderr

QY_BA_opt = result_lmfit.params['QY_BA'].value
std_dev_QY_BA = result_lmfit.params['QY_BA'].stderr


QY_AB_opt_min=fit_results[1].params['QY_AB'].value \
                    - fit_results[1].params['QY_AB'].stderr
QY_AB_opt_max=fit_results[2].params['QY_AB'].value \
                    - fit_results[2].params['QY_AB'].stderr

QY_BA_opt_min=fit_results[1].params['QY_BA'].value \
                    - fit_results[1].params['QY_BA'].stderr
QY_BA_opt_max=fit_results[2].params['QY_BA'].value \
                    - fit_results[2].params['QY_BA'].stderr
## Calculate error for QY_AB and QY_BA based on the error in the power 
error_QY_AB = max([QY_AB_opt-QY_AB_opt_min, QY_AB_opt_max-QY_AB_opt])
error_QY_BA = max([QY_BA_opt-QY_BA_opt_min, QY_BA_opt_max-QY_BA_opt])

## Print the results
print(f"Optimized QY_AB: {QY_AB_opt:.5f}" )
print(f"Error for QY_AB: {error_QY_AB:.5f} ")

print(f"Optimized QY_BA: {QY_BA_opt:.5f}" )
print(f"Error for QY_BA: {error_QY_BA:.5f} ")
##################################################
################## PLOT AND SAVE #################
##################################################
ToSave = input("Save plots (Yes or No)? ")
##################################################

## Integrate the rate equations with the optimized parameters
conc_opt= odeint(rate_equations, [initial_conc_A, initial_conc_B], time,
                  args=(QY_AB_opt, QY_BA_opt, epsilon_A_interp, 
                                        epsilon_B_interp, N, V))
PSS_A=conc_opt[-1,0]/(initial_conc_A+initial_conc_B)*100
PSS_B=conc_opt[-1,1]/(initial_conc_A+initial_conc_B)*100

print(f"At the PSS, {PSS_A:.2f} % of A and {PSS_B:.2f} % of B")
######################
fig = plt.figure(figsize=(8, 4),dpi=600,constrained_layout=True)
gs = gridspec.GridSpec(4, 2, figure=fig)
fig.suptitle(f'{wavelength_LED} nm: Integration Method\n \
              {datafilename}')

axresults_conc = fig.add_subplot(gs[0:3, 0])
axresults_Abs = fig.add_subplot(gs[0:3, 1])
axresults_res = fig.add_subplot(gs[3, 1])
axresults_notes = fig.add_subplot(gs[3, 0])

axresults_conc.set_title("Concentrations over time")

axresults_conc.plot(time, conc_opt[:,0], label="Species A")
axresults_conc.plot(time, conc_opt[:,1], label = "Species B")

axresults_conc.set_ylabel("Concentration (M)")
axresults_conc.yaxis.set_minor_locator(AutoMinorLocator(2))
axresults_conc.set_xlabel("Time (s)")

axresults_conc.legend()

##################################################
################ NOTES ################
axresults_notes.axis("off")
axresults_notes.text(0, 0, f"Starting Percentage A: {StartPercentage_A} %")
##################################################

##################################################
##### Get the fitted absorbance
total_abs_fit = conc_opt.dot(np.vstack([epsilon_A_interp, epsilon_B_interp]))
residuals=result_lmfit.residual.reshape((len(timestamps), 
                                          len(wavelengths_data))).T

#### Plot the experimental data and the fitted total absorbance curve
axresults_Abs.plot(timestamps, absorbance_values[index,:], linestyle='-', 
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
    savefilename = datafolder+"\\"+"QY-Integrated_"+str(wavelength_LED)+"nm_"+datafilename
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
