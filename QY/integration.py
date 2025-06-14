"""
Integrate over the emission spectrum of the LED

"""
import numpy as np
from scipy.integrate import trapezoid,odeint
from scipy.signal import savgol_filter
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt

import QY.constants as Constants
import data.experimental_parameters as ExpParams

##################################################

wavelength_low = None
wavelength_high = None
LEDindex_first = None
LEDindex_last = None

##################################################
###################### LED #######################
##################################################
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
    
    return emission_Intensity_proc, first_index, last_index, low, high
    ####################################

##################################################
################### ABSORBANCE ###################
##################################################
def Process_SpectralData(data_pd, wl_low, wl_high, wavelength_of_interest):
    """ Process Spectral Data (the spectra recorded during irradiation) """
    # Clean up the first column by splitting space-separated values and selecting the first one
    def clean_column(value):
        if isinstance(value, str):
            parts = value.split()  # Split the string by spaces
            return float(parts[0])  # Take the first part and convert to float
        return value  # If it's already a number, just return it

    data_pd.iloc[:, 0] = data_pd.iloc[:, 0].apply(clean_column)
    data = data_pd.to_numpy()  # Convert DataFrame to Numpy array

    ## Perform operations on the cleaned numeric data
    low = np.argmin(np.abs(data[:, 0].astype(float) - wl_low))
    high = np.argmin(np.abs(data[:, 0].astype(float) - wl_high))

    wavelengths_lowhigh = data[low:high, 0]
    absorbance_lowhigh = data[low:high, 1:]

    index = np.abs(wavelengths_lowhigh - wavelength_of_interest).argmin()

    return wavelengths_lowhigh, absorbance_lowhigh, index

def Interpolate_Epsilons(IrrSpectra_wavelengths,
                         eA_wavelengths, eA_values,
                         eB_wavelengths, eB_values,
                         LED_wavelengths, LED_intensity_proc):
    ## Interpolate ε_A and ε_B at each wavelength so that all data sets use the same 
    ## wavelengths data
    epsilons_R_interp = np.interp(IrrSpectra_wavelengths, eA_wavelengths, 
                                 eA_values)
    epsilons_P_interp = np.interp(IrrSpectra_wavelengths, eB_wavelengths, 
                                 eB_values)
    emission_interp = np.interp(IrrSpectra_wavelengths, LED_wavelengths,
                                LED_intensity_proc)
    return epsilons_R_interp, epsilons_P_interp, emission_interp

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
########## SOLVE RATE EQUATIONS AND PLOT RESULTS ##########
# StartPercentage_A = float(input("Starting Percentage of A: "))
# StartPercentage_A = float(100)
##################################################

def CreateParameters(absorbance_values, wavelengths_data,
                     e_R_inter, emission_inter):
    """
    Create the parameters needed for the following functions
    """
    #########################################
    ##### trying out: input starting percentages ######
    #########################################
    StartPercentage_R = float(100) ##!!! turn into optional parameter

    # print(f"Integration-CreateParameters===absorbance_values:{absorbance_values}\nand shape:{absorbance_values.shape}")
    # print(f"Integration-CreateParameters===e_A_inter:{e_A_inter}\nand shape:{e_A_inter.shape}")

    initial_conc_R_100 = trapezoid(absorbance_values[:,0],
                                   x=wavelengths_data) / trapezoid(e_R_inter,
                                                                   x=wavelengths_data)
    initial_conc_R = initial_conc_R_100/100*StartPercentage_R
    initial_conc_P_0 = 0
    initial_conc_P = initial_conc_P_0/100*(100-StartPercentage_R)
#########################################
#########################################
    # total_absorbance = absorbance_values.T  # array of total_absorbance values
    lambda_meters = wavelengths_data * 1e-9  ## Convert to meters

    ## Normalize the LED emission spectrum to ensure the area under the curve is 1
    normalized_emission = emission_inter / trapezoid(emission_inter, lambda_meters)

    return initial_conc_R, initial_conc_P, lambda_meters, normalized_emission

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
    k_BA = ExpParams.k_BA
    
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

def objective_function(params, lambda_meters, 
                       init_conc_A, init_conc_B,
                       time, experimental_data, 
                       epsilon_A_lambda, epsilon_B_lambda, 
                       N_lambda, V):
    QY_AB = params['QY_AB'].value
    QY_BA = params['QY_BA'].value

    concentrations_ode = odeint(rate_equations, [init_conc_A, init_conc_B], time,
                               args=(lambda_meters, QY_AB, QY_BA, epsilon_A_lambda, 
                                     epsilon_B_lambda, N_lambda, V))
    
    # Matrix multiplication using dot product
    total_absorbance_ode=concentrations_ode.dot(np.vstack([epsilon_A_lambda,
                                                            epsilon_B_lambda]))

    ##!!! ADD WARNING IF: 
        ## total_absorbance_ode (which is derived from #timestamps)
        ## and experimental_data (which is derived from #spectra)
        ## lengths do not match
    return total_absorbance_ode - experimental_data.T

def MinimizeQYs(I0_list,
                emission_norm, lambda_meters, 
                init_conc_A, init_conc_B,
                timestamps, absorbance_values,
                e_A_inter, e_B_inter,
                V):
    
    h = Constants.h
    c = Constants.c
    Avogadro = Constants.Avogadro
    
    ### loop over the powers in the list 
    fit_results=[]
    for total_power in I0_list:
        # power_at_each_wavelength = emission_norm * total_power * 1e-6 ## Photon flux in mmol/s per wavelength # for microwatt
        power_at_each_wavelength = emission_norm * total_power * 1e-3 ## Photon flux in mmol/s per wavelength # for milliwatt
        N = power_at_each_wavelength / (h * c / lambda_meters) / Avogadro * 1000
        
        ## Add wavelength-specific parameters to the lmfit Parameters object
        params = Parameters()
        params.add('QY_AB', value=0.5, min=0, max=1) ## with boundaries
        params.add('QY_BA', value=0.5, min=0, max=1) ## with boundaries
        
        ## Minimize the objective function using lmfit
        ## Fitting method: Levenberg-Marquardt (default)
        result_lmfit = minimize(objective_function, 
                                params, args=(lambda_meters,
                                              init_conc_A, init_conc_B,
                                              timestamps, absorbance_values, 
                                              e_A_inter, e_B_inter, 
                                              N, V))

        fit_results.append(result_lmfit)
    return N, fit_results
