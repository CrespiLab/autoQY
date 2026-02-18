"""
Integrate over the emission spectrum of the LED

"""
import numpy as np
from scipy.integrate import trapezoid,odeint
from scipy.signal import savgol_filter, find_peaks, peak_widths
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt

import QY.constants as Constants
import data.experimental_parameters as ExpParams
import data.calc_settings as CalcSettings

##################################################

wavelength_low = None
wavelength_high = None
LEDindex_first = None
LEDindex_last = None

##################################################
###################### LED #######################
##################################################
def Process_LEDemission(wavelengths_LED, intensity_LED):
    ''' Apply smoothing and remove negative values in LED emission data '''
    message_blcorr = None
    message = None
    emission_Intensity_proc = None
    
    try:
        emission_Intensity_smoothed = savgol_filter(intensity_LED, 12,3) ## Smoothing
        ##!!! REMOVE SMOOTHING (AND REPLACE BY PEKARIAN FIT)

        if CalcSettings.BaselineCorrection_LED == "ON":
            message_blcorr, emission_Intensity_proc = BaselineCorrection_LED(wavelengths_LED,emission_Intensity_smoothed)
        elif CalcSettings.BaselineCorrection_LED == "OFF":
            emission_Intensity_proc = emission_Intensity_smoothed
            message_blcorr = "No baseline correction of LED emission spectrum."
        emission_Intensity_proc[emission_Intensity_proc[:]<0] = 0 ## removal of negative values
    except Exception as e:
        message = f"Failed to process LED emission spectrum: {e}."

    return message_blcorr, message, emission_Intensity_proc

def BaselineCorrection_LED(wavelengths, intensities):
    try:
        # Step 1: Find peaks
        peaks, _ = find_peaks(intensities, height=np.max(intensities)*0.5)
        main_peak_index = peaks[np.argmax(intensities[peaks])]
        main_peak_x = wavelengths[main_peak_index]
        
        # Step 2: Estimate FWHM of the main peak
        results_half = peak_widths(intensities, [main_peak_index], rel_height=0.5)
        fwhm = results_half[0][0] * (wavelengths[1] - wavelengths[0])  # convert index width to nm
        
        # Step 3: Define adaptive window size
        window_width = 10 * fwhm ##!!! make this an option in the GUI?
        mask = (wavelengths < main_peak_x - window_width) | (wavelengths > main_peak_x + window_width)
        
        baseline_x = wavelengths[mask]
        baseline_y = intensities[mask]
        
        # Step 4: Fit linear baseline
        coeffs = np.polyfit(baseline_x, baseline_y, deg=0) ## vertical offset
        baseline_fit = np.polyval(coeffs, wavelengths)
        
        # Step 5: Apply correction
        intensities_baselinecorrected = intensities - baseline_fit

        message = f"Performed baseline correction of LED Emission spectrum. Exclusion window: ±{window_width:.2f} nm around the peak at {main_peak_x} nm."
    except Exception as e:
        intensities_baselinecorrected = None
        message = f"Failed to perform baseline correction: {e}"
    
    return message, intensities_baselinecorrected

def LEDemission_WavelengthLimits(wavelengths_LED, intensity_LED_proc, threshold):
    ''' Obtain wavelength limits of LED emission spectrum according to a set intensity threshold'''
    #### Cutting spectrum to match emission of LED 
    above_threshold_indices = np.where(intensity_LED_proc > threshold)[0] 
    
    #### Find the first and last indices
    first_index = above_threshold_indices[0]
    last_index = above_threshold_indices[-1]
    
    low = wavelengths_LED[first_index]
    high = wavelengths_LED[last_index]
    
    return low, high
    ####################################

def Epsilons_WavelengthLimits(epsilon_A_wavelengths, epsilon_B_wavelengths):
    ''' Obtain limits of epsilons wavelength data '''
    
    
    epsilon_wl_low = min([epsilon_A_wavelengths[0],epsilon_B_wavelengths[0]])
    epsilon_wl_high = min([epsilon_A_wavelengths[-1],epsilon_B_wavelengths[-1]])
    
    if epsilon_wl_low < CalcSettings.wl_low:
        wavelength_low = CalcSettings.wl_low
    else:
        wavelength_low = epsilon_wl_low
    
    if epsilon_wl_high > CalcSettings.wl_high:
        wavelength_high = CalcSettings.wl_high
    else:
        wavelength_high = epsilon_wl_high
    
    return wavelength_low, wavelength_high

##################################################
################### ABSORBANCE ###################
##################################################
def Process_SpectralData(data_pd, wl_low, wl_high, wavelength_of_interest):
    """ Process Spectral Data (the spectra recorded during irradiation) """
    # Clean up the first column by splitting space-separated values and selecting the first one
    ##!!! Why is this done? Is it a remnant of a previous method of importing?
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

def Interpolate_Spectra(reference_wavelengths,
                         target_wavelengths, target_values):
    '''
    Interpolate all spectral data according to wavelengths data of Irradiation Spectra,
        so that all datasets use the same wavelengths data
    LED emission data also at the same time gets cut to the same wavelength range as all other datasets.

    Parameters
    ----------
    reference_wavelengths : TYPE
        DESCRIPTION.
    target_wavelengths : TYPE
        DESCRIPTION.
    target_values : TYPE
        DESCRIPTION.

    Returns
    -------
    target_interp : numpy array
        interpolated y-data according to reference wavelengths.
    '''
    target_interp = np.interp(reference_wavelengths,
                             target_wavelengths,
                             target_values)

    return target_interp

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
############## SOLVE RATE EQUATIONS ##############
##################################################

def CreateParameters(absorbance_values, wavelengths_data,
                     e_R_inter, emission_inter):
    """
    Create the parameters needed for the following functions
    """
    #########################################
    ##### trying out: input starting percentages ######
    #########################################
    StartPercentage_R = float(100) 
    initial_conc_R_100 = trapezoid(absorbance_values[:,0],
                                   x=wavelengths_data) / trapezoid(e_R_inter,
                                                                   x=wavelengths_data)
    initial_conc_R = initial_conc_R_100/100*StartPercentage_R
    initial_conc_P_0 = 0
    initial_conc_P = initial_conc_P_0/100*(100-StartPercentage_R)
#########################################
#########################################
    lambda_meters = wavelengths_data * 1e-9  ## Convert to meters

    ## Normalize the LED emission spectrum to ensure the area under the curve is 1
    normalized_emission = emission_inter / trapezoid(emission_inter, lambda_meters)

    return initial_conc_R, initial_conc_P, lambda_meters, normalized_emission

def CreateParameters_Conc(absorbance_values, wavelengths_data,
                          fractions_R, fractions_P,
                          e_R_inter, e_P_inter,
                          emission_inter):
    """
    Create the parameters needed for the subsequent functions used by the 
    ODEMethod = Concentrations
    """
    #########################################
    total_conc = trapezoid(absorbance_values[:,0], x=wavelengths_data)\
        / (fractions_R[0]*trapezoid(e_R_inter, x=wavelengths_data)\
           + fractions_P[0]*trapezoid(e_P_inter, x=wavelengths_data))
    
    concs_RP = (np.stack((fractions_R,fractions_P),axis=1))*total_conc ## make array of concentrations of R and P
    initial_conc_R, initial_conc_P = concs_RP[0]
    #########################################
    lambda_meters = wavelengths_data * 1e-9  ## Convert to meters

    ## Normalize the LED emission spectrum to ensure the area under the curve is 1
    normalized_emission = emission_inter / trapezoid(emission_inter, lambda_meters)

    return total_conc, initial_conc_R, initial_conc_P, concs_RP, lambda_meters, normalized_emission

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
    
    ### ----- Absorbance contributions -----
    ### Total absorbance spectrum (base-10)
    total_abs = (A * epsilon_A_lambda + B * epsilon_B_lambda) 

    # Fraction of light absorbed by A vs. B at each λ
    frac_A = (A * epsilon_A_lambda) / total_abs
    frac_B = (B * epsilon_B_lambda) / total_abs

    # Photon absorption rates (photons·s⁻¹·L⁻¹)
    R_A_abs = np.trapezoid(-(frac_A * (1 - 10**(-total_abs)) * N_lambda) / V, x=lambda_meters)
    R_B_abs = np.trapezoid((frac_B * (1 - 10**(-total_abs)) * N_lambda) / V, x=lambda_meters)

    # ----- Kinetics -----
    dAdt = QY_AB * R_A_abs + QY_BA * R_B_abs + k_BA * B
    dBdt = -dAdt

    return [dAdt, dBdt]

#########################################

def objective_function(params, 
                       lambda_meters, 
                       init_conc_A, init_conc_B,
                       time, experimental_Abs, 
                       epsilon_A_lambda, epsilon_B_lambda, 
                       N_lambda, V):
    '''
    Function that minimises the QYs according to the absorbance values.
    '''
    QY_AB = params['QY_AB'].value
    QY_BA = params['QY_BA'].value

    concentrations_ode = odeint(rate_equations, 
                                [init_conc_A, init_conc_B], 
                                time,
                               args=(lambda_meters, 
                                     QY_AB, QY_BA, 
                                     epsilon_A_lambda, epsilon_B_lambda, 
                                     N_lambda, V))
    
    # Matrix multiplication using dot product
    total_absorbance_ode=concentrations_ode.dot(np.vstack([epsilon_A_lambda,
                                                            epsilon_B_lambda]))

    return total_absorbance_ode - experimental_Abs.T

def objective_function_conc(params,
                            lambda_meters,
                            init_conc_A, init_conc_B,
                            time,  experimental_concs,       # ← NEW
                            epsilon_A_lambda, epsilon_B_lambda,
                            N_lambda, V):
    '''
    Function that minimises the QYs according to the concentrations.
    '''
    QY_AB = params['QY_AB'].value
    QY_BA = params['QY_BA'].value

    concs_ode = odeint(rate_equations,
                       [init_conc_A, init_conc_B],
                       time,
                       args=(lambda_meters, 
                             QY_AB, QY_BA,
                             epsilon_A_lambda, epsilon_B_lambda,
                             N_lambda, V))

    return (concs_ode - experimental_concs).ravel()   # flattened residual

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

def MinimizeQYs_Conc(I0_list,
                emission_norm, 
                lambda_meters, 
                init_conc_R, init_conc_P,
                timestamps, concs_RP,
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
        result_lmfit = minimize(objective_function_conc, 
                                params, 
                                args=(lambda_meters,
                                      init_conc_R, init_conc_P,
                                      timestamps, concs_RP,  ## concentrations instead of Abs
                                      e_A_inter, e_B_inter, 
                                      N, V))

        fit_results.append(result_lmfit)
    return N, fit_results