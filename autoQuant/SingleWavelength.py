# -*- coding: utf-8 -*-
"""
@author: Jorn
Date started: Sep 16th 2024
Based on: script by Anouk Volker from our publication

Calculate the QYs from the absorption data at a single wavelength. 

##!!!
TO DO:
[] Create functions for the necessary code
[] Remove unnecessary code

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import odeint
from lmfit import Model, Parameters

import autoQuant.ExpParams as ExpParams
import autoQuant.Constants as Constants

# #%% Defining main variables
# ############# DATA #############
# #### Replace with filepath and filenames of your own data ####
# log_file= r"log.csv"
# #################################
# #### 365 nm ####
# ##!!! APPLY DIFFERENT METHOD USING ISOSBESTIC POINT??

# # datafolder =r"C:\Users\jorst136\Documents\Postdoc\Projects\GeorgeWilliams_ITI\Experiments\365nm"
# # datafilename = "GW1_365nm_processed"

# # data_file = datafolder+"\\"+datafilename+".dat" 

# # wavelength_nm = 365 ## Wavelength of irradiation

# # ### epsilons ####
# # epsilon_A = 5791.02136 ## Molar absorptivity of A at wavelength in M-1 cm-1 ## GW
# # epsilon_B = 6062.44385 ## Molar absorptivity of B at wavelength in M-1 cm-1 ## GW

# # ### POWER ####
# # I0_avg = 2160.282943 ## Photon flux in microWatt
# # I0_err = 1.187980054 ## Error on photon fplux in microWatt

# #################################
# #### 395 nm ####
# # datafolder =r"C:\Users\jorst136\Documents\Postdoc\Projects\GeorgeWilliams_hydrazoneswitch\Experiments\395nm"
# # datafilename = "GW1_395nm_processed"

# # data_file = datafolder+"\\"+datafilename+".dat" 

# # wavelength_nm = 395 ## Wavelength of irradiation

# # #### epsilons ####
# # epsilon_A = 8563 ## Molar absorptivity of A at wavelength in M-1 cm-1 ## Jacob's
# # epsilon_B = 4467 ## Molar absorptivity of B at wavelength in M-1 cm-1 ## Jacob's

# # epsilon_A = 8606.128 ## Molar absorptivity of A at wavelength in M-1 cm-1 ## GW
# # epsilon_B = 5025.226 ## Molar absorptivity of B at wavelength in M-1 cm-1 ## GW

# # ### POWER ####
# # I0_avg = 2955.958474 ## Photon flux in microWatt
# # I0_err = 1.53969162 ## Error on photon fplux in microWatt

# #################################
# #### 455 nm ####
# datafolder =r"C:\Users\jorst136\Documents\Postdoc\Projects\GeorgeWilliams_hydrazoneswitch\Experiments\455nm"
# datafilename= "GW1_455nm_processed"

# data_file = datafolder+"\\"+datafilename+".dat" 

# wavelength_nm = 455 ## Wavelength of irradiation

# #### epsilons ####
# # epsilon_A = 412 ## Molar absorptivity of A at wavelength in M-1 cm-1 ## Jacob's
# # epsilon_B = 21.5 ## Molar absorptivity of B at wavelength in M-1 cm-1 

# epsilon_A = 445.5937 ## Molar absorptivity of A at wavelength in M-1 cm-1 ## GW
# epsilon_B = 37.92818 ## Molar absorptivity of B at wavelength in M-1 cm-1 ## GW

# #### POWER ####
# I0_avg = 4292.117467 ## Photon flux in microWatt
# I0_err = 2.432442109 ## Error on photon fplux in microWatt
# #################################
# ##################################################################
# ##################################################################


#%% Loading experimental data

# data = np.loadtxt(os.path.join(datafolder, data_file),skiprows=1,usecols=lambda x: x not in 1)

# data_df = pd.read_csv(os.path.join(datafolder, data_file),sep = '\t', usecols=lambda x: x not in ["Wavenumbers [1/cm]"])
# data = data_df.to_numpy()

# log = pd.read_csv(os.path.join(datafolder, log_file),sep = ",", decimal = ".", 
#                   skiprows = 1, header=None,)

# # Get time from the log file
# log_t=log[log[3] == 'Measure']
# t=log_t[2]
# timestamps=t.to_numpy()

# Get the absorbance data closest to the wavelength specified
def find_closest_index(array_2d, target_value):
    # Extract the first column of the 2D array
    first_column = array_2d[:, 0]

    # Find the index of the closest value to the target value
    closest_index = np.abs(first_column - target_value).argmin()

    return closest_index


def CreateParameters(absorbance_data,wavelengths_data,
                     epsilon_R):

    index=find_closest_index(absorbance_data, wavelengths_data)
    absorbance_values = absorbance_data[index, 1:]
    
    ## Initial concentration of A (based on the first experimental value)
    initial_conc_R = absorbance_values[0] / epsilon_R
    initial_conc_P = 0.0

    total_absorbance = absorbance_values.T  # array of total_absorbance values
    lambda_meters = wavelengths_data * 1e-9  ## Convert to meters

    return initial_conc_R, initial_conc_P, total_absorbance, lambda_meters



#%% CALCULATE QYs, PLOT AND SAVE
def rate_equations(concentrations, time, 
                   QY_AB, QY_BA,
                   epsilon_R, epsilon_P,
                   N, V):
    """
    Rate equations for the isomerization, as described in 10.1038/srep41145 
    (Stranius & Börjesson)

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

    Returns
    -------
    [dAdt, dBdt] : list
        List with the functions for the changes in concentration over time

    """
    k_BA = ExpParams.k_BA
    A, B = concentrations
    total_absorbance = A * epsilon_R + B * epsilon_P

    dAdt = (N/V) * ((1 - 10 ** (-total_absorbance ))/total_absorbance) \
    * (QY_BA * B * epsilon_P - QY_AB * A * epsilon_R) + k_BA * B

    dBdt = -dAdt

    return [dAdt, dBdt]


def absorbance_func(time, QY_AB, QY_BA, initial_concentrations): 
    """
    Calculating concentrations and then converting to absorbance to match 
    the output of the experimental data

    Parameters
    ----------
    initial_concentrations : tuple
        Concentrations of A and B
    time : numpy array
        Time 
    QY_AB : float
        Quantum Yield of A --> B
    QY_BA : float
        Quantum Yield of B --> A


    Returns
    -------
    total_absorbance_ode : function
        Ordinary Differential Equation to solve

    """
    initial_concentration_A, initial_concentration_B = initial_concentrations
    concentrations_ode = odeint(rate_equations, initial_concentrations, time, 
                                args = (QY_AB, QY_BA) )
    total_absorbance_ode = concentrations_ode[:, 0] * epsilon_A \
                            + concentrations_ode[:, 1] * epsilon_B
    return total_absorbance_ode

    
I0_list=[I0_avg, I0_avg+I0_err, I0_avg-I0_err]
fit_results=[]

# Solve for the QYs for each of photon fluxes
for I0 in I0_list:
    I0_watt = I0 * 1e-6                              # Convert to watts
    lambda_meters = wavelength_nm * 1e-9             # Convert to meters
    N = I0_watt/(h*c/lambda_meters)/Avogadro * 1000  # Photon flux in mmol/s
    
    # Creating a model from the absorbance function
    absorbance_model=Model(absorbance_func, 
                           independent_vars=['time', 'initial_concentrations'])
    
    params = Parameters()
    params.add('QY_AB', value=0.5, min=0, max=1) ## with boundaries
    params.add('QY_BA', value=0.5, min=0, max=1) ## with boundaries
    
    # Fitting the data with the parameters
    result_lmfit=absorbance_model.fit(absorbance_values, params, 
                                        time=timestamps, 
                                        initial_concentrations=initial_conc)
    print(f"Results for I0={I0} mW: ")
    print(result_lmfit.fit_report())
    print("\n")
    fit_results.append(result_lmfit)

#####################################################
############ RETRIEVE PARAMETERS ####################
#####################################################
# result_lmfit=fit_results[0]
# ## Extract the optimized parameters and their standard deviations
# QY_AB_opt = result_lmfit.params['QY_AB'].value
# std_dev_QY_AB = result_lmfit.params['QY_AB'].stderr

# QY_BA_opt = result_lmfit.params['QY_BA'].value
# std_dev_QY_BA = result_lmfit.params['QY_BA'].stderr


# QY_AB_opt_min=fit_results[1].params['QY_AB'].value \
#                 - fit_results[1].params['QY_AB'].stderr
# QY_AB_opt_max=fit_results[2].params['QY_AB'].value \
#                 - fit_results[2].params['QY_AB'].stderr

# QY_BA_opt_min=fit_results[1].params['QY_BA'].value \
#                 - fit_results[1].params['QY_BA'].stderr
# QY_BA_opt_max=fit_results[2].params['QY_BA'].value \
#                 - fit_results[2].params['QY_BA'].stderr
# # Calculate error for QY_AB and QY_BA based on the error in the power 
# error_QY_AB = max([QY_AB_opt-QY_AB_opt_min, QY_AB_opt_max-QY_AB_opt])
# error_QY_BA = max([QY_BA_opt-QY_BA_opt_min, QY_BA_opt_max-QY_BA_opt])

# # Print the results
# print(f"Optimized QY_AB: {QY_AB_opt:.5f}" )
# print(f"Error for QY_AB: {error_QY_AB:.5f} ")


# print(f"Optimized QY_BA: {QY_BA_opt:.5f}" )
# print(f"Error for QY_BA: {error_QY_BA:.5f} ")


# conc_opt = odeint(rate_equations, [initial_conc_A, initial_conc_B], 
#                   timestamps, args=(QY_AB_opt, QY_BA_opt))

# absorbance_fitted = absorbance_func(timestamps, QY_AB_opt, QY_BA_opt, 
#                                     initial_conc)
# PSS_A=conc_opt[-1,0]/(initial_conc_A+initial_conc_B)*100
# PSS_B=conc_opt[-1,1]/(initial_conc_A+initial_conc_B)*100

# print(f"At the PSS, {PSS_A:.2f} % of A and {PSS_B:.2f} % of B")

#%%

##################################################
################# PLOT and SAVE ##################
##################################################
# ToSave = input("Save plots (Yes or No)? ")
##################################################

## Plot the experimental data, the fitted total absorbance curve, and residuals
# fig = plt.figure(figsize=(8, 4),dpi=600,constrained_layout=True)
# gs = gridspec.GridSpec(4, 2, figure=fig)
# fig.suptitle(f'{wavelength_nm} nm: SingleWavelength Method\n \
#              {datafilename}')

# # fig, ax1 = plt.subplots(2, 1, figsize=(6, 4),dpi=600,
# #                         gridspec_kw={'height_ratios': [4,1]}, 
# #                         constrained_layout=True)

# axresults_conc = fig.add_subplot(gs[0:3, 0])
# axresults_Abs = fig.add_subplot(gs[0:3, 1])
# axresults_res = fig.add_subplot(gs[3, 1])

# axresults_Abs.plot(timestamps, absorbance_values, linestyle='-', color="#DCEEFF", 
#             linewidth=5, label='Experimental Data')
# axresults_Abs.plot(timestamps, absorbance_fitted, linestyle='--', color="#9ACCFF", 
#             label=f"Fitted Total Absorbance Curve\nQY_AB: {QY_AB_opt:.2f}"+u"\u00B1"+f"{error_QY_AB:.3f}\
#                 \nQY_BA: {QY_BA_opt:.2f}"+u"\u00B1"+f"{error_QY_BA:.3f}")
# axresults_Abs.set_title(f'Absorbance at {wavelength_nm}nm')
# axresults_Abs.set_xlabel('Time (s)')
# axresults_Abs.set_ylabel('Absorbance')
# axresults_Abs.legend()

# #### residuals
# axresults_res.plot(timestamps, result_lmfit.residual, color="#FF952A", 
#             label='Residuals')
# axresults_res.set_xlabel('Time (s)')
# axresults_res.set_ylabel('Residual')

# #################################################
# #### CONCENTRATION OVER TIME ####
# axresults_conc.plot(timestamps, conc_opt[:,0], label='Species A', 
#          color="#9ACCFF")
# axresults_conc.plot(timestamps, conc_opt[:,1], label='Species B', 
#          color="#FF952A")
# axresults_conc.set_title('Concentrations over Time')
# axresults_conc.set_xlabel('Time (s)')
# axresults_conc.set_ylabel('Concentration (M)')
# axresults_conc.legend()
# # plt.grid(True)

# #### SAVE PLOTS ####
# if ToSave == "Yes":
#     savefilename = datafolder+"\\"+"QY-SingleWavelength_"+str(wavelength_nm)+"nm_"+datafilename
#     plt.savefig(savefilename+".svg",bbox_inches="tight")
#     plt.savefig(savefilename+".png",bbox_inches="tight")
#     print("Saved plots")
# elif ToSave == "No":
#     print("Plots not saved")
# else:
#     print("Wrong ToSave input")

# ####################
# plt.show()
########################################
#%%
