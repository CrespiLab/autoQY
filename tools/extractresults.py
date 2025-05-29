# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:59:07 2024

@author: jorst136
"""
import numpy as np
from scipy.integrate import odeint
from QY.integration import rate_equations

def ExtractResults(fit_results):
    """ Extract optimised QYs and errors """
    result_lmfit_avg=fit_results[0] ## results using I0_avg
    result_lmfit_avgpluserr=fit_results[1] ## results using I0_avg+I0_error
    result_lmfit_avgminerr=fit_results[2] ## results using I0_avg-I0_error

    ##!!! CHANGE TO PROPAGATING THE ERRORS: std. devs and power errors
    
    ##!!! ALSO ADD THE ERROR ON THE EPSILONS

    ## Extract the optimized parameters and their standard deviations
    QY_AB_opt_avg = result_lmfit_avg.params['QY_AB'].value ## value using I0_avg
    std_dev_fit_QY_AB = result_lmfit_avg.params['QY_AB'].stderr ## error in the fit using I0_avg
    
    QY_BA_opt_avg = result_lmfit_avg.params['QY_BA'].value ## value using I0_avg
    std_dev_fit_QY_BA = result_lmfit_avg.params['QY_BA'].stderr ## error in the fit using I0_avg
    
    QY_AB_opt_min=result_lmfit_avgpluserr.params['QY_AB'].value \
                        - result_lmfit_avgpluserr.params['QY_AB'].stderr ## results using I0_avg+I0_err
    QY_AB_opt_max=result_lmfit_avgminerr.params['QY_AB'].value \
                        + result_lmfit_avgminerr.params['QY_AB'].stderr ## results using I0_avg-I0_err

    QY_BA_opt_min=result_lmfit_avgpluserr.params['QY_BA'].value \
                        - result_lmfit_avgpluserr.params['QY_BA'].stderr ## results using I0_avg+I0_err
    QY_BA_opt_max=result_lmfit_avgminerr.params['QY_BA'].value \
                        + result_lmfit_avgminerr.params['QY_BA'].stderr ## results using I0_avg-I0_err

    
    ## Calculate error for QY_AB and QY_BA based on the error in the power 
        ## The error in the fit is considered in the definition of the min and max values (see above)
    error_QY_AB = max([QY_AB_opt_avg - QY_AB_opt_min, QY_AB_opt_max - QY_AB_opt_avg])
    error_QY_BA = max([QY_BA_opt_avg - QY_BA_opt_min, QY_BA_opt_max - QY_BA_opt_avg])
    
    ## Print the results
    print(f"Optimized QY_AB: {QY_AB_opt_avg:.5f}" )
    print(f"Standard deviation of QY_AB in the fit using I0_avg: {std_dev_fit_QY_AB:.5f}")
    print(f"Error for QY_AB: {error_QY_AB:.5f} ")
    
    print(f"Optimized QY_BA: {QY_BA_opt_avg:.5f}" )
    print(f"Standard deviation of QY_BA in the fit using I0_avg: {std_dev_fit_QY_BA:.5f}")
    print(f"Error for QY_BA: {error_QY_BA:.5f} ")
    
    return QY_AB_opt_avg, QY_BA_opt_avg, error_QY_AB, error_QY_BA

def CalculateConcentrations(lambda_meters, 
                            init_conc_A, init_conc_B,
                            timestamps,
                            QY_AB_opt, QY_BA_opt,
                            e_A_inter, e_B_inter,
                            N, V):
    '''
    Calculate concentrations over time using optimised QYs
    '''
    
    ## Integrate the rate equations with the optimized parameters
    conc_opt = odeint(rate_equations, 
                      [init_conc_A, init_conc_B], 
                      timestamps,
                      args=(lambda_meters, QY_AB_opt, QY_BA_opt, 
                            e_A_inter, e_B_inter, 
                            N, V))
        
    PSS_A = conc_opt[-1,0]/(init_conc_A+init_conc_B)*100
    PSS_B = conc_opt[-1,1]/(init_conc_A+init_conc_B)*100
    
    print(f"At the PSS, {PSS_A:.2f} % of Reactant and {PSS_B:.2f} % of Product")
    return conc_opt, PSS_A, PSS_B

def GetFittedAbs(fit_results, conc_opt,
                 e_A_inter, e_B_inter,
                 timestamps,
                 wavelengths_data):
    ''' 
    Get the fitted absorbance
    '''
    result_lmfit=fit_results[0] ## results using I0_avg
    total_abs_fit = conc_opt.dot(np.vstack([e_A_inter, e_B_inter]))
    residuals = result_lmfit.residual.reshape((len(timestamps), 
                                              len(wavelengths_data))).T
    return total_abs_fit, residuals
##################################################
