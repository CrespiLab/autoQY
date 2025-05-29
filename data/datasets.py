# -*- coding: utf-8 -*-
"""
@author: jorst136
"""

I0_list = []
initial_conc_R = None
initial_conc_P = None
wavelengths_meters = None
normalized_emission = None
N = None
fit_results = []

def ListPowers(I0_avg,I0_err):
    return [I0_avg, I0_avg+I0_err, I0_avg-I0_err]