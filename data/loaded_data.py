# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 17:21:49 2024

@author: jorst136
"""
count = 0

line_positions = [] 

PowersAtCuvette = {}
ErrorsAtCuvette = {}

filename_LED = None

epsilons_R_wavelengths = None 
epsilons_R_values = None
epsilons_P_wavelengths = None
epsilons_P_values = None

epsilons_R_interp = None
epsilons_P_interp = None
emission_interp = None

timestamps = None

SpectralData_Full = None
SpectralData_Wavelengths = None
SpectralData_Absorbance = None

SpectralDataCut_Wavelengths = None
SpectralDataCut_Abs = None
SpectralDataCut_Index = None

SpectralData_AbsAtLEDw = None

LEDemission_wavelengths = None
LEDemission_intensity = None
LEDemission_intensity_proc = None

######## EPSILONS FOR SINGLEWAVELENGTH METHOD########
epsilon_R = 0
epsilon_P = 0

