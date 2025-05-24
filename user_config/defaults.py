# -*- coding: utf-8 -*-
"""
Created on Fri May 23 19:29:57 2025

@author: Jorn Steen
"""

########################### FORMAT OF LOG FILE ################################
# format_timestamps = "Default" # simple two-column format with timestamps
format_timestamps = "AHK" ## Crespi group format
###############################################################################

################# METHOD OF SOLVING DIFFERENTIAL EQUATIONS ####################
# CalculationMethod = "SingleWavelength" ## Using a single wavelength (TBA)
CalculationMethod = "Integration" ## Integration over LED emission spectrum
###############################################################################

######################### METHOD OF POWER INPUT ###############################
PowerMethod = "Manual" ## Manual input of power
# PowerMethod = "PowerProcessing" ## Enables loading of PowerProcessing module
###############################################################################

################## THRESHOLD FOR PROCESSING OF LED EMISSION ###################
threshold = 500
###############################################################################