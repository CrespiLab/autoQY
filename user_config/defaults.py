# -*- coding: utf-8 -*-
########################### FORMAT OF LOG FILE ################################
# format_timestamps = "Default" # simple two-column format with timestamps
format_timestamps = "AHK" ## Crespi group format
###############################################################################

################# METHOD OF SOLVING DIFFERENTIAL EQUATIONS ####################
CalculationMethod = "Integration" ## Integration over LED emission spectrum

ODEMethod = "Concentrations" ## using the concentrations obtained beforehand (new method)
# ODEMethod = "Emission" ## using the part of the spectrum that overlaps with LED emission only (old method)
###############################################################################

######################### METHOD OF POWER INPUT ###############################
PowerMethod = "Manual" ## Manual input of power
# PowerMethod = "PowerProcessing" ## Enables loading of PowerProcessing module
###############################################################################

################## THRESHOLD FOR PROCESSING OF LED EMISSION ###################
threshold = 500
###############################################################################

################ BASELINE CORRECTION OF LED EMISSION SPECTRA ##################
BaselineCorrection_LED = "ON"
###############################################################################

####################### AXIS LIMITS FOR PROCESSED DATA ########################
xlim_min_ProcessedData = 220
xlim_max_ProcessedData = 760
ylim_min_ProcessedData = -0.1
ylim_max_ProcessedData = 1.25

##################### WAVELENGTH RANGE OF PROCESSED DATA#######################
wl_low = 250
wl_high = 800

####################### LED Wavelength ########################
LEDw = 340

############################## COLOURS ########################################
CrespiColours = {'black': "#000000",
           'grey': "#808080",
           'greylight': "#c0c0c0",
           'greylighter': "#dddddd",
           'greylightest': "#f2f2f2",
           'blue': "#346aa9",
           'bluelight': "#7dadce",
           'bluelighter': "#e0edf6",
           'orangered': "#872f17",
           'orange': "#e16203",
           'orangelight': "#ffd8b0",
           'green': '#429130',
           'greenlight': '#d7e5c5'}

colours = CrespiColours # set default here
colours_plot_first = colours['blue']
colours_plot_last = colours['orange']
colours_plot_grey = colours['greylighter']
###############################################################################
