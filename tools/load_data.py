# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import data.calc_settings as CalcSettings

def GetTimestamps(LogFile):
    """ Obtain timestamps from .ahk log file """
    ##!!! should add a function that checks that the format is correct
    
    ## CSV not DAT
    if CalcSettings.format_timestamps == "AHK":
        ''' Return corrected timestamps:
            - generate time intervals from time logged between turning LED on and off,
            i.e. the actual irradiation time.
            - Re-create timestamps by doing a cumulative sum of the obtained intervals,
            after the addition of timestamps 0.0 s to the start of the array
        '''
        log = pd.read_csv(LogFile, sep = ",", decimal = ".")
        log_measure=log[log[log.columns[3]] == 'Measure'] ## Measure lines in log file
        log_LEDon=log[log[log.columns[3]] == 'LEDon']
        log_LEDoff=log[log[log.columns[3]] == 'LEDoff']
        
        measure = log_measure.iloc[:, [2]].to_numpy() ## array of Measure instances
        timestamps_LEDon = log_LEDon.iloc[:, [2]].to_numpy()
        timestamps_LEDoff = log_LEDoff.iloc[:, [2]].to_numpy()
        timestamps_LEDon = timestamps_LEDon[:len(timestamps_LEDoff)]
        intervals_OffMinusOn = timestamps_LEDoff - timestamps_LEDon
        
        timestamps = np.cumsum(np.insert(intervals_OffMinusOn, 0, 0.0)) ## cumulative sum; add 0 to start
        
        if len(measure) != len(timestamps): ## remove final element in case of extra set of LEDon-LEDoff lines in log file
            timestamps = timestamps[:len(measure)]
            print(f"Cut timestamps array to length of Measure array: {len(measure)}")
        else:
            pass
        print(f"timestamps len: {len(timestamps)}")
    elif CalcSettings.format_timestamps == "Default":
        log = pd.read_csv(LogFile,
                        sep = ",", decimal = ".", skiprows = 1, header=None,)
        log.columns = ["Measurement", "Timestamp (s)"]
        timestamps = log.iloc[:,1].to_numpy()
    return timestamps

def Import_SpectralData(FileFormat, file):
    """ 
    Import Absorbance data_file 
    .dat from Spectragryph is with Wavenumbers column
    """
    if FileFormat == "Spectragryph":
        data_pd = pd.read_csv(file, sep='\t', usecols=lambda x: x not in ["Wavenumbers [1/cm]"])  # in nanometers
    elif FileFormat == "Not":
        data_pd = pd.read_csv(file, delimiter=',', header=0)  # in nanometers ####ALFREDO: FIXED?

    # print(f"Integration-Import_SpectralData===data_pd:{data_pd}")

    # data_pd.columns = ['Wavelength [nm]', 'Absorbance'] ## rename columns
    
    data_pd_full = data_pd
    # print(f"Integration-Import_SpectralData===data_pd_full:{data_pd_full}")

    data_pd_wavelengths = data_pd['Wavelength [nm]']
    # print(f"LoadData===Import_SpectralData===data_pd_wavelengths:\n{data_pd_wavelengths}")
    data_pd_absorbance = data_pd.iloc[:,1:]
    # print(f"LoadData===Import_SpectralData===data_pd_absorbance:\n{data_pd_absorbance}\n \
          # data_pd_absorbance.shape: {data_pd_absorbance.shape}")
    # print(f"data_pd_absorbance[0]: {data_pd_absorbance[0]}")

    return data_pd_full, data_pd_wavelengths, data_pd_absorbance

def Import_Epsilons(FileFormat, 
                    X):
    """ 
    Epsilons datafiles are in a certain format
    TO DO for GUI:
        - select delimiter: ',' or '\t'
        - select; ignore Wavenumbers column or not
    Returns: wavelengths and intensities as numpy.ndarray
    """
    if FileFormat == "Spectragryph":
        epsilon_data = pd.read_csv(X, delimiter='\t',
                                     usecols = lambda x: x not in ["Wavenumbers [1/cm]"])
    elif FileFormat == "Not":
        epsilon_data = pd.read_csv(X, delimiter=',', 
                                     skiprows=1, usecols=[0,1]) 
        
    epsilon_data.columns = ['Wavelengths', 'Epsilons'] ## rename columns
    epsilon_wavelengths = epsilon_data['Wavelengths'].values
    epsilon_values = epsilon_data['Epsilons'].values
    return epsilon_wavelengths, epsilon_values

def Import_LEDemission(FileFormat, file_LEDemission_raw):
    ## Import (raw) emission LED_file
    if FileFormat == "Spectragryph":
        emission_data = pd.read_csv(file_LEDemission_raw, sep = '\t', usecols = lambda x: x not in ["Wavenumbers [1/cm]"]) # in nanometers
    elif FileFormat == "Not":
        emission_data = pd.read_csv(file_LEDemission_raw, delimiter=',') # in meters

    emission_data.columns = ['Wavelength [nm]', 'Intensity'] ## rename columns

    emission_wavelengths = emission_data['Wavelength [nm]'].values
    emission_Intensity = emission_data['Intensity'].values ## not normalised
    return emission_wavelengths, emission_Intensity

def check_not_empty(thing):
    ''' Check if a list or numpy array exists (i.e. is not empty) '''
    try:
        if type(thing) == list:
            if not thing:
                # print('empty list')
                thing_exists = False
            else:
                # print("this list exists")
                thing_exists = True
        elif type(thing) == np.ndarray:
            if thing.size == 0:
                # print('empty numpy array')
                thing_exists = False
            else:
                # print("this numpy array exists")
                thing_exists = True
        else:
            print("this thing is something else")
        return thing_exists
    except Exception as e:
        print(f"FAILED to check whether list or numpy array exists: {e}") ##!!! send this to message console somehow (need to make this module a class?)