import pandas as pd
import data.calc_settings as CalcSettings

# def GetTimestamps(LogFile):
def GetTimestamps(LogFile):
    """ Obtain timestamps from .ahk log file """
    ##!!! should add a function that checks that the format is correct
    
    ## CSV not DAT
    if CalcSettings.format_timestamps == "AHK":
        #print("Start of load_data.GetTimestamps")
        log = pd.read_csv(LogFile,
                        sep = ",", decimal = ".", skiprows = 1, header=None,)
        #print(f"log:\n{log}")
        log_measure=log[log[log.columns[3]] == 'Measure']
        measurement_timestamps = log_measure.iloc[:, [0, 2]]
        measurement_timestamps.columns = ["Measurement", "Timestamp (s)"]
        timestamps = measurement_timestamps.iloc[:,1].to_numpy()
    elif CalcSettings.format_timestamps == "Default":
        log = pd.read_csv(LogFile,
                        sep = ",", decimal = ".", skiprows = 1, header=None,)
        log.columns = ["Measurement", "Timestamp (s)"]
        timestamps = log.iloc[:,1].to_numpy()
    # else:

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
