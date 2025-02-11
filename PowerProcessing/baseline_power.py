# import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fit_func(x, *coeffs):  #not used
    return np.polyval(coeffs, x)

def choose_sections(line_positions):
    """
    Choose sections based on the line positions provided.
    This function expects a list of line positions (integers).
    """
    # Ensure enough positions are provided
    if len(line_positions) < 12:
        raise ValueError("Not enough positions selected. Please ensure you have defined 12 positions.")

    # Convert positions to integers
    line_positions = [int(pos) for pos in line_positions]

    # Assign positions to different sections
    sections = {
        "start_0": line_positions[0],
        "end_0": line_positions[1],
        "start_1": line_positions[2],
        "end_1": line_positions[3],
        "start_0_1": line_positions[4],
        "end_0_1": line_positions[5],
        "start_2": line_positions[6],
        "end_2": line_positions[7],
        "start_3": line_positions[8],
        "end_3": line_positions[9],
        "start_4": line_positions[10],
        "end_4": line_positions[11],
    }
    return sections


def baseline_correction(RefPower, x, sections):
    """
    Perform baseline correction for the selected sections without jacket and cuvette.
    """
    # Ensure that the sections contain valid indices (integers)
    if not all(isinstance(sections[key], int) for key in sections):
        raise ValueError("Section indices must be integers.")
    
    # Concatenate sections for baseline correction
    try:
        x_masked = np.concatenate((x[sections["start_0"]:sections["end_0"]], 
                                   x[sections["start_0_1"]:sections["end_0_1"]]))
        y_masked = np.concatenate((RefPower[sections["start_0"]:sections["end_0"]], 
                                   RefPower[sections["start_0_1"]:sections["end_0_1"]]))
    except KeyError as e:
        print(f"Section key error: {e}")
        return None, None
    
    # Remove NaN values
    x_masked = x_masked[~np.isnan(y_masked)]
    y_masked = y_masked[~np.isnan(y_masked)]

    # Fit polynomial to the baseline
    n = 3  # Degree of polynomial for baseline correction
    p0 = np.full(n + 1, 0.000000001)  # Initial guess of baseline coefficients
    popt, _ = curve_fit(fit_func, x_masked, y_masked, p0=p0)

    # Generate baseline and baseline-corrected data
    baseline = np.polyval(popt, x)
    baselined = RefPower - baseline

    return baselined, baseline