# -*- coding: utf-8 -*-
"""
Fit a UV-Vis mixture spectrum using two reference spectra on different wavelength grids.

Input:
- Three CSV/TSV files with two columns each:
  1. mixture.csv:    λ, A_mix
  2. eps1.csv:       λ, ε₁
  3. eps2.csv:       λ, ε₂

Assumptions:
- Beer–Lambert law applies: A_mix = α·ε₁ + β·ε₂.

Options
-------
--poly N   Polynomial degree N for baseline correction (default 2).  Use 0 to
           disable baseline correction entirely.
"""

import numpy as np
from scipy.optimize import nnls

import data.results as Results

def polynomial_baseline(x: np.ndarray, y: np.ndarray,
                        deg: int) -> np.ndarray:
    """Return polynomial baseline of degree `deg` (least-squares)."""
    if deg <= 0:
        return np.zeros_like(y)
    coeffs = np.polyfit(x, y, deg)
    return np.polyval(coeffs, x)

def CalculateFractions(spectra, wavelengths, eps1, eps2):
    '''
    Load absorption spectra recorded during irradiation, and epsilon spectra of both species.
    '''
    matrix_Abs_spectra = spectra.T  # shape = (n_spectra, n_points)
    matrix_epsilons = np.column_stack([eps1, eps2])  # shape = (n_points, 2)
    
    # Fit each spectrum
    f1_list, f2_list = [], []
    
    for i, spectrum in enumerate(matrix_Abs_spectra):
        # baseline = polynomial_baseline(wavelengths, spectrum, 0) ## polynomial: 0 (off)
        # A_corr = spectrum - baseline
        # coeffs, _ = nnls(matrix_epsilons, A_corr)

        ## no baseline correction
        coeffs, _ = nnls(matrix_epsilons, spectrum) 
        total = coeffs.sum()
        frac1 = coeffs[0] / total if total > 0 else 0.0
        frac2 = coeffs[1] / total if total > 0 else 0.0

        f1_list.append(frac1)
        f2_list.append(frac2)
    
    reconstructed_spectra = []
    
    for f1, f2 in zip(f1_list, f2_list):
        fit = f1 * eps1 + f2 * eps2
        reconstructed_spectra.append(fit)
    reconstructed_spectra = np.array(reconstructed_spectra)  # shape = (n_spectra, n_points)

    original_spectra = matrix_Abs_spectra # shape = (n_spectra, n_points)
    return f1_list, f2_list, reconstructed_spectra, original_spectra

def CalculateResiduals(original, reconstructed_epsilon, total_conc):
    ''' 
    Convert reconstructed spectra from epsilons to Abs using obtained total concentration
    Calculate residuals in Abs: Original - Reconstructed
    '''
    reconstructed_Abs = reconstructed_epsilon * total_conc
    residuals = original - reconstructed_Abs
    return reconstructed_Abs, residuals
    