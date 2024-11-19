========================================
Based on the publication: 
A. Volker, J. D. Steen, S. Crespi, A fiber-optic spectroscopic setup for isomerization quantum yield determination, Beilstein J. Org. Chem. 2024, 20, 1684–1692, DOI: 10.3762/bjoc.20.150.

========================================
Default data format:

Wavelength (nm)	Timestamp_0	Timestamp_1	Timestamp_2
wl_1	Abs_1	Abs_1	Abs_1
wl_2	Abs_2	Abs_2	Abs_2

========================================

The PowerProcessing module is used for processing the power data (e.g., measured with the Thorlabs software Optical Power Monitor (OPM))
Workflow:
0) Select the x-coordinates of the areas to be used for subsequent baseline correction (interactive plot)
1) Baseline correction using a polynomial
2) Retrieve value+standard deviation
3) Calculate power at cuvette
4) Calculate power average and error over triplicate measurements

The power data should be recorded in the following sequence:
1) no insulating jacket, no cuvette -- LED off
2) no insulating jacket, no cuvette -- LED on
3) no insulating jacket, no cuvette -- LED off
4) insulating jacket, cuvette with solvent -- LED off
5) insulating jacket, cuvette with solvent -- LED on
6) insulating jacket, cuvette with solvent -- LED off
