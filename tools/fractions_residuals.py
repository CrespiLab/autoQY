# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 2024
@author: Jorn Steen
Window for showing residuals of fractions
"""
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

from tools.plotting import MplCanvas
import data.loaded_data as LoadedData
import data.results as Results

from UIs.WindowFractionsResiduals import Ui_MainWindow

class WindowFractionsResiduals(QtWidgets.QMainWindow, Ui_MainWindow):
    """Class for FractionsResiduals module."""
    def __init__(self, parent):
        self.parent = parent
        super(WindowFractionsResiduals, self).__init__()
        
        # progress = pyqtSignal(str) ##!!! TRY TO SEND SIGNAL TO message_console
        ### print(f"number of plots: {num_plots}")
        
        self.setupUi(self)
        
        ##!!! ADD toolbar for zoom

        ## Generate several plots
        num_plots = Results.original_spectra.shape[0]
        for i in range(num_plots):
            canvas = self.create_plot(i,
                                      LoadedData.SpectralDataCut_Wavelengths,
                                      Results.original_spectra[i],
                                      Results.reconstructed_spectra_fractions_Abs[i],
                                      Results.fractions_residuals[i])
            canvas.setMinimumHeight(50)  # ensures proper spacing in scroll area
            self.verticalLayout.addWidget(canvas)

    def create_plot(self, index, wl, original, reconstructed, residuals):
        ##!!! MOVE TO plotting.py
        ''' 
        Original and reconstructed spectra on the top
        Residuals (difference between original and reconstructed spectra) on the bottom
        Creates one plot.
        '''
        fig = Figure(figsize=(8,4), dpi=200, constrained_layout=True)
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1], figure=fig)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        spectrum_number = index+1
        ax1.set_title(f"Spectrum {spectrum_number}")

        ax1.plot(wl,original, label="Original")
        ax1.plot(wl,reconstructed, label="Reconstructed from Fit")
        ax2.plot(wl,residuals, label="Original - Reconstructed")
        
        ax1.set_xticklabels([]) ## put this before set_xlim, otherwise it resets the xlim
        ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
        
        for i in [ax1,ax2]:
            i.yaxis.set_minor_locator(AutoMinorLocator(2))
        
        ax1.set_ylabel("Absorbance")
        ax2.set_ylabel("Residuals")
        ax2.set_xlabel("Wavelength (nm)")

        ax1.legend()
        ax2.legend()

        canvas = FigureCanvas(fig)
        
        return canvas
