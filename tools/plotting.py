import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LinearSegmentedColormap

import data.experimental_parameters as ExpParams
import data.calc_settings as CalcSettings
import user_config.defaults as Defaults

class MplCanvas(FigureCanvas):
    """Widget class to render plots."""
    def __init__(self, parent=None, idx=None):
        self.parent_window = parent
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)
        self.setParent(parent)
        #self.lines = [] 
        #self.boxes = []
        self.selected_line = None
        self.press_event = None
        self.line_positions = []

        self.colours = ["#000000","#808080","#346aa9","#7dadce","#e16203","#872f17"]

        self.idx = idx  # Store the idx for this canvas

        self.mpl_connect('button_press_event', self.on_press)
        self.mpl_connect('motion_notify_event', self.on_motion)
        self.mpl_connect('button_release_event', self.on_release)

    def plot_power(self, x, RefPower, filename):
        self.ax.plot(x, RefPower,alpha=0.5,color='black',label="Power")
        self.ax.set_ylabel("Power (mW)")
        self.ax.set_xlabel("Index")
        self.ax.set_title(filename.split("/")[-1])

        #############################
        ## Adding initial vertical lines
        numberoflines = 12
        red = ['red']
        blue = ['blue']
        colours = red*2+blue*2+red*4+blue*2+red*2
        self.lines = []
        for i in range(0,numberoflines):
            self.lines.append(self.ax.axvline(x=x[-1]/numberoflines*i,
                                              alpha=0.2, color=colours[i],
                                              picker=5))
        ## By setting the picker property on each line (picker=5), we make the lines selectable when clicked near them.
        #############################
        self.draw()

        ## Initialize line positions
        self.line_positions = [line.get_xdata()[0] for line in self.lines]
        #############################
        ## initalize coloured boxes
        self.boxes = []
        self.update_boxes()

    def update_boxes(self):
        """Update the shaded areas (boxes) between pairs of lines based on updated line positions."""
        # Clear existing boxes
        for box in self.boxes:
            box.remove()
        self.boxes = []

        red = ['red']
        blue = ['blue']
        colours = red + blue + red*2 + blue + red  # Extend colors if needed

        ## Ensure line_positions is initialized for the given idx
        if self.idx is not None:
            if len(self.parent_window.line_positions) <= self.idx:
                # Initialize with default positions if idx is out of range
                self.parent_window.line_positions.append([0] * 12)
            
            line_positions = self.parent_window.line_positions[self.idx]

        else:
            line_positions = self.line_positions

        ## Create new boxes based on updated line positions
        for i in range(0, len(line_positions) - 1, 2):
            box = self.ax.axvspan(line_positions[i], 
                                line_positions[i+1], 
                                color=colours[int(i/2)], alpha=0.1)
            self.boxes.append(box)

        ## Redraw the canvas to reflect changes
        self.draw()

    def on_press(self, event):
        if event.inaxes != self.ax:
            return

        ## Check if the mouse is near any line
        for line in self.lines:
            if line.contains(event)[0]:
                self.selected_line = line
                self.press_event = event
                break

    def on_motion(self, event):
        if self.selected_line is None or event.inaxes != self.ax:
            return

        # Calculate the new x position for the line
        dx = event.xdata - self.press_event.xdata
        new_x = self.selected_line.get_xdata()[0] + dx
        ## Update the line's position with a sequence of two identical values
        self.selected_line.set_xdata([new_x, new_x])
        self.press_event = event
        self.ax.figure.canvas.draw()

        # Update the corresponding label and line position variable with the new x position
        line_index = self.lines.index(self.selected_line)
        self.parent_window.update_display(self.idx, line_index, new_x)  # Pass self.idx here
        
        self.parent_window.line_positions[self.idx][line_index] = new_x  # Update using self.idx
        self.line_positions[line_index] = new_x

        ## Update the boxes
        self.update_boxes()

    def on_release(self, event):
        self.selected_line = None
        self.press_event = None


    def plot_baseline_correction(self, x, RefPower, baseline, baselined, sections,
                                 case, filename):
        """Plot baseline-corrected power for both cases."""
        self.fig.clear()  # Clear the entire figure
        
        # Use GridSpec to create a grid layout with two rows
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.4)  # Add space between the plots
        if case == 0:
            # First plot (ax1) will be the original power and baseline
            ax1 = self.fig.add_subplot(gs[0])

            ax1.set_title(f'{filename.split("/")[-1]}\nNo Jacket and No Cuvette')

            # Plot original power and baseline for no jacket, no cuvette
            ax1.plot(x[sections["start_0"]:sections["start_2"]],
                    RefPower[sections["start_0"]:sections["start_2"]],
                    linestyle='-', color="black", label='Power')
            ax1.plot(x[sections["start_0"]:sections["start_2"]],
                    baseline[sections["start_0"]:sections["start_2"]],
                    linestyle='--', color="red", label='Baseline')

            ax1.set_xlabel('Index')
            ax1.set_ylabel('Power (mW)')
            ax1.legend()

            # Second plot (ax2) will be the baseline-corrected power
            ax2 = self.fig.add_subplot(gs[1])
            ax2.plot(x[sections["start_0"]:sections["start_2"]], 
                    baselined[sections["start_0"]:sections["start_2"]],
                    color="#FF952A", label='Power, baseline corrected')
            ax2.set_xlabel('Index')
            ax2.set_ylabel('Power (mW)')
            ax2.legend()

        elif case == 1:
            ax1 = self.fig.add_subplot(gs[0])
            ax1.set_title(f'{filename.split("/")[-1]}\nWith Jacket and Cuvette with Solvent')

            # Plot original power and baseline for jacket, cuvette with solvent
            ax1.plot(x[sections["end_1"]:sections["end_4"]],
                    RefPower[sections["end_1"]:sections["end_4"]],
                    linestyle='-', color="black", label='Power')
            ax1.plot(x[sections["end_1"]:sections["end_4"]],
                    baseline[sections["end_1"]:sections["end_4"]],
                    linestyle='--', color="red", label='Baseline')

            ax1.set_xlabel('Index')
            ax1.set_ylabel('Power (mW)')
            ax1.legend()

            # Second plot (ax2) will be the baseline-corrected power
            ax2 = self.fig.add_subplot(gs[1])
            ax2.plot(x[sections["end_1"]:sections["end_4"]],
                    baselined[sections["end_1"]:sections["end_4"]],
                    color="#FF952A", label='Power, baseline corrected')
            ax2.set_xlabel('Index')
            ax2.set_ylabel('Power (mW)')
            ax2.legend()

        # Redraw the figure and adjust layout to prevent overlap
        self.fig.tight_layout()
        self.draw()


    def dynamic_range_x(self, wavelengths):
        ''' Dynamic x-axis range with some padding '''
        x_min, x_max = wavelengths.min(), wavelengths.max()
        x_padding = 0.05 * (x_max - x_min)  # 5% padding
        x_min_dynamic, x_max_dynamic = x_min - x_padding, x_max + x_padding
        return x_min_dynamic, x_max_dynamic
    
    def dynamic_range_y(self, absorbance):
        ''' Dynamic y-axis range with some padding '''
        y_min, y_max = absorbance.min(), absorbance.max()
        y_padding = 0.05 * (y_max - y_min)  # 5% padding
        y_min_dynamic, y_max_dynamic = y_min - y_padding, y_max + y_padding
        return y_min_dynamic, y_max_dynamic
    

    def plot_EpsilonsOnly(self, eR_wavelengths, eR_Abs, eP_wavelengths,eP_Abs):
        """Plot epsilons before interpolation."""
        self.fig.clear()  # Clear the entire figure
        self.fig.set_constrained_layout(True)
        
        # Use GridSpec to create a grid layout
        # gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.5)  # Two rows, evenly spaced, with some space between
        gs = gridspec.GridSpec(1, 1, figure=self.fig) 

        # Create first plot (Epsilons) in the first row
        ax1 = self.fig.add_subplot(gs[0])
        
        # Plot epsilons on the first subplot
        ax1.plot(eR_wavelengths, eR_Abs, label="Reactant",
                 color = Defaults.colours_plot_first)
        ax1.plot(eP_wavelengths, eP_Abs, label="Product",
                 color = Defaults.colours_plot_last)
        ax1.legend(fontsize=12)
        ax1.set_xlabel("Wavelength (nm)")
        # ax1.set_xlim(220, 650)
        ax1.set_title("Epsilons")
        ax1.set_ylabel(r"$\epsilon$ (M$^{-1}$ cm$^{-1}$)")

        self.draw()  # Redraw the canvas

    def PlotData_Full(self, wavelengths, absorbance):
        """Plot interpolated epsilons and LED emission."""
        self.fig.clear()  # Clear the entire figure
        self.fig.set_constrained_layout(True)
        
        # Use GridSpec to create a grid layout
        gs = gridspec.GridSpec(1, 1, figure=self.fig) 

        # Create first plot (Epsilons) in the first row
        ax1 = self.fig.add_subplot(gs[0]) 

        cmap = LinearSegmentedColormap.from_list('mylist',
                                                 [(0, Defaults.colours_plot_first),
                                                  (1, Defaults.colours_plot_last)])
        num_plots = absorbance.shape[1]
        colours_formap = cmap(np.linspace(0,1,num_plots))
        
        ## Plot spectra
        for i in range(0,len(absorbance.columns)):
            ax1.plot(wavelengths, absorbance[absorbance.columns[i]],
                     color = colours_formap[i])

        ax1.set_xlabel("Wavelength (nm)")
        # ax1.set_xlim(220, 650)
        # ax1.set_ylim(-0.05, 1.2)
        ax1.set_title("Spectra during Irradiation")
        ax1.set_ylabel(r"Absorbance")

        line_first = mlines.Line2D([], [], color=Defaults.colours_plot_first, marker='',
                                  markersize=15, label='Before irr.')
        line_last = mlines.Line2D([], [], color=Defaults.colours_plot_last, marker='',
                                  markersize=15, label='PSS')
        ax1.legend(handles=[line_first, line_last],
                   loc='best', fontsize=12)

        self.draw()  # Redraw the canvas

    def plot_LEDemission_full(self, LED_wavelength,LED_intensity):
        """Plot epsilons before interpolation."""
        self.fig.clear()  # Clear the entire figure
        self.fig.set_constrained_layout(True)
        
        # Use GridSpec to create a grid layout
        gs = gridspec.GridSpec(1, 1, figure=self.fig) 

        # Create plot
        ax1 = self.fig.add_subplot(gs[0])
        
        # Plot 
        ax1.plot(LED_wavelength, LED_intensity, 
                 label=f"{ExpParams.LEDw} nm",
                 color = Defaults.colours_plot_first)
        ax1.legend(fontsize=12)
        ax1.set_xlabel("Wavelength (nm)")
        # ax1.set_xlim(220, 650)
        ax1.set_title("LED emission")
        ax1.set_ylabel(r"Intensity")

        self.draw()  # Redraw the canvas

    def PlotData_Cut(self, absorbance_full, wavelengths_full,
                     absorbance_cut, wavelengths_cut,
                     em_wl, em_int,
                     # index_first, index_last,
                     em_int_proc):
        """Data (epsilons and spectra) cut according to LED Emission spectrum"""
        self.fig.clear()  # Clear the entire figure
        self.fig.set_constrained_layout(True)
        
        # Create a grid for the subplots
        gs = plt.GridSpec(2, 1, figure=self.fig)

        cmap = LinearSegmentedColormap.from_list('mylist',
                                                 [(0, Defaults.colours_plot_first),
                                                  (1, Defaults.colours_plot_last)])
        num_plots = absorbance_cut.shape[1]
        colours_formap = cmap(np.linspace(0,1,num_plots))

        ###########################################
        ###### Spectral data unprocessed and processed ######
        ###########################################

        # Add subplots using self.fig.add_subplot and GridSpec
        ax1 = self.fig.add_subplot(gs[0])

        for i in range(0,len(absorbance_full.columns)): ## pandas dataframe
            ax1.plot(wavelengths_full, absorbance_full[absorbance_full.columns[i]],
                     color = Defaults.colours_plot_grey)

        ## plot cut spectra
        for i in range(0,absorbance_cut.shape[1]): ## numpy array
            ax1.plot(wavelengths_cut, absorbance_cut.T[i],
                     color = colours_formap[i])

        ########### INSET ############
        ax1_inset = ax1.inset_axes([0.62,0.22,0.35,0.7])
        for i in range(0,absorbance_cut.shape[1]): ## numpy array
            ax1_inset.plot(wavelengths_cut, absorbance_cut.T[i],
                     color = colours_formap[i])

        inset_x_min, inset_x_max = self.dynamic_range_x(wavelengths_cut)
        inset_y_min, inset_y_max = self.dynamic_range_y(absorbance_cut)
        ax1_inset.set_xlim(inset_x_min, inset_x_max)
        ax1_inset.set_ylim(inset_y_min, inset_y_max)
        ###############################
        xlim_min = wavelengths_full.values[0] ## xlimits of full absorption spectra
        xlim_max = wavelengths_full.values[-1]
        
        # Set common x-axis properties for both subplots
        ax1.set_xlabel("Wavelength (nm)")
        ax1.set_xlim(xlim_min, xlim_max)

        # Customize the first subplot (Spectral Data)
        ax1.set_title("Spectral Data")
        ax1.set_ylabel("Absorbance")

        line_first = mlines.Line2D([], [], color=Defaults.colours_plot_first, marker='',
                                  markersize=15, label='Before irr.')
        line_last = mlines.Line2D([], [], color=Defaults.colours_plot_last, marker='',
                                  markersize=15, label='PSS')
        ax1.legend(handles=[line_first, line_last],
                   loc='best')

        ################################################################
        ############ LED emission unprocessed and processed ############
        ################################################################
        ax2 = self.fig.add_subplot(gs[1])
        ax2.set_xlabel("Wavelength (nm)")
        ax2.set_xlim(xlim_min, xlim_max)

        ## Plot LED emission data in the second subplot (ax2)
        if CalcSettings.BaselineCorrection_LED == "ON":
            bl_corr = "baseline correction applied"
        elif CalcSettings.BaselineCorrection_LED == "OFF":
            bl_corr = "no baseline correction applied"
        
        ax2.plot(em_wl, em_int, label="Untreated",
                 color = Defaults.colours_plot_first)
        ax2.plot(wavelengths_cut, em_int_proc, ## em_int_proc is the interpolated (and appropriately cut) LED emission data
                label=f"Smoothed,\n{bl_corr},\nremoved negative values,\ncut to appropriate range",
                color = Defaults.colours_plot_last)
        
        ax2.legend(fontsize=8)
        ax2.set_title(f"LED Emission ({ExpParams.LEDw} nm)")
        ax2.set_ylabel("Intensity")

        self.draw()  # Redraw the canvas


    def PlotFractions(self, wavelengths, spectra, 
                      reconstructed_spectra,
                      fractions_1, fractions_2,):
        ''' Plot fractions of Reactant and Product obtained from calculation
        '''
        self.fig.clear()  # Clear the entire figure
        self.fig.set_constrained_layout(True)

        gs = gridspec.GridSpec(2, 2, figure=self.fig)
        cmap = LinearSegmentedColormap.from_list('mylist',
                                                 [(0, Defaults.colours_plot_first),
                                                  (1, Defaults.colours_plot_last)])
        
        A_matrix_raw = spectra.T  # shape = (n_spectra, n_points)
        num_plots = len(A_matrix_raw)
        colours_formap = cmap(np.linspace(0,1,num_plots))
        
        # Top-left: Original spectra
        ax1 = self.fig.add_subplot(gs[0, 0])
        for i in range(len(fractions_1)):
            ax1.plot(wavelengths, A_matrix_raw[i], 
                     # label=f"{spectrum_names[i]}", 
                     color = colours_formap[i])
        ax1.set_title("Original Spectra")
        ax1.set_xlabel("Wavelength (nm)")
        ax1.set_ylabel("Absorbance")
        # ax1.set_xlim(xlim_left,xlim_right)
        # ax1.set_ylim(ylim_Abs_low,ylim_Abs_high)
        # ax1.grid(True)
        # ax1.legend(fontsize="x-small")
        
        # Top-right: Reconstructed spectra
        ax2 = self.fig.add_subplot(gs[0, 1])
        for i in range(len(fractions_1)):
            ax2.plot(wavelengths, reconstructed_spectra[i], 
                     # label=f"{spectrum_names[i]}", 
                     color = colours_formap[i])
        ax2.set_title("Reconstructed Spectra from Fit")
        ax2.set_xlabel("Wavelength (nm)")
        ax2.set_ylabel("Absorbance")
        # ax2.set_xlim(xlim_left,xlim_right)
        # ax2.grid(True)
        # ax2.legend(fontsize="x-small")
        
        # Bottom-left: Fractions vs spectrum index
        ax3 = self.fig.add_subplot(gs[1, 0])
        x = np.arange(len(fractions_1))
        ax3.plot(x, fractions_1, "o-", label="Reactant",
                 color = Defaults.colours_plot_first)
        ax3.plot(x, fractions_2, "o-", label="Product",
                 color = Defaults.colours_plot_last)
        
        ax3.set_title("Fitted Fractions per Spectrum")
        ax3.set_xlabel("Spectrum Index")
        ax3.set_ylabel("Fraction")
        ax3.set_ylim(-0.05, 1.05)
        # ax3.grid(True)
        ax3.legend()
        
        # Hide bottom-right subplot (not used)
        self.fig.add_subplot(gs[1, 1]).axis('off')
        
        # self.fig.tight_layout()
        self.draw() # Redraw the canvas

    def PlotFractionsResiduals(wl, original, reconstructed, residuals):
        ''' 
        Original and reconstructed spectra on the top
        Residuals (difference between original and reconstructed spectra) on the bottom
        Creates one plot.
        '''
        ##!!! ADD FUNCTION FROM fractions_residuals.py HERE


    def PlotResults(self, 
                    LEDwavelength,
                    timestamps, 
                    conc_opt,
                    absorbance,
                    index,
                    total_abs_fit,
                    residuals,
                    QY_AB_opt, QY_BA_opt, 
                    error_QY_AB, error_QY_BA,
                    ODEMethod,
                    SaveResults = "No",
                    SaveFileName = ""):
        
        ######################
        self.fig.clear()  # Clear the entire figure

        if SaveResults == "Yes": # re-define fig for savefig
            self.fig = Figure(figsize=(8, 4),dpi=600,constrained_layout=True)
            gs = gridspec.GridSpec(4, 2, figure=self.fig)
        else:
            # gs = gridspec.GridSpec(4, 2, figure=self.fig, hspace = 0.4, wspace=0.3)
            gs = gridspec.GridSpec(4, 2, figure=self.fig)
            self.fig.set_constrained_layout(True)
            self.fig.suptitle(f'{LEDwavelength} nm irradiation: {ODEMethod} Method')
        
        ##################################################
        axresults_conc = self.fig.add_subplot(gs[0:3, 0])
        axresults_Abs = self.fig.add_subplot(gs[0:3, 1])
        axresults_res = self.fig.add_subplot(gs[3, 1])
        axresults_notes = self.fig.add_subplot(gs[3, 0])
        
        ##################################################
        ############### CONCENTRATIONS ###################
        ##################################################
        axresults_conc.set_title("Conversion over time")
        
        full_conc=conc_opt[0,0]+conc_opt[0,1]
        percentage_reactant=conc_opt[:,0]/full_conc*100
        percentage_product=conc_opt[:,1]/full_conc*100
        
        axresults_conc.plot(timestamps, percentage_reactant, 
                            color = Defaults.colours_plot_first,
                            label = "Reactant")
        axresults_conc.plot(timestamps, percentage_product, 
                            color = Defaults.colours_plot_last,
                            label = "Product")
        axresults_conc.set_ylabel("Fraction (%)")
        axresults_conc.set_xlabel("Time (s)")
        
        ##################################################
        ###################### NOTES #####################
        axresults_notes.axis("off")
        ##################################################
        
        ##################################################
        ################## ABSORBANCE ####################
        ##################################################
        #### Plot the experimental data and the fitted total absorbance curve
        axresults_Abs.plot(timestamps, absorbance[index,:], linestyle='-', 
                    color=Defaults.colours_plot_first, linewidth=8, alpha=0.5, 
                    label='Experimental Data')
        axresults_Abs.plot(timestamps, total_abs_fit[:,index], linestyle='--', 
                           color=Defaults.colours_plot_last, 
                           label="Fit:\n"+r"$QY_{fwd}$"+f": {QY_AB_opt}% "+u"\u00B1 "+f"{error_QY_AB}%"+"\n"+
                           r"$QY_{bwd}$"+f": {QY_BA_opt}% "+u"\u00B1 "+f"{error_QY_BA}%")

        axresults_Abs.set_title('Fit vs Experimental')
        axresults_Abs.set_xticklabels([]) ## put this before set_xlim, otherwise it resets the xlim
        axresults_Abs.set_ylabel(r"Absorbance at $\lambda$"f"$_{{{ExpParams.LEDw}}}$")
        
        #################### RESIDUALS ###################
        axresults_res.plot(timestamps, residuals[index,:], color=Defaults.colours_plot_last, label='Residuals')
        axresults_res.set_xlabel('Time (s)')
        axresults_res.set_ylabel('Residuals Abs')
        ##################################################
        
        for i in [axresults_conc, axresults_Abs, axresults_res]:
            i.xaxis.set_minor_locator(AutoMinorLocator(2))
            i.yaxis.set_minor_locator(AutoMinorLocator(2))

        for i in [axresults_conc, axresults_Abs]:
            # i.legend(loc='upper right',frameon=False)
            i.legend(loc='best',frameon=False)
            
        ##################################################
        #################### SAVE ########################
        ##################################################
        if SaveResults == "Yes":
            savefile_svg = SaveFileName+".svg"
            savefile_png = SaveFileName+".png"
            self.fig.savefig(savefile_svg,bbox_inches="tight")
            self.fig.savefig(savefile_png,bbox_inches="tight")
            print("Saved plots")
        elif SaveResults == "No":
            print("Plots not saved")
        else:
            print("Wrong ToSave input")

        self.draw()  # Redraw the canvas
        ##################################################

    def PlotResults_Conc(self, 
                     LEDwavelength,
                     timestamps, 
                     conc_opt, 
                     conc_exp, 
                     wavelengths, 
                     epsilon_R, epsilon_P, 
                     QY_AB_opt, QY_BA_opt, 
                     error_QY_AB, error_QY_BA,
                     ODEMethod,
                     SaveResults = "No",
                     SaveFileName = ""):
        
        ######################
        self.fig.clear()  # Clear the entire figure

        if SaveResults == "Yes": # re-define fig for savefig
            self.fig = Figure(figsize=(8, 4),dpi=600,constrained_layout=True)
            gs = gridspec.GridSpec(4, 6, figure=self.fig)
        else:
            gs = gridspec.GridSpec(4, 6, figure=self.fig)
            self.fig.set_constrained_layout(True)
            self.fig.suptitle(f'{LEDwavelength} nm irradiation: {ODEMethod} Method')
        
        axconc_values = self.fig.add_subplot(gs[0:2, 0:3])
        axconc_res_abs = self.fig.add_subplot(gs[2:4, 0:3])
        axAbs_spectra = self.fig.add_subplot(gs[0:2, 3:6])
        axAbs_res_abs = self.fig.add_subplot(gs[2:4, 3:6])
        axAbs_res_colourbar = self.fig.add_subplot(gs[2:4, 5])
        axAbs_res_colourbar.set_axis_off()
        
        ##################################################
        ############### CONCENTRATIONS ###################
        ##################################################        
        
        axconc_values.set_title("Concentration over time")
        
        ### Residuals ###
        residuals = conc_exp - conc_opt
        
        ### Top: data vs fit ###
        axconc_values.plot(timestamps, conc_exp[:,0], 'o', label='Reactant (exp.)', 
                           color=Defaults.colours_plot_first, alpha=0.4)
        axconc_values.plot(timestamps, conc_exp[:,1], 'o', label='Product (exp.)',
                           color=Defaults.colours_plot_last, alpha=0.4)
        axconc_values.plot(timestamps, conc_opt[:,0], '-', label='Reactant (fit)',
                           color=Defaults.colours_plot_first)
        axconc_values.plot(timestamps, conc_opt[:,1], '-', label='Product (fit)',
                           color=Defaults.colours_plot_last)
        axconc_values.set_ylabel("Concentration (mol/L)")
        axconc_values.legend()
        
        ### Bottom: residuals ###
        axconc_res_abs.plot(timestamps, residuals[:,0], 'o-', label='Reactant',
                            color=Defaults.colours_plot_first)
        axconc_res_abs.plot(timestamps, residuals[:,1], 'o-', label='Product',
                            color=Defaults.colours_plot_last,)
        axconc_res_abs.axhline(0, color='k', lw=1)
        axconc_res_abs.set_xlabel("Time (s)")
        axconc_res_abs.set_ylabel("Residuals (Exp - Fit)")
        axconc_res_abs.legend()
        
        ##################################################
        ################## ABSORBANCE ####################
        ##################################################

        
        axAbs_spectra.set_title(r"Absorption spectra (conc*$\epsilon$)")
        
        ## Absorbance ##
        A_exp = np.array([c[0] * epsilon_R + c[1] * epsilon_P for c in conc_exp])
        A_fit = np.array([c[0] * epsilon_R + c[1] * epsilon_P for c in conc_opt])
        A_res = A_exp - A_fit

        for i, (Ae, Af) in enumerate(zip(A_exp, A_fit)):
            axAbs_spectra.plot(wavelengths, Ae, 'o', alpha=0.3, 
                               label=r'Experimental' if i==0 else "")
            axAbs_spectra.plot(wavelengths, Af, '-', alpha=0.6, 
                               label=r"Fit:"+"\n"+
                               r"$QY_{fwd}$"+f": {QY_AB_opt}% "+u"\u00B1 "+f"{error_QY_AB}%"+"\n"+
                               r"$QY_{bwd}$"+f": {QY_BA_opt}% "+u"\u00B1 "+f"{error_QY_BA}%" if i==0 else "")
        
        # axAbs_spectra.set_xlabel("Wavelength (nm)")
        axAbs_spectra.set_ylabel("Absorbance")
        axAbs_spectra.legend()
        
        #################### RESIDUALS HEATMAP ###################
        im1 = axAbs_res_abs.imshow(A_res,
                        extent=[wavelengths.min(), wavelengths.max(),
                                timestamps[-1], timestamps[0]],
                        aspect='auto', cmap='bwr',
                        vmin=-np.max(np.abs(A_res)), vmax=np.max(np.abs(A_res)))
        axAbs_res_abs.set_xlabel("Wavelength (nm)")
        axAbs_res_abs.set_ylabel("Time (s)")
        self.fig.colorbar(im1, ax=axAbs_res_colourbar, label="Abs Residuals (Exp - Fit)")
        
        ##################################################
        #################### SAVE ########################
        ##################################################
        if SaveResults == "Yes":
            savefile_svg = SaveFileName+".svg"
            savefile_png = SaveFileName+".png"
            self.fig.savefig(savefile_svg,bbox_inches="tight")
            self.fig.savefig(savefile_png,bbox_inches="tight")
            print("Saved plots")
        elif SaveResults == "No":
            print("Plots not saved")
        else:
            print("Wrong ToSave input")

        self.draw()  # Redraw the canvas
        ##################################################


##################################################
    
