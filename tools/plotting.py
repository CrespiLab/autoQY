import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
import numpy as np

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

        # Ensure line_positions is initialized for the given idx
        if self.idx is not None:
            if len(self.parent_window.line_positions) <= self.idx:
                # Initialize with default positions if idx is out of range
                self.parent_window.line_positions.append([0] * 12)
            
            line_positions = self.parent_window.line_positions[self.idx]
            
        else:
            line_positions = self.line_positions

        # Create new boxes based on updated line positions
        for i in range(0, len(line_positions) - 1, 2):
            box = self.ax.axvspan(line_positions[i], 
                                line_positions[i+1], 
                                color=colours[int(i/2)], alpha=0.1)
            self.boxes.append(box)

        # Redraw the canvas to reflect changes
        self.draw()



    def on_press(self, event):
        if event.inaxes != self.ax:
            return

        # Check if the mouse is near any line
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
        #print(self.parent_window.line_positions)
        # Update the corresponding label and line position variable with the new x position
        line_index = self.lines.index(self.selected_line)
        self.parent_window.update_display(self.idx, line_index, new_x)  # Pass self.idx here
        self.parent_window.line_positions[self.idx][line_index] = new_x  # Update using self.idx
        self.line_positions[line_index] = new_x
        #print(self.parent_window.line_positions)
        #print('test')
        ## Update the boxes
        self.update_boxes()


    def on_release(self, event):
        self.selected_line = None
        self.press_event = None


#    def plot_sections(self, x, RefPower, sections, filename):
#        """
#        Plot the power measurements with highlighted sections.
#        """
#        self.ax.clear()  # Clear previous plot
#        self.ax.set_title(filename)
#        self.ax.plot(x, RefPower, label="Power", color="black")
#
#        # Highlighting sections
#        self.ax.axvspan(sections["start_0"], sections["end_0"], alpha=0.1, color='red')
#        self.ax.axvspan(sections["start_1"], sections["end_1"], alpha=0.1, color='blue')
#        self.ax.axvspan(sections["start_0_1"], sections["end_0_1"], alpha=0.1, color='red')
#        self.ax.axvspan(sections["start_2"], sections["end_2"], alpha=0.1, color='red')
#        self.ax.axvspan(sections["start_3"], sections["end_3"], alpha=0.1, color='blue')
#        self.ax.axvspan(sections["start_4"], sections["end_4"], alpha=0.1, color='red')
#
#        self.ax.set_ylabel("Power (mW)")
#        self.ax.set_xlabel("Index")
#        self.ax.legend()
#
#        # Add vertical lines to mark the sections
#        self.lines = []  # Reset the lines list
#        for key in ["start_0", "end_0", "start_1", "end_1", "start_0_1", "end_0_1", "start_2", "end_2", "start_3", "end_3", "start_4", "end_4"]:
#            line = self.ax.axvline(x=sections[key], linestyle='--', color='green')
#            self.lines.append(line)  # Store the lines
#        print(f"Generated {len(self.lines)} lines for section markers.")
#        self.draw()  # Render the plot
#
#
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

                
    def Plot_LEDemission_Processed(self, absorbance, wavelengths,
                                em_wl, em_int,
                                index_first, index_last,
                                em_int_proc):
        self.fig.clear()  # Clear the entire figure
        
        # Create a grid for the subplots
        gs = plt.GridSpec(2, 1, figure=self.fig, height_ratios=[1, 1], hspace=0.4)

        # Add subplots using self.fig.add_subplot and GridSpec
        ax1 = self.fig.add_subplot(gs[0])

        # print(f"absorbance.shape: {absorbance.shape}")


        for i in range(0,absorbance.shape[1]):
            ax1.plot(wavelengths, absorbance.T[i])
            # print(i)

        # Set common x-axis properties for both subplots
        ax1.set_xlabel("Wavelength (nm)")
        ax1.set_xlim(220, 650)

        # Customize the first subplot (Spectral Data)
        ax1.set_title("Spectral Data")
        ax1.set_ylabel("Absorbance")
        ax1.set_ylim(-0.05, 2.5)

        ax2 = self.fig.add_subplot(gs[1])
        ax2.set_xlabel("Wavelength (nm)")
        ax2.set_xlim(220, 650)

        # Plot LED emission data in the second subplot (ax2)
        ax2.plot(em_wl, em_int, label="Untreated (in this .py file at least)")
        ax2.plot(em_wl[index_first:index_last], em_int_proc[index_first:index_last],
                label="Smoothed, removed zeroes\nand applied threshold")
        ax2.legend(fontsize=8)
        ax2.set_title("LED Emission")
        ax2.set_ylabel("Intensity")

        # Manually adjust figure layout to remove excess padding
        self.fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)  # Adjust margins as needed
        self.fig.set_size_inches(10, 8)  # Adjust the size if needed

        # Automatically adjust subplots to prevent overlap
        self.fig.tight_layout()

        # Display the figure
        # plt.show()
        self.draw()  # Redraw the canvas


    def plot_Epsilons(self, wavelengths, e_A_inter, e_B_inter, emission_inter):
        """Plot interpolated epsilons and LED emission."""
        self.fig.clear()  # Clear the entire figure
        
        # Use GridSpec to create a grid layout
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.5)  # Two rows, evenly spaced, with some space between

        # Create first plot (Epsilons) in the first row
        ax1 = self.fig.add_subplot(gs[0])
        
        # Plot epsilons on the first subplot
        ax1.plot(wavelengths, e_A_inter, label="Reactant")
        ax1.plot(wavelengths, e_B_inter, label="Product")
        ax1.legend(fontsize=12)
        ax1.set_xlabel("Wavelength (nm)")
        ax1.set_xlim(220, 650)
        ax1.set_title("Epsilons, interpolated")
        ax1.set_ylabel(r"$\epsilon$ (M$^{-1}$ cm$^{-1}$)")

        # Create second plot (LED Emission) in the second row
        ax2 = self.fig.add_subplot(gs[1])
        ax2.plot(wavelengths, emission_inter)
        ax2.set_xlabel("Wavelength (nm)")
        ax2.set_xlim(220, 650)
        ax2.set_title("LED emission, interpolated")
        ax2.set_ylabel("Intensity")

        self.draw()  # Redraw the canvas

    def PlotResults(self, LEDwavelength,
                    timestamps, 
                    conc_opt,
                    absorbance,
                    index,
                    total_abs_fit,
                    residuals,
                    QY_AB_opt, QY_BA_opt, 
                    error_QY_AB, error_QY_BA,
                    CalculationMethod,
                    SaveResults = "No",
                    SaveFileName = ""):
        
        ######################
        self.fig.clear()  # Clear the entire figure

        if SaveResults == "Yes": # re-define fig for savefig
            self.fig = Figure(figsize=(8, 4),dpi=600,constrained_layout=True)

        gs = gridspec.GridSpec(4, 2, figure=self.fig)
        
        # self.fig.suptitle(f'{LEDwavelength} nm irradiation: Integration Method')
        self.fig.suptitle(f'{LEDwavelength} nm irradiation: {CalculationMethod} Method')
        
        ##################################################
        axresults_conc = self.fig.add_subplot(gs[0:3, 0])
        axresults_Abs = self.fig.add_subplot(gs[0:3, 1])
        axresults_res = self.fig.add_subplot(gs[3, 1])
        axresults_notes = self.fig.add_subplot(gs[3, 0])
        
        ##################################################
        ############### CONCENTRATIONS ###################
        ##################################################
        axresults_conc.set_title("Concentrations over time")
        
        axresults_conc.plot(timestamps, conc_opt[:,0], color = self.colours[2], label = "Reactant")
        axresults_conc.plot(timestamps, conc_opt[:,1], color = self.colours[4], label = "Product")
        
        axresults_conc.set_ylabel("Concentration (M)")
        # axresults_conc.yaxis.set_minor_locator(AutoMinorLocator(2))
        axresults_conc.set_xlabel("Time (s)")
        
        # axresults_conc.legend()
        
        ##################################################
        ###################### NOTES #####################
        axresults_notes.axis("off")
        # axresults_notes.text(0, 0, f"Starting Percentage A: {StartPercentage_A} %")
        ##################################################
        
        #### Plot the experimental data and the fitted total absorbance curve
        axresults_Abs.plot(timestamps, absorbance[index,:], linestyle='-', 
                    color=self.colours[4], linewidth=8, alpha=0.5, label='Experimental Data')
        axresults_Abs.plot(timestamps, total_abs_fit[:,index], linestyle='--', color=self.colours[2], 
                    label=f"Fit:\nQY_RtoP: {QY_AB_opt:.3f} "+u"\u00B1 "+f"{error_QY_AB:.3f}\
                    \nQY_PtoR: {QY_BA_opt:.3f} "+u"\u00B1 "+f"{error_QY_BA:.3f}")
        axresults_Abs.set_title('Fit vs Experimental')
        # axresults_Abs.set_xlabel('Time (s)')
        axresults_Abs.set_xticklabels([]) ## put this before set_xlim, otherwise it resets the xlim
        
        axresults_Abs.set_ylabel(r'Absorbance at $\lambda_{irr}$')
        # axresults_Abs.legend()
        
        #################### RESIDUALS ###################
        axresults_res.plot(timestamps, residuals[index,:], color=self.colours[4], label='Residuals')
        axresults_res.set_xlabel('Time (s)')
        axresults_res.set_ylabel('Residuals Abs')
        ##################################################
        
        for i in [axresults_conc, axresults_Abs, axresults_res]:
            i.xaxis.set_minor_locator(AutoMinorLocator(2))
            i.yaxis.set_minor_locator(AutoMinorLocator(2))

        for i in [axresults_conc, axresults_Abs]:
            i.legend(frameon=False)
            
        ##################################################
        #################### SAVE ########################
        ##################################################
        if SaveResults == "Yes":
            savefilename = SaveFileName+f"_{CalculationMethod}"
            self.fig.savefig(savefilename+".svg",bbox_inches="tight")
            self.fig.savefig(savefilename+".png",bbox_inches="tight")
            print("Saved plots")
        elif SaveResults == "No":
            print("Plots not saved")
        else:
            print("Wrong ToSave input")

        self.draw()  # Redraw the canvas
        ##################################################

##################################################
    