Plotting functions
==================

	* Plotting non Fitted data:	:meth:`self.Plot_RAW()<plot_func.TA.Plot_RAW>`
	* Plotting Fitted data:		:meth:`self.Plot_fit_output()<plot_func.TA.Plot_fit_output>`
	* Interative Plotting RAW and Fitted:	:meth:`self.Plot_Interactive()<plot_func.TA.Plot_Interactive>`
	* Adjust fonts in plots:		:meth:`pf.changefonts()<plot_func.changefonts>`

One core function of this tool is to create plots for evaluation and
publication. Internally there are a number of separate functions that 
create each plot type (see below). The methods Plot_RAW and Plot_fit_output 
wrap the parameter into the object and simplify their use. Two additional functions
provide additional features. Both "Save_Plots" and "Save_Powerpoint" are 
calling both plot functions and dump their output into separate figure files or two 
slides of a power point file.

Common to both plotting function is that either a single plot can be called by giving
the plotting parameter or a series of plots (Default) by giving a list of number with 
e.g. "range(3)".

Most of the plot defining parameter (like, what for which wavelength the kinetic 
is plotted or at what times the kinetics are extracted are defined by the
:ref:´Plot shaping options without influence on the fitting´.

Plot_RAW
--------

:meth:`self.Plot_RAW()<plot_func.TA.Plot_RAW>` plots all raw figures. The different figures can be called
separately or with a list of plots (standard) e.g. plotting=range(4)
call plots 0-3, plotting=1 a single plot. The plots have the following
numbers: 0 - Matrix, 1 - Kinetics, 2 - Spectra, 3 - SVD. The plotting
can take all parameter from the "ta" object. See:
:meth:`self.Plot_RAW()<plot_func.TA.Plot_RAW>`

Plot_fit_output
---------------

:meth:`self.Plot_fit_output()<plot_func.TA.Plot_fit_output>` plots the fit results. For this is uses the data
contained in the shaped and cut datasets that were used for the fit,
including all rebinning or temporal restrictions. The figures can be
called separately or with a list of plots (standard)
The plotting function takes all parameter from the object.

	:meth:`self.Plot_fit_output()<plot_func.TA.Plot_fit_output>`

**Contents of the plots**

   #. DAC contains the assigned spectra for each component of the fit. For
      a modelling with independent exponential decays this corresponds to
      the "Decay Associated Spectra" (DAS). For all other models this
      contains the "Species Associated Spectra" (SAS). According to the
      model the separate spectra are labeled by time (process) or name, if
      a name is associated in the fitting model.
   #. summed intensity. All wavelength of the spectral axis are summed for
      data and fit. 
   #. plot kinetics for selected wavelength 
   #. plot spectra at selected times
   #. plots matrix (measured, modelled and error Matrix). The parameter are
      the same as used for the corresponding RAW plot with the addition of
      "error_matrix_amplification" which is a scaling factor multiplied
      onto the error matrix.
   #. concentrations. In the progress of the modelling/fitting a matrix is
      generated that contains the relative concentrations of the species
      as function of time. 

This function is a convenience function and is suppose to be used in
conjunction with the object and the embedded parameter (see above). The
use of qt as backend allows the easy customization of the plots via the
GUI. If the plots are saved as "svg" they can easily be adjusted in
inkscape or similar afterwards.
For more details see: :meth:`self.Plot_fit_output()<plot_func.TA.Plot_fit_output>`

Plot shaping options without influence on the fitting
-----------------------------------------------------

In addition to the general shaping parameter from section :ref:`Data shaping settings that affect the fits`
a number of parameter only affect one or multiple of the plots but not the fitting of the data.

* 	The plotting of the kinetics is governed by the selection of the wavelength in the list **rel_wave** 
	and the width of each **wavelength_bin**
* 	The plotting of the spectra is governed by the selection of the timepoint in the list  **rel_time** 
	and potentially a percentual binning around this time-point with **time_width_percent**. If this is set to 0
	then the measured timepoint is used. 
*	The intensity (color) in the 2 plots as well as the height of the y-axis is determined by the **intensity_range** 
	parameter that can be set symmetric or a-symmetric for best representation. With **log_scale** 
	This intensity can be scaled logarithmic and **error_matrix_amplification** only amplifies the intensity of the 
	difference matrix (measured-fitted) in the 2d plots
* 	The color scheme can be set very flexible using the Matplotlib palets, or a manually provided color scheme 
	(e.g. university colors)
*	The titles of all plots are chosen either by the filename or can be given flexible in each plotting functions 
	through the title parameter. All the plots can be automatically saved if **save_figures_to_folder** is set to True,
	Which is useful for fast surveys, otherwise the method :meth:`self.Save_Plots()<plot_func.TA.Save_Plots>` 
	stores all plots (see :ref:`Data Export and Project Saving`). The axis labels are accessible via the **baseunit** 
	and the Fonts are accessible via the function :meth:`pf.changefonts()<plot_func.changefonts>`
*	The parameter **equal_energy_bin** can be set to a value which results in that the spectral plots are shown in enqual energy 
	bins. This is useful for tracking vibrations and such. As of version 6.7.1 this is only happing for the RAW plotting. 

interactive Plotting
---------------------

Interactive plot function that allows the interactive slicing of both time and wavelength. The main parameter of the object apply



extended Raw plotting
---------------------

:meth:`self.Plot_raw()<plot_func.Plot_raw>` is an extended function. All the parameters are 
accessible (and need then to be set manually). This function also plots a single 
or multiple plots bzt setting the "plotting" parameter. 

There are even more detailed manipulations possible by using the
separate plot functions:
 
	* for plotting kinetics at fixed wavelength: :func:`pf.plot1d()<plot_func.plot1d>`
	* for plotting spectra at selected times :func:`pf.plot_time()<plot_func.plot_time>` 
	* for plotting the data matrix :func:`pf.plot2d()<plot_func.plot2d>`
	* for plotting the 3 fit data matrix :func:`pf.plot2d_fit()<plot_func.plot2d_fit>`
	* for the SVD plots. :func:`pf.SVD()<plot_func.SVD>` 
	
Each of the functions allows to hand in an axis and thus plot multiple things
