Main Tasks overview
====================

This menu is a shortcut to the main function and tasks that are performed during an analysis.
In general one opens one or multiple data Files and after defining a number of shaping parameter, that 
e.g. set the axis limits or correct the arrival time of different wavelength,  plots various graphs.
Different measurements or Fits can be compared and the results saved in various forms.

* :ref:`Opening of data`

	* Open single file: 			:meth:`pf.TA()<plot_func.TA.__init__>`
	* Open many files: 				:meth:`pf.GUI_open()<plot_func.GUI_open>`
	* Combine many scans:			:meth:`pf.Summarize_scans()<plot_func.Summarize_scans>`

* :ref:`Shaping of Data`

	* Background correction:			:meth:`self.Background()<plot_func.TA.Background>`
	* Filter bad data:				:meth:`self.Filter_data()<plot_func.TA.Filter_data>`
	* Correct arrival time (Chirp)	:meth:`self.Cor_Chirp()<plot_func.TA.Cor_Chirp>` 
	* :ref:`Data shaping settings that affect the fits`
	* :ref:`Plot shaping options without influence on the fitting`

* :ref:`Plotting functions`

	* Plotting non Fitted data:	:meth:`self.Plot_RAW()<plot_func.TA.Plot_RAW>`
	* Plotting Fitted data:		:meth:`self.Plot_fit_output()<plot_func.TA.Plot_fit_output>`
	* Adjust fonts in plots:		:meth:`pf.changefonts()<plot_func.changefonts>`

* :ref:`Fitting, Parameter optimization and Error estimation`

	* Fitting data:				:meth:`self.Fit_Global()<plot_func.TA.Fit_Global>`

* :ref:`Comparative plotting`

	* Compare spectra:				:meth:`self.Compare_at_time()<plot_func.TA.Compare_at_time>`
	* Compare kinetics:				:meth:`self.Compare_at_wave()<plot_func.TA.Compare_at_wave>`
	* Compare calculated spectra (SAS or DAS):	:meth:`self.Compare_DAC()<plot_func.TA.Compare_DAC>`

* :ref:`Data Export and Project Saving`

	* Copy project					:meth:`self.Copy()<plot_func.TA.Copy>`
	* Save Project as hdf5			:meth:`self.Save_project()<plot_func.TA.Save_project>`
	* Save Plots					:func:`self.Save_Plots()<plot_func.TA.Save_Plots>`
	* Save Plots as Powerpoint		:func:`self.Save_Powerpoint()<plot_func.TA.Save_Powerpoint>`
	* Save/export data as ascii/text	:meth:`self.Save_data()<plot_func.TA.Save_data>`

.. figure:: _static\\structure.png

