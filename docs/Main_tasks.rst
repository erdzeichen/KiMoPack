Main Tasks overview
====================

This menu is a shortcut to the main function and tasks that are performed during an analysis.
In general one opens one or multiple data Files and after defining a number of shaping parameter, that 
e.g. set the axis limits or correct the arrival time of different wavelength,  plots various graphs.
Different measurements or Fits can be compared and the results saved in various forms.

* :ref:`Opening of data`

	* Open single file: 			:meth:`self.TA.__init__`
	* Open many files: 				:meth:`plot_func.GUI_open`
	* Combine many scans:			:meth:`plot_func.Summarize_scans`

* :ref:`Shaping of Data`

	* Background correction:			:meth:`plot_func.TA.Background`
	* Filter bad data:				:meth:`plot_func.TA.Filter_data`
	* Correct arrival time (Chirp)	:meth:`plot_func.TA.Cor_Chirp` 
	* :ref:`Data shaping settings that affect the fits`
	* :ref:`Plot shaping options without influence on the fitting`

* :ref:`Plotting functions`

	* Plotting non Fitted data:	:meth:`plot_func.TA.Plot_RAW`
	* Plotting Fitted data:		:meth:`plot_func.TA.Plot_fit_output`
	* Adjust fonts in plots:		:meth:`plot_func.changefonts`

* :ref:`Fitting, Parameter optimization and Error estimation`

	* Fitting data:				:meth:`plot_func.TA.Fit_Global`

* :ref:`Comparative plotting`

	* Compare spectra:				:meth:`plot_func.TA.Compare_at_time`
	* Compare kinetics:				:meth:`plot_func.TA.Compare_at_wave`
	* Compare calculated spectra (SAS or DAS):	:meth:`plot_func.TA.Compare_DAC`

* :ref:`Data Export and Project Saving`

	* Copy project					:meth:`plot_func.TA.Copy`
	* Save Project as hdf5			:meth:`plot_func.TA.Save_project`
	* Save Plots					:func:`plot_func.TA.Save_Plots`
	* Save Plots as Powerpoint		:func:`plot_func.TA.Save_Powerpoint`
	* Save/export data as ascii/text	:meth:`plot_func.TA.Save_data`

.. figure:: _static\\structure.png

