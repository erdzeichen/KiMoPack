Data Export and Project Saving
==============================

Save_Plots
-------------------

Save Plots						:meth:`self.Save_Plots()<plot_func.TA.Save_Plots>`

Convenience function that calls both "Plot_RAW" and if possible
"Plot_fit_output" and saves the plots. Effectively this is intented to
be used at the end of a quick exploratory work to capture a status. The
parameter are indentical to the two plotting functions and are handed
through.

Save_Powerpoint
--------------------

Save Plots as Powerpoint		:meth:`self.Save_Powerpoint()<plot_func.TA.Save_Powerpoint>`

Convenience function that calls both "Plot_RAW" and if possible
"Plot_fit_output" and saves the plots as "png". 

Then it creates a power point file with one slide for the RAW plots and 
one slide for the Fits.
Effectively this is intented to be used at the end of the a quick
exploratory work to capture a status and create a quick presention
slide. The parameter are intentical to the plotting functions and are
handed through. The additional switches save_RAW and save_Fit are
convenience to allow for faster processing.

If the "savetype" contains 'png', 'svg', or 'pdf' then a summary file is created 
that looks close to the powerpoint file.

Saving of the  project
--------------------------------------

Save Project as hdf5			:meth:`self.Save_project()<plot_func.TA.Save_project>`

This software allows the saving of the TA project as a HDF5 file that
can be reloaded. The HDF5 file contains all the set parameter as well as
the fit results (if applicable) and the raw data. To reduce the space
consumption we save ta.ds_ori and the parameter that are used to create
ta.ds. If manual changes were made to ta.ds, these have to be stored
externally. As only the obvious errors are filtered in ta.Filter_data
this is can savely replace the original data File. We are also saving
the arrival-time (chirp) correction in the file and restore the chirp
corrected data ta.ds during the import. The import function understands
the file type and re-creates the object.

The one limitation to this method is the external fit function. If an
external **ta.mod** is used, the save function stores the name and the
documentation string of this function as a string. So after reloading of
the analysis object the external function will have to be set with
ta.mod=imported_function. The parameter of the fit are however stored.
Only the filename and the path of the file can be changed during saving
of the project. If left empty the path and filename of the original
ASCII file is used.

Save ASCII data
---------------------------

Save/export data as ascii/text	:meth:`self.Save_data()<plot_func.TA.Save_data>`

This is a convenient function to export the data for use with other
plotting programs, the chirp corrected data, all the slices defined by
ta.rel_wave and ta.rel_time for both the fits and a the RAW data. The
external options include: save_RAW,save_Fit while there is an
automatic that recognizes if for example fit data is present, this
switch allows the manual selection which datasets are stored.
save_slices selected if the slices defined by ta.rel_wave and
ta.rel_time are saved save_binned this switch chooses if the chirp
corrected and rebinned dataset (ta.ds with ta.wavelength_nm_bin) is
saved. If the ta.wavelength_nm_bin is None, this saves the chirp
corrected RAW data. filename sets the basis filename that is used for
all the files path this can be a full path or a simple string,
defining a folder relative to the folder in ta.path. If the folder
does not exist, it will be created, if it exists a file with exactly
the same name will be overwritten without confirmation sep defines the
separator user to separate different values. Standard is a "TAP". A good 
choice would also be a space or a comma, unless you are
located in one of the countries that uses commas for decimal points.
Decimals will be separated with "dots".

This function by default also dumps a text file with the fit results
