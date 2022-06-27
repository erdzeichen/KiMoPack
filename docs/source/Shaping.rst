Shaping of Data
===============

In the Following sections we discuss the parameter and values that are used to
filter and shape measured data. In general all loaded data read is
stored in the un-altered "ta.ori_ds". A second matrix for the same data
is created named ta.ds.

Only the function "Filter_data" works on the dataset ds_ori and ds.
The chirp correction creates a new ds from ds_ori.
The background correction is applied to ta.ds (unless a
specific matrix is given´). All other parameter are only applied during a
function and do not alter ta.ds. That means that in each of the
plotting/fitting functions a local copy of the ta.ds is created (using
the function sub_ds) to which all the shaping is applied.

**The intended work flow is:** 

	#. loading of data :ref:`Opening of data`
	#. Filtering of data :ref:`Bad data Filter`
	#. (optional) chirp correction :ref:`Arrival time correction`
	#. (optional) Background correction :ref:`Background subtraction`
	#. setting of parameter for fit :ref:`Data shaping settings that affect the fits`
	#. setting of parameter for plot :ref:`Plot shaping options without influence on the fitting`
	#. all the plotting/fitting. :ref:`Plotting functions` and :ref:`Fitting, Parameter optimization and Error estimation`
	#. saving/exporting :ref:`Data Export and Project Saving`

The point 5 and 6 (Parameter) can be easily changed many times and a new plot/fit 
generated. Important to not is that the parameter are stored with the object.
This means that a parameter that is explicitly set, will stay until it is 
overwritten or the object is fresh loaded. So if e.g. commenting out a certain
parameter does not return the value to its Default "None". The Parameter needs 
to be set explicitly to::

	ta.intensity_range = 3e-3 
	#ta.intensity_range = 3e-3 #no effect
	ta.intensity_range = None

Often it is faster to reload the object by choosing "run all above" in the 
Notebook.

Bad data Filter
---------------

Filter bad data:				:meth:`self.Filter_data()<plot_func.TA.Filter_data>`

In some cases there are bad data points or other strange things. NA
values will normally be replaced by 0 during import and all data is
converted into floats during import. In many recording software
(including Pascher instruments) a specific **value** is used to indicate
that something went wrong. This function filters everything bigger than
this value as error. Real "NaN" values are filtered during the import of Data.
There is the option to either drop the times that contain bad values or to replace
bad values with a specific value. There is the option to put a uppervalue, lowervalue
or a single value that is then used for upper and (as negative) for the lower value.

If the filtering does not work, a manual way of filtering is 
ta.ds (the chirp corrected data) or ta.ds_ori[ta.ds_ori>20]=0 is the classical way to filter 

Arrival time correction
-----------------------

Correct arrival time (Chirp)	:meth:`self.Cor_Chirp()<plot_func.TA.Cor_Chirp>` 

*Cor_Chirp* is a powerful Function to correct for a different arrival times of
different wavelength (sometimes call chirp). 

In general if a file is opened for the first time this function is opening 
a plot and allows the user to select a number of points, which are then 
approximated with a 4th order polynomial and finally to select a point 
that is declared as time zero. The observed window as well as the intensities 
and the colour map can be chosen to enable a good correction. Here a fast 
iterating colour scheme such as "prism" is often a good choice. In all of the 
selections a left click selects, a right click removes the last point and 
a middle click (sometime appreviated by clicking left and right together) 
finishes the selection. If no middle click exists, the process
automatically ends after max_points (40 preset).

Note that scattercut, bordercut and intensity_range can be used to adjust intensity.

After the first run the polynom is stored in self.fitcoeff, a new matrix 
calculated from self.ds_ori that is stored as self.ds and a file stored in the 
same location as the original data. The second time the function *Cor_Chirp* is 
run the function will find the file and apply the chirp correction automatically.

If one does want to re-run the chirp correction the function *Man_Chirp* does
not look for this file, but creates after finishing a new file.

Alternatively the polynom or a filename can be given that load a chirp correction
(e.g. from a different run with the same sample).
The function *Cor_Chirp* selects in the order: 

	#. "fitcoeff"
	#. "other files"
	#. "stored_file"
	#. call Man_Chirp (clicking by hand)

Correct arrival time (Chirp)	:meth:`self.Cor_Chirp()<plot_func.TA.Cor_Chirp>` 
Manual overwrite arrival time correction 	:meth:`self.Man_Chirp()<plot_func.TA.Man_Chirp>` 

Background subtraction
----------------------

Background correction:			:meth:`self.Background()<plot_func.TA.Background>`

This tool is one of two ways to remove a flat background from the data (typically seen before t=0). 
This tool averages for each measured  wavelength separately the values from 'lowlimit' to 'uplimit' and 
subtracts it from the data. The low and uplimit can be set 
anywhere to substract any background. (so one could e.g. substract a product 
instead) It is important to note that many problems during measurements might
be visible in the data before time zero. So I recommend to first
plot without background correction and only after this inspection 
apply the background correction. 
The fit function has its own way to calculcate and apply a background 
That could be used instead (but making the fit less stable) 

Data shaping settings that affect the fits
------------------------------------------

in general the data is handled in each of the plotting/fitting functions
separately. In each function a copy of the matrix with the limitation
below is created. 
A number of parameter cut and potentially rebin the raw measured data and as such affet the fit. 
The typical workflow would therefore be to adjust these parameter before the fitting stage using the 
RAW plotted fits as a feedback.

* Cut the outside limits of the spectrum: *Bordercut*
* Blank one or multiple regions in the spectrum (e.g. suppress scatter) *Scattercut*
* Cut the outside of the time axis: *timelimits*
* Blank one or multiple temporal regions (e.g. around t=0) *ignore_time_region*
* rebin the temporal axis (useful for e.g. steady state long term UV-vis data) *time_bin*
* rebin the spectral axis (useful for prism based spectrometer) *wave_nm_bin* 

For further details and examples see: :meth:`self.__make_standard_parameter()<plot_func.TA.__make_standard_parameter>`
or e.g. the general plotting function :meth:`pf.plot_raw()<plot_func.plot_raw>`.

The parameter that only change the plots are discussed in :ref:`Plot shaping options without influence on the fitting`
