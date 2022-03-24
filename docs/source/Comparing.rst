Comparative plotting
====================

The comparative plotting is an important tool to compare different
measurements (e.g. different conditions), different fits or steady state 
spectra. In general you can compare different kinetics 
(at one or multiple fixed wavelength) with :ref:`Compare_at_wave` , 
compare different spectra at one or multiple given timepoints with
:ref:`Compare_at_time` and compare the extracted spectra (decay associated 
or species associated) with :ref:`Compare_DAC`. The essential idea is that the
parameters in the project that contains the comparative plotting are used
for all the plots. So the ta.bordercut or ta.intensity is used for all
the plot, independent of e.g. the intensity in the other projects.

New is that the compare functions take "other" as a parameter, which can
be either a single or multiple projects (TA - objects). These projects
need to be loaded into the program. Loading a project can be done by 
having them open from prior import or analysis (e.g. when comparing different fits)
and then using 

:meth:`self.Copy()<plot_func.TA.Copy>` 

More usual other (saved) projects will be
opened with the function 

:meth:`pf.GUI_open()<plot_func.GUI_open>` 

. See :ref:`Opening multiple files`
for more information on that.

As this comparision very often contains a lot of files the images are
automatically saved using the filenames and the wavelength/time points.
The images are however open and if the standard approach using QT was
used can be manipulated using the GUI tools. So is e.g. the conversion
into a linear time-axis currently not implemented, but can easily be
achieved by changing the axis in the QT GUI.

A very important function provided by this set of tools is the comparision against 
other spectra. So can for examples be reference spectra (e.g. UV-vis) be added to 
the plots.

Normalization and Scaling
----------------------------------------

An important option is the normalization in a certain window that
applies for both Compare_at_time and Compare_at_wave. Very often data
needs to be normalized before it can be compared to anything e.g. to
the size of the ground state bleach or an excited state absorption.
Here I offer the normalization in a certain window in time and space.
In this window a value in the "ta" and then each of the "other"
projects is calculated. The intensity of each in the other projects
(but not of the "ta" project) is then mutliplied by the quotient of
this value in this specific window. This means e.g. that even if the
kinetics is plotted for 600nm the normalization can be performed at
1-2ps at 450-500nm. This is very useful to plot e.g. the efficiency of
product formation in the study of catalytic processes. For this
normalization a "window" needs to be defined in the order::

 [Start time, End time, Start wavelength, End Wavelength]

Care should be take to where this normalization is performed. The region
just around t=0 is dangerous due to the artifacts, that do not represent
real values. If external values are suppose to be used for scaling, the
individual intensities can be manipulated. For each of the loaded
projects ta.ds is the matrix that contains the data. With::

	"ta.ds*=1.1"

could for example the individual intensity be raised by 10%. But be aware 
that with this action you are changing the data in the object. The original 
data ta.ds_ori remains unchanged. If you save and reload the data, the intensity
will revert to the originally measured value.

Compare_at_time
-------------------

This function plots multiple spectra into the same figure at one or multiple given timepoints (rel_time) and 
allows for :ref:`Normalization and Scaling` 

Very useful to compare the spectra for different solvents or quenchers, or e.g. different fits. 
The ta.time_width_percent parameter defines if this is a single time 
(if time_width_percent = 0) or an integrated window.

A normalization window can be given at which all the plotted curves are normalized to. 
This window does not have to be in the plotted region. See :ref:`Normalization and Scaling`
		
Very often one would like to compare the measured spectra at a certain
time to an external spectrum (e.g. spectro-electro-chemistry or steady
state absorption). This can be done by loading a specific spectrum into
a DataFrame and handing this data Frame to the comparision function. The
function can also be used to plot e.g. the measured spectra vs. an
external spectrum without giving any "other" Projects. 

For more information, details on the parameter and examples see:

:meth:`self.Compare_at_time()<plot_func.TA.Compare_at_time>`

Compare_at_wave
--------------------

This function plots multiple kinetics into the same figure at one or
multiple given wavelength (rel_wave) and  allows for :ref:`Normalization and Scaling`  
Very useful to compare the kinetics for different quencher concentrations or
pump powers, or e.g. different fits. The parameter width or the general ta.wavelength_bin 
defines the width of the spectral window that is integrated and shown.

A normalization window can be given at which all the plotted curves are normalized to. 
This window does not have to be in the plotted region. See :ref:`Normalization and Scaling`

Often multiple wavelength are to be plotted, and if at the same time
many projects go into the same plot, things tend to get messy. As the
files are saved separately this approach proofed to be useful.

For more information, details on the parameter and examples see:

:meth:`self.Compare_at_wave()<plot_func.TA.Compare_at_wave>`

Compare_DAC
----------------

This is a convenience function to plot multiple extracted spectra (DAS
or species associated) into the same figure or into a separate figure
each. Other should be ta.plot_func objects (loaded or copied). By
standard it plots all into the same window. If all project have the same
number of components one can activate "separate_plots" and have each
separated (in the order created in the projects).

The "Spectra" parameter allows as before the inclusion of an external
spectrum. Others is optional and I use this function often to compare
species associated spectra with one or multiple steady state spectra.

For more information, details on the parameter and examples see:

:meth:`self.Compare_DAC()<plot_func.TA.Compare_DAC>`
