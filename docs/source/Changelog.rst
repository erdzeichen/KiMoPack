*********
Changelog
*********
In this changelog I will only mention the new capabilities. All other version jumps are bug fixes

.. _release-7.12.0:
7.12.0
========
Added the option to return the figures from the common plotfunction setting the option "return_figures_handles" will in Plot_RAW and Plot_Fit_output trigger that a dictionary of the plots is returned.
Added X-ray workflow
Added the general option pf.halfsize=False If set to true all figures will be half the size (better for inline plotting)

.. _release-7.11.0:
7.11.0
========
Add the option to use Cor chirp for only time shifting
Add Streak camera workflow

.. _release-7.10.0:
7.10.0
========
Added time delay for printing fit output to avoid overflow of notebook output, has big impact on speed for more complex model fits.


.. _release-7.7.0:
7.7.0
========
Added X-ray support
added novel global variable that permits half size plots (nicer in inline mode)


.. _release-7.5.0:
7.5.0
========
Changed the dimension of the dc in all the models. This sped up the speed by nearly a factor of 10. There has been so many smaller changes over the last versions that I go to 7.5

.. _release-7.4.0:
7.4.9
========
Cleaned up previous updates.



.. _release-7.3.0:
7.3.6
========
Introduced that Print_results can dump the result parameter into a file.


7.3.0
========

introduced: sub_sample, pulse_sample  into the fitting (not yet documented) but they allow to temporarily add time points to the model (that willl then be removed again) 
			Pulse_sample is needed when the pulse is not in the modelled data (e.g. when the time_limit or ignore_time_regions is set
			sub_sample divided the times and is useful if the measure data is to sparse in time.
			if the parameter "sub_steps" is present this will be used to define the number of sub_steps in the iterative sampling


7.2.17
========

The to sparsely measured datapoints can now be sub-sampled with ta.Fit_Global(sub_sample=10)
The intensity is now proper if the modelled points do not include the pump pulse


7.2.5
=======

On Windows the fit can now be interrupted if "q" is pressed

7.2.1
=======

Add write_paras as Fit_Global option. This will print the params continuously to the output. For now this serves as a way to interrupt the fitting and not loose the results

7.1.20
========

Added Check_Chirp function that allows to check (look at) a chirp vector vs the data

7.1.2
========

Add a description call to the TA object meaning that "ta()" will give you a mini instruction

7.1.1
========

* implemented the option to define the external spectrum as explicit or as guidance spectrum
* Explicit Ground state added

7.0.3
========

* Allows mutli data fit with sam DAS without real fitting (for parameter check)
* move from qt to tk in all notebooks

7.0.0
========

Change to tk version as default (no qt needed)