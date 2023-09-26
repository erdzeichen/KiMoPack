*********
Changelog
*********
In this changelog I will only mention the new capabilities. All other version jumps are bug fixes

.. _release-7.2.14:

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