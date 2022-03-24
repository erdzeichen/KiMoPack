Fitting, Parameter optimization and Error estimation
=======================================================

Fitting data:				:meth:`self.Fit_Global()<plot_func.TA.Fit_Global>`

One of the main functions of this program is to perform Global analysis
of one or multiple datasets. The fitting function is in its current
implementation embedded in the TA object and uses the parameter control
options of an lmfit parameter object as an essential tool. (my thanks to Matthew
Newville and colleagues for creating this phantastic tool) [M. Newville,
T. Stensitzki, D. B. Allen, A. Ingargiola, 2014. DOI:
10.5281/ZENODO.11813.]. To excecute an Optimization/Fit these essential
steps have to be followed (assuming that "ta" is your object):

Model
   you have to define a model by setting ta.mod="internal model name"
   for one of the internal models, of by setting
   ta.mod=external_function to an external function. For the internal
   models three standards are implemented: "exponential", "consecutive"
   and "full_consecutive". (see below). The external function will
   receive a vector named "times" and parameter pandas.DataFrame named
   "pardf". It is expected to return a DataFrame with the times as index
   and in the columns an expression of the relative concentrations. The 
   external function can reload whatever data is desired. The names of the 
   columns are used in all plots and reports. (see below a description 
   of the example functions defined in "plot_func_function_library.py")

Parameter
   For handling of Parameters I am using the lmfit Parameter objects to
   have a flexible and fast Parameter handling toolset. The general
   steps are: create a Parameter object (or use an existing parameter
   object) set starting values and (optional) limits, relative
   expressions and vary = True (if the parameter is to be optimised) (see below
   for more details

Trigger the fitting
   To trigger the fitting the function ta.Fit_Global() is called. The
   fitting function will display its results on the screen and write
   them into the TA object. first it will create an parameter object
   ta.par_fit that contains the optimum parameter and can be used later,
   second it writes the result dictionary ta.re. In this dictionary it
   adds a number of results and parameters that are useful (see below
   for further details) A number of error catching routines are build in
   the fitting software. However it does sometimes (rare) come to
   crashes. Then please adjust the parameter and or read the error
   message. Important to notice is that this is a refinement of
   parameter using an algorithm. Limiting the region where parameter can
   be and choosing good starting values is essential to obtain reliable
   results. (see below)

The general Modelling/Fitting process happens with the same approach for
internal models as well as for externally provided functions. Other
fitting processes like the optimization of the arrival time correction
are discussed further down in this document.

#. First a copy of the Data-Matrix ta.ds is created with the shaping
   parameters

#. Then a Matrix is created that represents the fractional population of
   each species (or process in case of the paral model). This Matrix
   contains one entry for each timepoint and represents the kinetic
   model based upon the starting parameter. (see below for a description
   of the models). This model formation can by done by using a build in
   or a user supplied function. (handled in the function "pf.build_c")

#. Then the process/species associated spectra for each of the species
   is calculated using the linalg.lstsq algorithm from numpy
   (https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html)

#. From the convoluted calculated species concentrations and spectra a
   calculated matrix is formed (handled in the function "pf.fill_int")

#. The difference between calculated and measured spectra is calculated,
   point-wise squared and summed together. (function "err_func")

#. This difference is minimized by iterating 2-4 with changing
   parameters using an optimization algorithm (generally nelder-mead
   simplex)

#. Finally in a last run of 2-5 the final spectra are calculated (using
   the "final" flag) and the optimized parameter, the matrixes
   ("A"-measured, "AC" - calculated, "AE" - linear error), spectra
   (always called DAS) the concentrations (called "c") are written in
   the dictionary "ta.re" together with a few representations and other
   fit outputs. The optimized parameter are also written into ta.par_fit
   (as an parameter object) that can be re-used as input into further
   optimization steps.

The choice to use lmfit parameters allows the flexible handling of
constraints and freezing/releasing of parameter. (see below). Additional
fitting options include the fitting of the parameter of the arrival time
correction, the use of advanced optimization algorithms, and the
possibility to fit multiple datasets together.

Description of models
------------------------

internal kinetic models
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The internal kinetic models are created to be highly flexible. The
program recognizes all parameter that are "rates" by recognizing the
names with "k" followed by one or multiple digits. The program starts at
"k0" and then uses "k1","k2",... in an un-interrupted series. If one key
is missing the program stops searching for more. However that offers the
fast and easy possibility to add additional components to the fit by
adding new parameters with an increasing number. All rates are to be
given as rates and are assumed to be in the same units as is the
measured matrix. So if the measured data is in picoseconds, the rates
are in :math:`ps^-1`. See section :ref:`Data shaping settings that affect the fits` 
for more information.

**ta.mod="exponential"** In this model the data is represented by
independent exponential decays. For each component the a symmetric
response function is formed (error function) using the the parameter
"resolution" as characteristic width (corresponding to 2 x
:math:`\sigma`) with with "t0" as the point of 50\ :math:`\%` rise from
0. At the time 1 x resolution the signal is 97.8\ :math:`\%`. This rise
is convoluted with an exponential decay function per given decay
constant. The integral signal of the error function is 1, to which all
the decays are set (assuming that they do not decay faster than 2x
resolution) The extracted Spectra are commonly called decay associated
spectra (DAS). Their relative intensity is corresponding with the
exponential pre-factors for single exponential fits to the data (for
each wavelength). Additional parameter that can be set are "background"
(assuming a species that is constantly "1") and "infinite" a species
that is rising with the response function and then constant "1".

**ta.mod="consecutive"** This fit assumes a model of consecutive
exponential decays. A response function with "t0" = 50\ :math:`\%` rise
is formed that rises symmetric to :math:`2\sigma \approx 98\%` (of 1) at
1 x the parameter "resolution" followed by A->B->C consecutive decay.
This particular model uses a preudo approach to this fit to speed up the
calculations. The parameter are optimised by modelling an "exponential"
model (see above) followed by a single step of a "true consecutive"
decay (see below). This approach is quite representative unless there
are fast components of the order of the response function involved in
the process and the different processes are clearly separated (each rate
one of magnitude separated). Additional parameter that can be set are
"background" (assuming a species that is constantly "1") and "infinite"
a species that is with the last decay constant to a constant "1" and not
decaying.

**ta.mod="full_consecutive"** This fit assumes a model of consecutive
exponential decays. A response function with "t0" = 50\ :math:`\%` rise
is formed that rises to :math:`2\sigma \approx 98\%` (of 1), at 1x the
parameter "resolution" followed by A->B->C consecutive decay. This model
is formed by a stepwise integrated differential equation and represents
the "true" sequential model. The "rise" is simulated by sampleing a true
gaussian function and adding the appropiate fraction to the excited
state. Arbitrary pulse/response shapes can be sampled in the advanced
functions. Additional parameter that can be set are "background"
(assuming a species that is constantly "1") and "infinite" a species
that is with the last decay constant to a constant "1" and not decaying.

external kinetic models as defined in example file "plot_func_function_library.py"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

External model functions can be written and used to create the matrix of
populations. The external function will receive a vector named "times"
and a pandas.DataFrame with the columne "value" named "pardf". It is
expected to return a DataFrame with the times as index and in the
columns the an expression of the relative concentrations. The parameters
have a name and a float value. The external function can however load
whatever other data is required. We have for example modelled
spectro-electro-chemistry data by reading the current from cyclic
voltametry and using the value to extract a number representing the
concentration of a certain species. An important feature of external
functions are that columns in the DataFrame can be labeled with names.
These names will be used in the plots and significantly improve the work
with complex models. The parameter that are given to the functions can
be name arbitrarily but must match:: 

	[a-z_][a-z 0-9_]* 

meaning, it must start with a letter and can only contain (small) letters, 
numbers and "_". Important is that in order allow for some of the parameter
settings (see section :ref:'Setting of Fit parameter') if the first
letter is a "k" and the second string is a number the parameter will be
interpreted as a rate. (recognition is done by::

	bool(re.match(re.compile(’[k]'̣), parameter_name[:2])

In the file "plot_func_function_library.py" we provide a number of
useful example functions that show how to model a variety of different
cases. In general there is no restriction to what type of function can
be modelled here, but all these examples are based upon transient
absorptions spectroscopy. In these examples we model the instrument
response by a gaussian pulse. In general, any pulse shape could be
loaded from an external file. In these examples we sample the
differential change of a concentration by writing the differential for
the dynamics. The excitation is then represented by sampling the
gaussian pulse and "raising" a certain fraction of the molecules into
the excited state. As the gaussian used here is normalized to have the
integral of "1", the total initial concentration is "1" and the numbers
in this matrix representative with a "fractional population". Each step
in the code is documented and the code can be adapted easily to a wide
variety of problems. The example functions provided are

	* 	**manual_consecutative**
		An example how a stepwise consecutative decay would look like.
	*	**Square_dependence**
		An example in which the pumping "pulse" is scaled by a parameter and
		a non linear decay is included (e.g. recombination)
	*	**gaussian_distribution**
		A model where a substance is excited into a state, followed by a
		linear decays step into a state that decays with a distribution of
		rates (here assumed gaussian) into a final, non decaying state. These
		type of complex decays are often observed in protein folding

**Usage of external functions:** To use an external function, this
function needs to be handed to ta.mod. For an external function this
means that it has to be imported, and then handed to ta.mod. In the
example below we import an external module (the file
"plot_func_function_library.py") as "func" and then use from this
external module the function "Square_dependence".
All the models are extensively documented in the function library. These
functions can load any external file with additional information. It is
highly recommended to use the versatile parameter setting (see below)
to adjust models. E.g. can a certain kinetic pathway be disabled by
setting its rate to "0" and using the option "vary=False" to lock it.
(see below in the parameter section).

It is highly recommended to use the docstring (description) directly
below the definition of the model to describe what it does. This string
is stored with ta.Save_Project and should be sufficient to identify the
model. Also if all the species are labeled (label the columns of the
returning DataFrame) These names are used throughout the plotting
functions. (please see examples for more explanation)

Remark: Importing an external function happens in python only if it has
not already been imported. So if the fitting function is adapted, either
the whole notebook/console needs to be restart, or (better) the function
should be reloaded. I recommend to use the function "reload" from the
"importlib" for this purpose (see the example below) This should happen
before the function is handed to ta.mod (as shown in the workflow
notebook).

Setting of Fit parameter
----------------------------

The fit parameter are a crucial point for achieving meaningful results
from an optimization. In general three different types need to be
chosen, first the model (see ) then if the rate parameter (necessarily
call k0, k1, k2, ..) will be handed into the fitting function as they
are or in log space. (ta.log_fit) and finally the parameter themselves.
The **log_fit** option can be important as it brings widely separated
rats into the very similar numerical range, simplifying the function of
the simplex optimizer. In this program all rates are limited to be above
0 independent if they are handled linearly or in log. This happens in
the begin of the fit function, here all "rates" are identified that have
the name "ki" with i =0-99 and then their lower limit is set to zero
(unless they have already a lower limit >0).

The parameter are handled as a lmfit Parameter object. Inside the
fitting function this object is converted into a pandas Dataframe that
is handed to the function generating the time dependent
"concentrations".

	*	initialize
		The lmfit parameter object needs to be initialized with
		"ta.par=lmfit.Parameters()". In the fitting function I convert the
		parameter object into a DataFrame and back on several places. A
		function par_to_pardf and pardf_to_par does this conversion. All the
		parameter set are available through the ta.par object and can thus be
		given to other fits. After the fit there is a new object calles
		ta.par_fit that contains the optimized fit results. So if you would
		like to re-use the old results ta.par=ta.par_fit accomplishes this.
	*	add parameter
		Each parameter must
		have a name from::
	   
			[a-z_][a-z0-9_]\* 
		  
		(starting with a letter and then
		letters and stars and "_"). In the included models (see
		:ref:'Description of models') parameters like "background"
		and "infinite" trigger the inclusion of e.g. the background or a non
		decaying component. Other parameters should be initiated with a value
		that has to be of type float (number). Important to not is that the
		code automatically recognizes parameter that have the name "ki" with
		i and element of 0-99 as a rate. These rates are brought in and out
		of logspace with the switch "log_fit". All other names can be freely
		chosen. I highly recommend to do this for the sake of structure. e.g.
		a "threshold" could be named as such
	*	add or set
		New parameter are "added" to the parameter object. Existing
		parameter can be "set" to a certain value. If Set is used any of
		the additional/optional things can be set alone.
	*	limits min and max
		Optional is the settings of limits (**min** and **max**). If a limit
		is set the parameter will stay within the limit, even if a starting
		value outside the limit is given! Important to note is that after
		each optimization that included limits the results should be checked
		if the limits were reached. (the printed output states the limits).
		Limits are very important for the more advanced optimizers like AMPGO
		(see section :ref:trigger-the-fit. The parameter tunneling
		uses these limits as guidelines.
	*	Vary=True/False
		Very useful is the option "vary=True/False". This switch freezes the
		parameter, or allows it to be optimized by the algorithm. In the
		progress of an analysis one often freezes a parameter to develop a
		stable model and releases this parameter later. Particular the
		parameter "t0" which is in my models the starting point and
		"resolution", which is in my models the instrument response function
		are parameter that are often frozen in the beginning. Fitting with
		them enabled significantly extends the duration for finding a stable
		fit. Often I first plot the function with the starting parameter,
		temporarily setting all parameter to vary=False with the trick below,
		to then step by step enable the optimization, while the starting
		parameter are adapted.
	* 	expr
		An advanced option is the setting of expressions. This are relations
		to other parameter. e.g. expr=’k0’ sets the value of the current
		parameter always the same as "k0". The values are always given as
		string so expr=’1-k0’ sets the value to 1 - the value of "k0". Please
		see the documentation of lmfit for further details

Very useful trick to set temporarily set all parameter to vary=False to
test e.g. starting conditions and then enable the optimization of a
single parameter. As here the "set" is used, the parameter can be
initially added with a different value. (see workflow notebook for
further examples).

storing of fit results
----------------------

	*	ta.par
		always contains the initial fit parameter (parameter object)
	*	ta.par_fit
		contains the fit results and can be directly re-used with
		ta.par=ta.par_fit (parameter object)
	*	ta.re[’fit_results_rates’]
		contains the fit results in a neatly formated DataFrame in the form
		of rates
	*	ta.re[’fit_results_times’]
		contains the fit results in a neatly formated DataFrame in the form
		of decay times (1/rates)
	*	ta.re[’fit_output’]
		Is the results oject of the fit routine. It can be called and then
		shows details like number of iterations, chi\ :math:`^2`, fit
		conditions and a lot more. This object is stored after a fit but is
		NOT saved by ta.Save_Project!

Trigger the Fit
----------------

The Fitting process is triggered by calling the function "Fit_Global".
if the parameter were set as part of the object that contains this Fit
(as is usually the case with ta.par), than just calling the function
without any other parameters is a good choice. Internally the Fit
function is making a copy of the parameter and shapes the data, then it
optimized single or multiple datasets. As standard it uses the Nelder
Mead Simplex algorithm to minimize the error values defined by the
function pf.err_func and pf.err_func_multi. Currently the maximum
iterations are hard-coded to be max 10000. I have not needed more than
1000 for any well defined problem. The optimizer can be changes to
"Ampgo" that offers an advanced "tunneling" algorithm for checking for
global minima. Important for this to work properly all optimizing
parameter need "min" and "max" definitions. Parameter can additionally
be given via the parameter and module input at this stage, but in
general it is better to define them as part of the ta object. The
pf.err_func and pf.err_func_multi recognise if an internal or an
external fitting model is to be used by checking if "ta.mod" (or the
here given "mod") are strings or something else (in which case it
assumes it is an external function). Additional modules from https://lmfit.github.io/lmfit-py/fitting.html 
can be easily implemented.

See section 
:ref:`external kinetic models as defined in example file "plot_func_function_library.py"`
for examples how to define those. The fitting process is in all cases
the same. Advanced options include the use of fit_chirp that runs
multiple iterations of chirp fitting and global fitting iterative (to a
maximum of fit_chirp_iterations), or the multi_project module (see
below). In general the dump_paras can be used to write into the working
directory a file with the current fitting parameter and the optimum
achieved fitting parameter. This is intended for long and slow
optimizations to keep a record of the fits even if the fitting process
did not finish.

:meth:`self.Fit_Global()<plot_func.TA.Fit_Global>`

Fitting multiple measured files at once
-----------------------------------------

To fit multiple projects the fit function needs to get a number of
projects. These can of course be opened with a hand written loop. A
cleaner way is to either use the Gui function to open a list of
files. :ref:`Opening multiple files` 
As each file needs a chirp correction and these things I
recommend to use saved projects (hdf5 files) for this purpose. Please see
the function documentation for further details. In general this function is 
fitting each of the projects separately, but using the same parameter. This means 
that in general a new (different) DAS is calculated for each of the measurements. 

To work with the same DAS for the measured and calculated matrices need to be
concatenated before the fitting. A convenient way to do this is to use the 
shaped matrices (potentially with different scattercuts) and the concentrations 
that are created after an individual fit::
	
# without scaling::

	A_con=ta_list[0].re['A']
	c_con=ta_list[0].re['c']
	
# With scaling::

	A_con=pd.concat([ta_list[0].re['A'],ta_list[1].re['A']*0.5,ta_list[2].re['A']*0.4])
	c_con=pd.concat([ta_list[0].re['c'],ta_list[1].re['c'],ta_list[2].re['c']])

# using the "fill_int" function calculates just the DAC::

	re=pf.fill_int(A_con,c_con,final=True)
	
Alternatively one could give this concatenated Matrix to the Fit_Global and construct 
an assembled external function.

Error Estimation
----------------

Estimating errors correctly is based on estimating the validity of the full set of optimized parameter for this we use the F-statistics of the single or combined datasets to define a cutoff value. At the cutoff value the combined Chi^2 is so much larger than the minimum Chi^2 that this can not be explained statistically anymore. Practically this corresponds to making the "Null hypothesis" that all parameters are zero and if the difference of Chi^2 is statistically significant, the coefficients improve the fit
the f-statistics compares the number of 
"fitted parameter"=number of species*number of spectral points + number of kinetic parameter
"free points"=number of datasets*number of spectral points*number of time points - fitted parameter
within the target quality, meaning, what fraction do my variances need to have, so that I'm 100% * target_quality sure that they are different from zero
This is done in the function :meth:`plot_func.s2_vs_smin2`. The confidence level then defines the cutoff value. For each (varied) parameter a separate optimization is performed, that attempts to find the upper and lower bound at which the total error of the re-optimized globally fitted results reaches the by F-statistics defined confidence bound. Careful, this option might run for very long time. Meaning that it typically takes 50 optimization per variable parameter (hard coded limit 200) The confidence level is to be understood that it defines the e.g. 0.65 * 100% area that the parameter with this set of values is within this bounds.

Iterative Fitting
------------------

as the fit results are written into the parameter ta.par_fit the fit can be very conveniently 
triggered in an iterative fashion. This is particularly useful for refining the chirp. 
The initially achieved optimal kinetic parameters are used as starting parameter for each
global fit after the chirp optimization. e.g. a 5 times iterative improvement can be achieved with::

	for i in range(5):
		start_error=ta.re['error']
		ta.par=ta.par_fit
		ta.Fit_Global(fit_chirp=True)
		if not ta.re['error'] < start_error:break 


