Opening of data
==========================================

A key challenge in using a non graphical programming software is to locate and open files. 
This tool provides a mixed interface to solve this challenge and in general allows three different ways to import data.

Each of the following function has a "Gui" keyword that can trigger a standard file 
opening dialogue. Alternatively the filenames can be written together with an (optional) 
path argument. 
If the analysis uses the from us provided workflow notebooks, then we suggest that a fresh notebook
is used for each particular analysis and that the notebook is copied close to the data. 

Three different pathways of importing data are offered:

1. All import functions provide a wide variety of options to adapt for data formats. If a particular option is missing and is desired, please contact the developers via email or raise an issue on github. We will consider to do so but most likely instead provide you with a function for option 2. The formats of the example file have the spectral information as the first row, the time vector as first entrance of each of the following rows and are separated by tab. Files of this type can be read without any further adaption (using the standard parameter). typical options include the transposing of columns, the conversion of time and energy vectors or the providing of external files for energy and times.

2. All import function have the option of providing an external import function. (from version 7.8.0 onwards) this function gets the filename and should return a dataframe to KiMoPack. We provide a function library that contains the formats of befriended groups. If you would like help to develop an import function then please contact the developers.

3. The two main import function have a "ds" parameter to which a Pandas DataFrame can be given. Thus the user might simply import and shape the file and then hand it over to KiMoPack.  

Opening single file and creating TA object
------------------------------------------

Open single file: 			:meth:`pf.TA()<plot_func.TA.__init__>`

Typical use of the tool is based upon an object containing all
parameter, functions and the data. This object is created by importing a
data file using this function.

The filename can be either provided as a string, with the (optional) path to other folders given.
Using the keyword "gui" instead of a filename opens the graphical file selection interface.
The function can either open a text style format (using any file ending) or the internally used "hdf5" file format. 
The latter is exclusively used as internal storage formate that stores the complete project including the RAW data.

After import of either filetype the missing parameter in the "TA" object are set with the 
:meth:`self.__make_standard_parameter()<plot_func.TA.__make_standard_parameter>` function. 

Opening multiple files
----------------------------

Open many files: 			:meth:`pf.GUI_open()<plot_func.GUI_open>`

Sometimes multiple files are to be opened. Typical use examples are the options to compare different 
measurements or analysis runs. This function provides a convenient way to create a list of opened projects. 
One can

	* open a gui and select multiple saved projects, which are returned as a list
	* given a list of project names to open them
	* open all files in a given folder 

The general behavior is selected by the first parameter (project_list)
For more details see the examples in :meth:`pf.GUI_open()<plot_func.GUI_open>`

Opening and handling single scans
----------------------------------

Combine many scans:			:meth:`pf.Summarize_scans()<plot_func.Summarize_scans>`

Typically the experiments consists of a larger number of scans that are combined into a single experimental file.
The function "Summarize_scans" reads, selects and eventually combines a
series of single scans with a bunch of useful options. The essential idea is
that for each scan one or two numbers are generated through integration of the intensity 
in a temporal and spectral window. This single number is plotted as function against the scan number. 
Then either a list of numbers or a GUI is used to select the scans that are 
removed from the average. A number of opening and selection options are given.

Observe the new automatic filter options of summarize_scans

This function could also be used to combine a number of different experiments.
