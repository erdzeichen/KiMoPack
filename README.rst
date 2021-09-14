.. image:: https://readthedocs.org/projects/kimopack/badge/?version=latest
	:target: https://kimopack.readthedocs.io/en/latest/?badge=latest
	:alt: Documentation Status
	
.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
	:target: http://www.gnu.org/licenses/gpl-3.0
	:alt: License: GPL v3
	
.. image:: https://anaconda.org/erdzeichen/kimopack/badges/installer/conda.svg   
	:target: https://conda.anaconda.org/erdzeichen
	
.. image:: https://badge.fury.io/py/KiMoPack.svg
    :target: https://badge.fury.io/py/KiMoPack
	
.. image:: https://anaconda.org/erdzeichen/kimopack/badges/latest_release_date.svg   
	:target: https://anaconda.org/erdzeichen/kimopack
	


KiMoPack
==========

KiMoPack is a project for the handling of spectral data measure at
multiple time-points. The current design is optimised for the use with
optical transient absorption data, but it has been successfully adapted
for the use with transient x-ray emission and spectro-electro chemistry
data.

It focuses on the main tasks an experimentator has
Loading and shaping of experiments, plotting of experiments, comparing of experiments,
analysing experiments with fast and/or advanced fitting routines and saving/exporting/presenting 
the results. 

For typical use a series of juypter notebooks are provided that guide 
through the a number of different use scenarios, and are suggesting the 
parameter that are typically set.

Please note that the uploaded version is currently a dummy that will be replaced with the fully functional version at the moment of publication
===============================================================================================================================================

Installion
----------

The basis of the program is a module called "plot_func.py" that contains all the necessary functions and classes. 
We provide a series of jupyter based work flow packages that guide the user through a series of typical tasks 
during the analysis of optical transient absorption data and that we strongly recommend.
The files can be downloaded from the github directory https://github.com/erdzeichen/KiMoPack and manually installed (added to the path).
Alternatively we recommend the usage of the usual python install commands "pip" or if the distribution is using the Anaconda
package manager, The conda type installation. For both please open a command line (e.g. using "cmd" in windows) and execute the following commands. 
The Jupyter notebooks are copied during the install process. The notebooks can also be downloaded from the github server https://github.com/erdzeichen/KiMoPack.

Install and update using "pip":

.. code-block:: text

    $ pip install KiMoPack

Install and update using "conda" from the channel erdzeichen:

.. code-block:: text

    $ conda install -c erdzeichen kimopack
	
After publication we will also add zenodo and a conda-forge branch



Links
-----

	* Documentation: https://kimopack.readthedocs.io/
	* PyPI Releases: https://pypi.org/project/KiMoPack/
	* Source Code: https://github.com/erdzeichen/KiMoPack
	* Issue Tracker: https://github.com/erdzeichen/KiMoPack/issues
	* Website: https://www.chemphys.lu.se/research/projects/kimopack/

