.. image:: https://readthedocs.org/projects/kimopack/badge/?version=latest
	:target: https://kimopack.readthedocs.io/en/latest/?badge=latest
	:alt: Documentation Status
	
.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
	:target: http://www.gnu.org/licenses/gpl-3.0
	:alt: License: GPL v3
	
.. image:: https://anaconda.org/erdzeichen/kimopack/badges/version.svg  
	:target: https://conda.anaconda.org/erdzeichen
	
.. image:: https://badge.fury.io/py/KiMoPack.svg
    :target: https://badge.fury.io/py/KiMoPack

.. image:: https://anaconda.org/erdzeichen/kimopack/badges/latest_release_date.svg   
	:target: https://anaconda.org/erdzeichen/kimopack
	
.. image:: https://mybinder.org/badge_logo.svg		  
	:target: https://mybinder.org/v2/gh/erdzeichen/KiMoPack/HEAD

.. image:: https://zenodo.org/badge/400527965.svg
   :target: https://zenodo.org/badge/latestdoi/400527965

KiMoPack
==========

KiMoPack is a project for the handling of spectral data measure at
multiple time-points. The current design is optimised for the use with
optical transient absorption data, but it has been successfully adapted
for the use with transient x-ray emission and spectro-electro chemistry
data.

It focuses on the main tasks an experimentator has:
Loading and shaping of experiments, plotting of experiments, comparing of experiments,
analysing experiments with fast and/or advanced fitting routines and saving/exporting/presenting 
the results. 

For typical use a series of juypter notebooks are provided that guide 
through the a number of different use scenarios, and are suggesting the 
parameter that are typically set.

Installation
--------------

The basis of the program is a module called "plot_func.py" that contains all the necessary functions and classes. 
We provide a series of jupyter based work flow packages that guide the user through a series of typical tasks 
during the analysis of optical transient absorption data and that we strongly recommend.
The files can be downloaded from the github directory https://github.com/erdzeichen/KiMoPack and manually installed (added to the path).
Alternatively we recommend the usage of the usual python install commands "pip" or if the distribution is using the Anaconda
package manager, The conda type installation. For both please open a command line (e.g. using "cmd" in windows) and execute one of the following commands. 
The notebooks can also be downloaded from the github server https://github.com/erdzeichen/KiMoPack.

Install using "pip":

.. code-block:: text

    $ pip install KiMoPack 

Upgrade if already installed:

    $ pip install KiMoPack -U

Install and update using "conda" from the channel erdzeichen:

.. code-block:: text

    $ conda install -c erdzeichen kimopack

Hint: the pip version is usually actueller than the conda version
	
Best usage
-----------
While KiMoPack is a python library, we facilitate its use with Jupyter notebooks. For the typical analysis tasks we have developed a series of Notebooks that guide through the tasks.\n 
These notebooks can be downloaded from https://github.com/erdzeichen/KiMoPack/tree/main/Workflow_tools or by command line. 

To do that start any console (under windows e.g. type "cmd" and hit enter). In the console you then start python by typing "python" and hit enter, lastly you import Kimopack and run a function that downloads the files for you by typing "import KiMoPack; KiMoPack.download_all()" This downloads the notebooks and tutorials from github for you. If you instead use "import KiMoPack; KiMoPack.download_notebooks()" then only the workflow tools are downloaded.
Please copy one of these notebooks into your data analysis folder and rename them to create a analysis log of your session. For more information please see the publication https://doi.org/10.1021/acs.jpca.2c00907, the tutorial videos, or the tutorial notebooks under https://github.com/erdzeichen/KiMoPack/tree/main/Tutorial_Notebooks_for_local_use. 
	
Citation
------------
We have written and submitted a paper introducing the toolbox under https://doi.org/10.1021/acs.jpca.2c00907

Links
-----

	* Publication: https://pubs.acs.org/doi/10.1021/acs.jpca.2c00907
	* Documentation: https://kimopack.readthedocs.io/
	* PyPI Releases: https://pypi.org/project/KiMoPack/
	* Source Code: https://github.com/erdzeichen/KiMoPack
	* Issue Tracker: https://github.com/erdzeichen/KiMoPack/issues
	* Website: https://www.chemphys.lu.se/research/projects/kimopack/
	* Zenodo: https://zenodo.org/badge/latestdoi/400527965
	* Tutorial videos: https://www.youtube.com/channel/UCmhiK0P9wXXjs_PJaitx8BQ
