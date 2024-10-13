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

.. image:: https://colab.research.google.com/assets/colab-badge.svg
	:target: https://colab.research.google.com/github/erdzeichen/KiMoPack/blob/main/Tutorial_Notebooks/KiMoPack_tutorial_0_Introduction.ipynb

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
We recommend to use a package manager to install the program.  

Install using "pip":

.. code-block:: text

    $ pip install KiMoPack 

Install and update using "conda" from the channel erdzeichen:

.. code-block:: text

    $ conda install -c erdzeichen kimopack

Hint: the pip version is usually more recent than the conda version
The files can also be downloaded from the github directory https://github.com/erdzeichen/KiMoPack or zenodo (see below)

These commands are installing only KiMoPack and the absolutely needed dependencies. However, there are several modules that work better if additional packages are installed. In general one should try to install all packages at the same time. Additional packages that I generally recommend are h5py and tables (for saving files), python-pptx (for saving power point slides) and keyboard (Window only, for interrupting the fits). Quite useful is also nbopen that allows you to double click on the notebook files. nbopen requires an activation at the end.

(leave away keyboard for Linux!)
.. code-block:: text

    $ pip install KiMoPack h5py tables nbopen python-pptx 
	(windows) python -m nbopen.install_win
	(Linux) python3 -m nbopen.install_xdg
	(MacOS) Clone the repository and run ./osx-install.sh

Upgrade if already installed:

.. code-block:: text

    $ pip install KiMoPack -U


In general it is a good idea to create a local environment to install files in python if you are using python for many tasks. In a local environment only the packages that are needed are installed, which usually avoids that conflicts can appear. It is very easy to do that. 

Under Windows: open the anaconda command prompt or power shell (type anaconda under windows start) 
Under Linuxs: open a console

.. code-block:: text

	$ conda create --name kimoPack
	$ conda activate kimokack
	$ pip install KiMoPack h5py tables keyboard nbopen python-pptx

Or if you also want make sure to have a later version of python	

.. code-block:: text

	$ conda create --name kimopack python=3.11 ipython jupyterlab jupyter
	$ conda activate kimopack
	$ pip install KiMoPack h5py tables keyboard nbopen python-pptx


Error: insufficient rights: If one of the installs complains (error) that the user does not has sufficient rights, this installation can be done attaching "--user"

.. code-block:: text

	$ conda create --name kimoPack
	$ conda activate kimokack
	$ pip install KiMoPack h5py tables keyboard nbopen python-pptx --user

Error: pytables:
	in some versions I have been running in a problem with pytables when loading saved data. 
	Using the conda forge version solved this problem for me 

.. code-block:: text
	
	conda install -c conda-forge pytables  

Best usage
-----------
While KiMoPack is a python library, we facilitate its use with Jupyter notebooks. For the typical analysis tasks we have developed a series of Notebooks that guide through the tasks.\n 
These notebooks can be downloaded from https://github.com/erdzeichen/KiMoPack/tree/main/Workflow_tools or by command line. 

You can try either of these "lazy" oneliners

.. code-block:: text

	ipython -c "import KiMoPack; KiMoPack.download_notebooks()"
	python -c "import KiMoPack; KiMoPack.download_notebooks()"
	python3 -c "import KiMoPack; KiMoPack.download_notebooks()"

If none of these work then start any console (under windows e.g. type "cmd" and hit enter). In the console you then start python by typing "python" and hit enter, lastly you import Kimopack and run a function that downloads the files for you by typing "import KiMoPack; KiMoPack.download_all()" This downloads the notebooks and tutorials from github for you. If you instead use "import KiMoPack; KiMoPack.download_notebooks()" then only the workflow tools are downloaded.
Please copy one of these notebooks into your data analysis folder and rename them to create a analysis log of your session. For more information please see the publication https://doi.org/10.1021/acs.jpca.2c00907, the tutorial videos, or the tutorial notebooks under https://github.com/erdzeichen/KiMoPack/tree/main/Tutorial_Notebooks_for_local_use. 
	
Citation
------------
We have published a paper introducing the toolbox under https://doi.org/10.1021/acs.jpca.2c00907

Links
-----
	* Overview talk: I gave a recent overview talk at the LaserLab Europe meeting: https://youtu.be/z9QqVLFWYrs
	* Publication: https://pubs.acs.org/doi/10.1021/acs.jpca.2c00907
	* Documentation: https://kimopack.readthedocs.io/
	* PyPI Releases: https://pypi.org/project/KiMoPack/
	* Source Code: https://github.com/erdzeichen/KiMoPack
	* Issue Tracker: https://github.com/erdzeichen/KiMoPack/issues
	* Website: https://www.chemphys.lu.se/research/projects/kimopack/
	* Zenodo: https://zenodo.org/badge/latestdoi/400527965
	* Tutorial videos: https://www.youtube.com/channel/UCmhiK0P9wXXjs_PJaitx8BQ
