import urllib3
import sys
import os
import pathlib
import shutil
from pathlib import Path

def check_folder(path = None, current_path = None, filename = None):
	'''Helper function using robust path determination.\n 
	In any case if a valif file name is given it is attached to the total path\n
	The path can be string or windows/linux path or pure path or byte type paths.\n
	paths that do not exists (including parents) are created\n
	1. if path is given absolute, it is returned\n_colors
	2. if path is a string (relative) the current_path + path is returned.\n
	3. if current_path is not absolute or None, the current working directory is assumed as path.\n
	4. IF all is None, the current working directory is returned
	
	Parameters
	-----------
	
	path : str, purePath, absolute or relative, optional
		the final part of the path used
	
	current_path : None, str, purePath, absolute, optional
		path that sits before the "path variable, is filled with current working directory if left None
	
	filename: None, str, optional
		attached after path and returned if not None
		
	'''
	
	if isinstance(path,bytes):
		path = '%s'%path
	if path is not None:
		path = pathlib.Path(path)

	if isinstance(current_path, bytes):
		current_path = '%s'%current_path
	if current_path is not None:
		current_path=pathlib.Path(current_path)
		
	if isinstance(filename, bytes):
		filename='%s'%filename
	if filename is not None:
		filename = pathlib.Path(filename)
	if path is None:
		if current_path is None:
			directory = Path.cwd()
		elif current_path.is_absolute():
			directory=current_path
		else:
			print('attention, current_path was given but not absolute, replaced by cwd')
			directory = Path.cwd()
	elif path.is_absolute():
		directory = path
	else:
		if current_path is None:
			directory = Path.cwd().joinpath(path)
		elif current_path.is_absolute():
			directory = current_path.joinpath(path)
		else:
			print('attention, current_path was given but not absolute, replaced by cwd')
			directory = Path.cwd().joinpath(path)
	directory.mkdir( parents=True, exist_ok=True)
	if filename is None:
		return directory
	else:
		return directory.joinpath(filename)
def download_notebooks():
	'''function loads the workflow notebooks into the active folder'''
	http = urllib3.PoolManager()
	list_of_tools=['TA_Advanced_Fit.ipynb',
					'TA_comparative_plotting_and_data_extraction.ipynb',
					'TA_Raw_plotting.ipynb',
					'TA_Raw_plotting_and_Simple_Fit.ipynb',
					'TA_single_scan_handling.ipynb',
					'Function_library_overview.pdf',
					'function_library.py',
					'import_library.py']
	print('Now downloading the workflow tools')
	for f in list_of_tools:
		url = "https://raw.githubusercontent.com/erdzeichen/KiMoPack/main/Workflow_tools/%s"%f
		print('Downloading Workflow Tools/%s'%f)
		with open(check_folder(path = 'Workflow_tools', current_path = os.getcwd(), filename = f), 'wb') as out:
			r = http.request('GET', url, preload_content=False)
			shutil.copyfileobj(r, out)
def download_all():
	''' function loads workflow notebooks and example files and tutorials'''
	download_notebooks()
	http = urllib3.PoolManager()
	list_of_tools=['TA_Advanced_Fit.ipynb',
					'TA_comparative_plotting_and_data_extraction.ipynb',
					'TA_Raw_plotting.ipynb',
					'TA_Raw_plotting_and_Simple_Fit.ipynb',
					'TA_single_scan_handling.ipynb',
					'Function_library_overview.pdf',
					'function_library.py',
					'import_library.py']
	print('Now downloading the workflow tools and tutorials')
	list_of_example_data=['sample_1_chirp.dat',
							'Sample_2_chirp.dat',
							'sample_1.hdf5',
							'sample_2.hdf5',
							'Sample_1.SIA',
							'Sample_2.SIA']
	print('Now downloading the example files')
	for f in list_of_example_data:
		url = "https://raw.githubusercontent.com/erdzeichen/KiMoPack/main/Workflow_tools/Data/%s"%f
		print('Downloading Workflow Tools/Data/%s'%f)
		with open(check_folder(path = 'Workflow_tools'+os.sep+'Data', current_path = os.getcwd(), filename = f), 'wb') as out:
			r = http.request('GET', url, preload_content=False)
			shutil.copyfileobj(r, out)
	print('Now downloading zipfile with tutorials')
	url = "https://raw.githubusercontent.com/erdzeichen/KiMoPack/main/Workflow_tools/Tutorial_Notebooks_for_local_use.zip"
	with open(check_folder(path = 'Tutorials', current_path = os.getcwd(), filename = "Tutorial_Notebooks_for_local_use.zip"), 'wb') as out:
		r = http.request('GET', url, preload_content=False)
		shutil.copyfileobj(r, out)