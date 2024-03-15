import numpy as np
import pandas as pd
import os



def Eric_data(filename):
	#define the path to the files
	path_to_wl=os.path.dirname(os.path.realpath(__file__)) 	# using this files path as start, type the path according to your system!
	wl=os.sep.join([path_to_wl,'Data','Practical','wl.txt'])			# combine the filepath with the name of the file
	time=os.sep.join([path_to_wl,'Data','Practical','time_ps.txt'])		# combine the filepath with the name of the file
	
	#do the general reading
	wl=pd.read_csv(wl,header=None).values					# using pandas reading function to get the values
	time=pd.read_csv(time,header=None).values		# using pandas reading function to get the values
	data=pd.read_csv(filename,sep='\s+',header=None)
	
	#do the shaping
	data.index=np.squeeze(wl)								#convert to single dimension and set as index
	data.columns=np.squeeze(time)
	data=data.T
	data.columns.name='nm'
	data.index=(data.index.values/1000.)+0.82
	return data, 'differential absorption in mOD', 'ns' # return data, data_type, units