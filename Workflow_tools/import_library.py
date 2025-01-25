import numpy as np
import pandas as pd
import os

def NRL(filename):
	ds=pd.read_csv(filename, sep=',', index_col=0, header=None)
	ds.columns=ds.iloc[0,:]
	ds.drop(0,inplace=True)
	ds.index=ds.index.astype(float)
	ds.columns=ds.columns.astype(float)
	ds.sort_index(inplace=True,axis=0)
	ds.sort_index(inplace=True,axis=1)
	ds.astype('float')
	ds=ds.T
	ds.index.name='Time in ps'
	ds.columns.name='Wavelength in nm'
	return ds,'differential absorption in mOD','ps'

def Uppsala(filename):
	df=pd.read_csv(filename,index_col=0,sep=',').T
	#df=df.fillna(0)
	df=df.dropna(axis=1)
	df.index=df.index.astype('float')
	df.columns=df.columns.astype('float')
	df.sort_index(inplace=True)
	df.sort_index(inplace=True,axis=1)
	df.astype('float')
	df.index.name='Time in ps'
	df.columns.name='Wavelength in nm'
	return df,'differential absorption in mOD','ps'


def Ivan_horse(filename):
	#print(filename)
	import scipy.constants as const
	ds=pd.read_csv(filename,sep='\t',index_col=0)
	ds.index=ds.index.astype(float)
	ds.columns=ds.columns.astype(float)
	#ds.index=ds.index.values/1e8
	#ds.columns=ds.columns.values/1e5
	#ds=ds.apply(lambda x:x*ds.index.values)
	#ds=ds.apply(lambda x:x*ds.columns.values,axis=1)
	#per_photon=const.h*const.c/(485e-9)
	#ds=	ds*per_photon
	ds.index.name='Fluence in Photons/cm2 s'
	ds.columns.name='Repetitions rate in Hz'
	ds.sort_index(inplace=True)
	ds.sort_index(inplace=True,axis=1)
	return ds.T,'PLQY'

def cor_streak_lund(filename):
	'''This reads the file of the "DAC" that is exported as corrected file from the streak camera software ending'''
	ds=pd.read_csv(filename,sep='\t',index_col=0)
	ds.columns.name="nm"
	ds.index.name="Time in ps"
	data_type="Emission intensity"
	baseunit="ps"
	ds.index=ds.index.astype(float)
	ds.columns=ds.columns.astype(float)
	ds.sort_index(inplace=True,axis=1)
	ds.sort_index(inplace=True,axis=0)
	return ds,data_type,baseunit

def streak_Lund(filename):
	'''This is reading the filetype that is saved by the streak camera software as "dat" type'''
	code=str(filename).split(os.sep)[-1]#split of the path
	code=code.split('.')[0]				#split of the dot and fileending
	code=code.split('-')
	ds=pd.read_csv(filename,sep='\t',header=None)
	n_times=len(ds.index.values)
	n_waves=len(ds.columns)
	times={	't6':2000,
			't5':1000,
			't4':500,
			't3':200,
			't2':100,
			't1':50}
	for i in [4,3,2,5]:#position can change
		try:
			times=times[code[i]]
			times=np.linspace(0,times,n_times)
			break
		except:
			pass
	for i in [6,5,7,8]:
		if code[i][0]=='w':
			center=code[i][1:]
			waves=np.linspace(float(center)-60,float(center)+75,n_waves)
			break
		else:
			continue
	ds.index=times
	ds.index=ds.index.astype(float)
	ds.columns=waves
	ds.columns=ds.columns.astype(float)
	ds.index.name='Time in ps'
	ds.columns.name='Wavelength in nm'
	return ds,'emission intensity'
	
def Amine_func(filename):
	df=pd.read_csv(filename,sep='\t',header=None)
	wavelength=pd.Series(np.linspace(343.33,656.03,512))
	time=pd.Series(np.linspace(0,50.500,512))
	df.columns=wavelength
	df.index=time
	df.index=df.index.astype(float)
	df.columns=df.columns.astype(float)
	df.index.name='Time in ns'
	df.columns.name='Wavelength in nm'
	return df,'differential absorption','ns'
	
def Lund_colors():
	cols={
	'black':[0,0,0],  		
	'wine_red':[152,30,50],	
	'brown':[153,102,51],	
	'dark_blue':[0,0,128],	
	'gold':[233,131,0],
	'darker grey':[146,139,129],
	'green':[85,118,48],
	'light_grey':[203,199,191],
	'light blue grey':[189,203,197],
	'light green grey':[199,210,138],
	'pink':[233,196,199], 
	'blue':[185,211,220],  		
	'light_green':[173,202,184],  		
	'yellow':[214,210,196],  		
	'light_brown':[191,184,175],  		
	'pastel pink':[219,173,177],  
	'pastel blue':[164,196,207],  		
	'pastel green':[153,190,167],  		
	'pastel yellow':[203,197,169],  		
	'pastel brown':[180,168,154],
	'white':[255,255,255],	
	}
	cols=pd.DataFrame(cols).T/255.
	cols.columns=['R','G','B']
	return cols
