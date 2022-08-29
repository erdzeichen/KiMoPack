import numpy as np
import pandas as pd

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

def streak_Lund(filename):
    code=filename.split('-')
    ds=pd.read_csv(filename,sep='\t',header=None)
    n_times=len(ds.index.values)
    n_waves=len(ds.columns)
    times={'t6':np.linspace(0,2000,n_times),
           't5':np.linspace(0,1000,n_times),
           't4':np.linspace(0,500,n_times),
           't3':np.linspace(0,200,n_times),
           't2':np.linspace(0,100,n_times),
           't1':np.linspace(0,50,n_times)}
    for i in [4,3,2,5]:
        try:
            times=times[code[i]]
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