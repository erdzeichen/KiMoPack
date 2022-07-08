import numpy as np
import pandas as pd

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