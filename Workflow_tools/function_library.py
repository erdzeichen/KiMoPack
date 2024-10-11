import pandas
import scipy
import numpy as np
from numpy import power,log10,shape,exp
from scipy.special import erf
import warnings 
# suppress warnings for overflow
warnings.filterwarnings('ignore') 

FWHM=2.35482

def gauss(t,sigma=0.1,mu=0,scale=1):
	y=np.exp(-0.5*((t-mu)**2)/sigma**2)
	y/=sigma*np.sqrt(2*np.pi)
	return y*scale
def rise(x,sigma=0.1,begin=0):
	'''	 my own implementation of the instrument response function. 
		 Based upon an error function from 0 to 1. 
		 Sigma is the width (after which it has 50%) 
		 and begin is 10% of height'''
	return (erf((x-sigma-begin)*np.sqrt(2)/(sigma))+1)/2

def osc_split_sin_cos(times,pardf,comp=0):
	'''returns a damped osciallation that used the "comp" parameter as index
	if comp=0 then
	f0 is frequency and the only mandatory entry
	A0 is amplitude (default=1)
	S0 fraction of cos while 1-S0 is the fraction of sin'''
	if not 'f%i'%comp in list(pardf.index.values):
		raise ValueError('frequency f%i is the minimum required parameter'%comp)
	if 'S%i'%comp in list(pardf.index.values):
		S=pardf['S%i'%comp]
		if S>1 or S<0:
			raise ValueError('fraction S%i must be between 1 and 0'%comp)
	else:
		S=1
	oscil1=S*np.cos(2*np.pi*times/pardf['f%i'%comp])
	oscil2=(1-S)*np.sin(2*np.pi*times/pardf['f%i'%comp])
	oscil1=pandas.Series(oscil1,index=times)
	oscil2=pandas.Series(oscil2,index=times)
	oscil1.name='os%i_cos'%comp
	oscil2.name='os%i_sin'%comp
	return oscil1,oscil2

def oscil_comp(times,pardf,max_oscil=None):
	'''counts how many osciallations are in there (assuming max 10 as guess) 
	if more need to be fitted the parameter max_oscil must be set
	then uses the function osc_split_sin_cos to calculate oscialltions and lets them decay 
	with the parameter tk_i where "i" is the same index as f_i the frequencies'''
	if max_oscil is None: #If we don't know how many osciallations we have we might want to find that number
		for a in range(20):# lets assume we have max 20 osciallations and find out how many
			if 'f%i'%a in list(pardf.index.values):
				max_oscil=a
				if 'tk%i'%a in list(pardf.index.values):
					continue
				else:
					pardf['tk%i'%a]=0
			else:
				break
	count_oscil=max_oscil+1
	listen=['tk%i'%a for a in range(count_oscil)]
	params=pardf[listen].values.reshape(1,-1)
	inner=np.tile(times-pardf['t0'],(count_oscil,1)).T*params.astype(float)
	c=exp(-1*inner)
	c[(times-pardf['t0'])<0]=1
	c*=np.tile(rise(x=times,sigma=pardf['resolution'],begin=pardf['t0']),(count_oscil,1)).T
	c=pandas.DataFrame(c,index=times)
	c.columns=['os%i_cos'%a for a in range(count_oscil)]
	c.index.name='time'
	for a in range(count_oscil):
		oscil1,oscil2=osc_split_sin_cos(times,pardf,comp=a)
		c['os%i_sin'%a]=c.loc[:,'os%i_cos'%a]*oscil2.values
		c.loc[:,'os%i_cos'%a]=c.loc[:,'os%i_cos'%a]*oscil1.values
	c.fillna(method='backfill')
	return c
	
def exponential(times,pardf):
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10
	param=[]
	pardf.sort_index(inplace=True)
	for key in pardf.index.values:
		if key[0] == 'k':#its a time parameter
			param.append(float(pardf[key]))
	t0=float(pardf['t0'])
	resolution=float(pardf['resolution'])
	c=np.exp(-1*np.tile(times-t0,(len(param),1)).T*param)
	c[(times-t0)<0]=1
	c*=np.tile(rise(x=times,sigma=resolution,begin=t0),(len(param),1)).T
	c=pandas.DataFrame(c,index=times)
	c.index.name='time'
	if 'infinite' in list(pardf.index.values):
		c['infinite']=rise(x=times,sigma=resolution,begin=t0)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
		#if 'infinite' in list(pardf.index.values):#not sure yet
		#c['GS']=-c.iloc[:,:-1].sum(axis=1)
	if 'background' in list(pardf.index.values):
		c['background']=1
	return c
	
def consecutive(times,pardf):
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10
	param=[]
	pardf.sort_index(inplace=True)
	for key in pardf.index.values:
		if key[0] == 'k':#its a time parameter
			param.append(float(pardf[key]))
	t0=float(pardf['t0'])
	resolution=float(pardf['resolution'])
	n_decays=len(param)
	g=gauss(times,sigma=resolution/FWHM,mu=t0)
	if 'infinite' in list(pardf.index.values):
		infinite=True
		n_decays+=1
	else:
		infinite=False
	decays=param
	c=np.zeros((len(times),n_decays),dtype='float')
	if 'explicit_GS' in list(pardf.index.values):
		GS=True
		gs=np.zeros((len(times),1),dtype='float')
	else:
		GS=False
	for i in range(1,len(times)):
		dc=np.zeros(n_decays,dtype='float')
		dt=(times[i]-times[i-1])/(sub_steps)
		c_temp=c[i-1,:]
		for j in range(int(sub_steps)):
			for l in range(0,n_decays):
				if l>0:
					if infinite:
						if l<(n_decays-1):
							dc[l]=decays[l-1]*dt*c_temp[l-1]-decays[l]*dt*c_temp[l]
						else:
							dc[l]=decays[l-1]*dt*c_temp[l-1]
					else:
						dc[l]=decays[l-1]*dt*c_temp[l-1]-decays[l]*dt*c_temp[l]
				else:
					if infinite and n_decays==1:
						dc[l]=g[i]*dt
					else:
						dc[l]=g[i]*dt-decays[l]*dt*c_temp[l]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])					
		c[i,:] =c_temp
		if GS and not infinite:
			gs[i]=-c[i,:].sum()
		elif GS:
			gs[i]=-c[i,:-1].sum()
	c=pandas.DataFrame(c,index=times)
	c.index.name='time'
	c.columns=['species %i'%i for i in range(len(c.columns))]
	if infinite:
		labels=list(c.columns.values)
		labels[-1]='Non Decaying'
		if 'background' in list(pardf.index.values):
			c['background']=1
	else:
		if 'background' in list(pardf.index.values):
			c['background']=1
	if GS:
		c['GS']=gs
	return c

def consec_oscil(times,pardf):
	c=pandas.concat([consecutive(times,pardf),oscil_comp(times,pardf)],axis=1)
	return c

def manconsec_oscil(times,pardf):
	c=pandas.concat([manual_consecutive(times,pardf),oscil_comp(times,pardf)],axis=1)
	return c

def expon_oscil(times,pardf):
	c=pandas.concat([exponential(times,pardf),oscil_comp(times,pardf)],axis=1)
	return c



def manual_consecutive(times,pardf):								
	'''we have 1MLCT,3MLCT,3MC state with 
	consecutive decays, hint this function could have been created with the linear exp model
	'''
	c=np.zeros((len(times),3),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10 													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		#dc=np.zeros((3,1),dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dc=np.zeros(3,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt					#excite a small fraction with g[i] and decay with 'k0'
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]	#form with "k0" and decay with "k1"
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]	#form with "k1" and decay with "k2"
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=[1,2,3]								#this is optional but very useful. The species get names that represent some particular states
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c														#return the concentrations to the global fitting
	
def Square_dependence(times,pardf):
	'''initial A then two paths to C one over B and one direct, last paramter is the ratio'''
	c=np.zeros((len(times),3),dtype='float')						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0'])*pardf['f0']	#creating the gaussian pulse that will "excite" our sample. Note the additional fraction with "f0". This fraction is breaking the normalization but allows the pump fluence to be added
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(3,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-2*pardf['k2']*dt*c_temp[0]**2+g[i]*dt #excite a small fraction with g[i] and decay with 'k0' linear, two state decay with the square dependence and "k2"
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]	#form with "k0" and decay with "k1"
			dc[2]=pardf['k2']*dt*c_temp[0]**2						#one single part of c[2] is formed from the non linear combination of two c[0]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A*','B*','C']										#this is optional but very useful. The species get names that represent some particular states
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c														#return the concentrations to the global fitting
	
def gaussian_distribution(times,pardf):
		#first attempt, we have one decay, then the gauss, then f0 is the spread of the gauss in rate
		#so G - with pulse to A to gauss B to intermediate C and back to G
		decays=3
		spread_shape=gauss(np.linspace(-1,1,91),sigma=1/3.,mu=0)					#this vector holds the fresh distribution
		rate_spread=pardf['k1']+np.linspace(-3*pardf['rate_spread']/FWHM,3*pardf['rate_spread']/FWHM,91) #this vector has the same shape as the distribution and 
																										#holds one rate entry in the spread_shape
		spread=np.zeros(spread_shape.shape) 															#initially there is nothing in the spread matrix
		c=np.zeros((len(times),decays),dtype='float') 													#here I define number of concentrations
		g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0'])									#this is my pump
		if 'sub_steps' in list(pardf.index.values):
			sub_steps=pardf['sub_steps']
		else:
			sub_steps=10 																					#We sample with 10 steps per measured timestep
		for i in range(1,len(times)):
			dc=np.zeros(c.shape[1],dtype='float')													# this contains the usual concentration differences
			c_temp=c[i-1,:]																				#load the previous concentrations (absolute)
			dt=(times[i]-times[i-1])/(sub_steps)														#create the momentary timestep	
			for j in range(int(sub_steps)):
				dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt 												#C0 is filled by the pump and decays with k0				
				spread+=spread_shape*pardf['k0']*dt*c_temp[0]-spread_shape*rate_spread*dt				#whatever decays from C0 flows into the distribution, 
																										#important on the new stuff is distributed with the gaussian, 
																										#each unit has its own flowing out rate
				dc[1]=0																					#set to 0 because we do use the matrix later to record the change
				dc[2]=(spread_shape*rate_spread*dt).sum()												#whatever flows out of the C1 (the distrubution) is collected into 
				c_temp=c_temp+dc
				c_temp[c_temp<0]=0
				#for b in range(c.shape[1]):
				#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])										#check that nothing will be below 0 (concentrations)
			c[i,:] =c_temp
			c[i,1] =spread.sum()																		#here we fill the record matrix with the sum of the units
		c=pandas.DataFrame(c,index=times)
		c.index.name='time'
		if 'explicit_GS' in list(pardf.index.values):
			c['GS']=-c.sum(axis=1)
		if 'background' in list(pardf.index.values):													#we still might want to have a backgraound
			c['background']=1
		c.columns=['initial','gauss_populated','final']													#this is optional but very useful. The species get names that represent 
																										#some particular states
		return c

def P12(times,pardf):								
	'''P12'''
	c=np.zeros((len(times),3),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(3,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k2']*dt*c_temp[0]+g[i]*dt					
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]+pardf['k2']*dt*c_temp[0]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','Inf']									    #this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P13(times,pardf):								
	'''P13'''
	c=np.zeros((len(times),4),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(4,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k4']*dt*c_temp[0]-pardf['k2']*dt*c_temp[0]+g[i]*dt					
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[0]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k1']*dt*c_temp[1]+pardf['k3']*dt*c_temp[2]+pardf['k4']*dt*c_temp[0]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P14(times,pardf):								
	'''P14'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k6']*dt*c_temp[0]-pardf['k0']*dt*c_temp[0]-pardf['k2']*dt*c_temp[0]-pardf['k4']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[0]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k4']*dt*c_temp[0]-pardf['k5']*dt*c_temp[3]
			dc[4]=pardf['k6']*dt*c_temp[0]+pardf['k1']*dt*c_temp[1]+pardf['k3']*dt*c_temp[2]+pardf['k5']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P21(times,pardf):								
	'''P21'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k4']*dt*c_temp[0]-pardf['k0']*dt*c_temp[0]-pardf['k2']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[0]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k4']*dt*c_temp[0]+pardf['k1']*dt*c_temp[1]+pardf['k3']*dt*c_temp[2]-pardf['k5']*dt*c_temp[3]
			dc[4]=pardf['k5']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P22(times,pardf):								
	'''P22'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[1]-pardf['k4']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[1]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k4']*dt*c_temp[1]-pardf['k5']*dt*c_temp[3]
			dc[4]=pardf['k1']*dt*c_temp[1]+pardf['k3']*dt*c_temp[2]+pardf['k5']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P23(times,pardf):								
	'''P23'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k2']*dt*c_temp[1]-pardf['k4']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[1]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k4']*dt*c_temp[1]-pardf['k5']*dt*c_temp[3]
			dc[4]=pardf['k1']*dt*c_temp[0]+pardf['k3']*dt*c_temp[2]+pardf['k5']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P24(times,pardf):								
	'''P24'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k2']*dt*c_temp[0]-pardf['k4']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[0]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k1']*dt*c_temp[1]+pardf['k3']*dt*c_temp[2]-pardf['k5']*dt*c_temp[3]
			dc[4]=pardf['k4']*dt*c_temp[0]+pardf['k5']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c	
		
def P31(times,pardf):								
	'''P31'''
	c=np.zeros((len(times),4),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(4,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k3']*dt*c_temp[0]+g[i]*dt					
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]
			dc[3]=pardf['k2']*dt*c_temp[2]+pardf['k3']*dt*c_temp[0]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c
	
def P32(times,pardf):								
	'''P32'''
	c=np.zeros((len(times),4),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(4,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k2']*dt*c_temp[0]+g[i]*dt					
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]+pardf['k2']*dt*c_temp[0]-pardf['k3']*dt*c_temp[2]
			dc[3]=pardf['k3']*dt*c_temp[2]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c
		
def P33(times,pardf):								
	'''P33'''
	c=np.zeros((len(times),4),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(4,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt					
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]-pardf['k3']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]
			dc[3]=pardf['k3']*dt*c_temp[1]+pardf['k2']*dt*c_temp[2]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','Inf']									#this is optional but very useful. The species get names that represent some particular states
	c.loc[:,'GS']=c.loc[:,'GS'].values*(-1)
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P34(times,pardf):								
	'''P33'''
	c=np.zeros((len(times),3),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10 													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(3,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt					
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[1]
			dc[2]=pardf['k2']*dt*c_temp[1]+pardf['k1']*dt*c_temp[1]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c
		
def P41(times,pardf):								
	'''P41'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10 													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k4']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]+pardf['k4']*dt*c_temp[0]-pardf['k2']*dt*c_temp[2]-pardf['k5']*dt*c_temp[2]
			dc[3]=pardf['k2']*dt*c_temp[2]-pardf['k3']*dt*c_temp[3]
			dc[4]=pardf['k5']*dt*c_temp[2]+pardf['k3']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def P42(times,pardf):								
	'''P42'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]-pardf['k3']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]
			dc[3]=pardf['k2']*dt*c_temp[2]+pardf['k3']*dt*c_temp[1]-pardf['k4']*dt*c_temp[3]
			dc[4]=pardf['k4']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c
		
def P43(times,pardf):								
	'''P43'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k3']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]
			dc[3]=pardf['k2']*dt*c_temp[2]+pardf['k3']*dt*c_temp[0]-pardf['k4']*dt*c_temp[3]
			dc[4]=pardf['k4']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c	
		
def P44(times,pardf):								
	'''P44'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]-pardf['k4']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]
			dc[3]=pardf['k2']*dt*c_temp[2]-pardf['k3']*dt*c_temp[3]
			dc[4]=pardf['k4']*dt*c_temp[1]+pardf['k3']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c	

def P45(times,pardf):								
	'''P45'''
	c=np.zeros((len(times),5),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(5,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-pardf['k4']*dt*c_temp[0]+g[i]*dt
			dc[1]=pardf['k0']*dt*c_temp[0]-pardf['k1']*dt*c_temp[1]
			dc[2]=pardf['k1']*dt*c_temp[1]-pardf['k2']*dt*c_temp[2]
			dc[3]=pardf['k2']*dt*c_temp[2]-pardf['k3']*dt*c_temp[3]
			dc[4]=pardf['k4']*dt*c_temp[0]+pardf['k3']*dt*c_temp[3]
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['A','B','C','D','Inf']									#this is optional but very useful. The species get names that represent some particular states
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c

def ABC_model(times,pardf):								
	'''Classical ABC model for solids, -A*n^1-2*B*n^2-3*C*n^3, A=k0, B=k1, C=k2, k3 single charge returns to excited state 
	note that the recombination looses 2 excitations (so B might be slightly different), C is auger recombination where the single charge has spectrum if a parameter with name [Auger] is present, which is invisible if not. Both Background and an non decaying (damaged) spectrum are implented and triggered by including of 'infinite' and 'background' as parameter '''
	c=np.zeros((len(times),3),dtype='float') 						#creation of matrix that will hold the concentrations
	g=gauss(times,sigma=pardf['resolution']/FWHM,mu=pardf['t0']) 	#creating the gaussian pulse that will "excite" our sample
	if 'sub_steps' in list(pardf.index.values):
		sub_steps=pardf['sub_steps']
	else:
		sub_steps=10  													#defining how many extra steps will be taken between the main time_points
	for i in range(1,len(times)):									#iterate over all timepoints
		dc=np.zeros(3,dtype='float')							#the initial change for each concentration, the "3" is representative of how many changes there will be
		dt=(times[i]-times[i-1])/(sub_steps)						# as we are taking smaller steps the time intervals need to be adapted
		c_temp=c[i-1,:]												#temporary matrix holding the changes (needed as we have sub steps and need to check for zero in the end)
		for j in range(int(sub_steps)):
			dc[0]=-pardf['k0']*dt*c_temp[0]-2*pardf['k1']*dt*c_temp[0]^2-3*pardf['k2']*dt*c_temp[0]^3+pardf['k3']*dt*c_temp[1]+g[i]*dt
			dc[1]=pardf['k2']*dt*c_temp[0]-pardf['k3']*dt*c_temp[1]
			dc[2]=pardf['k0']*dt*c_temp[0]+2*pardf['k1']*dt*c_temp[0]^2+2*pardf['k2']*dt*c_temp[0]^3
			c_temp=c_temp+dc
			c_temp[c_temp<0]=0			
			#for b in range(c.shape[1]):
			#	c_temp[b] =np.nanmax([(c_temp[b]+float(dc[b])),0.])		#check that nothing will be below 0 (concentrations)
		c[i,:] =c_temp												#store the temporary concentrations into the main matrix
	c=pandas.DataFrame(c,index=times)								#write back the right indexes
	c.index.name='time'												#and give it a name
	c.columns=['Excited_state','Auger','Inf']									#this is optional but very useful. The species get names that represent some particular states										#background always there (flat)
	if not 'auger' in list(pardf.index.values):
		c.drop('Auger',axis=1,inplace=True)
	if not 'infinite' in list(pardf.index.values):
		c.drop('Inf',axis=1,inplace=True)
	if 'explicit_GS' in list(pardf.index.values):
		c['GS']=-c.sum(axis=1)
	if 'background' in list(pardf.index.values):					#optional but usefull, allow the keyword "background" to be used to fit the background in the global analysis
		c['background']=1											#background always there (flat)
	return c
