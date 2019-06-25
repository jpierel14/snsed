import glob,os,sys,scipy,sncosmo,math
import numpy as np
from scipy.interpolate import interp1d,interp2d
from astropy.table import Table,hstack
from astropy.io import ascii
from scipy.interpolate import UnivariateSpline as spl
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from .utils import __dir__,_defaultSpec
from .helpers import _find_nearest

'''
sne=np.loadtxt(os.path.join('data','spectra','snlist'),dtype='str')
files=glob.glob(os.path.join('data','spectra','*.lnw'))
types=[]
for sn in sne:
	with open(os.path.join('data','spectra',sn+'.lnw')) as theFile:
		for line in theFile:
			line=line.split()
			if ('-' in line[-3] and line[-3][-4:]=='norm'):
				types.append((sn,line[-3][:line[-3].rfind('-')]))
			elif '-' not in line[-3]:
				types.append((sn,line[-3]))
			break
np.savetxt(os.path.join('data','spectra','snTypes.ref'),types,fmt='%s')
'''

def createAveSpec(snType,specList,phases,waves):
	finalSpec=Table()
	minP=np.min(phases)
	maxP=np.max(phases)
	#finalSpec['wave']=waves
	#for p in phases:
	#	finalSpec[str(p)]=np.zeros(len(waves))
	allSpec=None
	for spec in specList:
		with open(spec) as f:
			for line in f:
				line=line.split()
				break
		data=np.loadtxt(spec,skiprows=int(line[4])+2)
		time=data[0][1:]
		if len(np.unique(time))!=len(time):
			print(os.path.basename(spec)+" has duplicate columns, skipping...")
			continue
		data=data[1:][:]
		tempSpec=Table(data,names=['wave']+[str(x) for x in time])

		if not allSpec:
			allSpec=Table()
			allSpec['wave']=tempSpec['wave']
		else:
			for col in tempSpec.colnames:
				if col != 'wave':
					if float(col)<=maxP:
						if float(col)>=minP:
							if col not in allSpec.colnames:
								allSpec[col]=tempSpec[col]
							else:
								allSpec[col]=np.average([allSpec[col],tempSpec[col]],axis=0)
					else:
						break
	
					
	x=np.sort([float(p) for p in allSpec.colnames if p !='wave'])
	
	y=allSpec['wave']
	
	
	z=[]
	for i in range(len(waves[waves<y[0]])):
		z.append(np.zeros(len(phases[phases<x[0]])+len(phases[phases>x[-1]])+len(allSpec.colnames)-1))

	for row in allSpec:
		z.append(np.append([0 for i in range(len(phases[phases<x[0]]))],np.append([row[str(p)] for p in x],[0 for i in range(len(phases[phases>x[-1]]))]))) 
	for i in range(len(waves[waves>y[-1]])):
		z.append(np.zeros(len(phases[phases<x[0]])+len(phases[phases>x[-1]])+len(allSpec.colnames)-1))

	x=np.append(phases[phases<x[0]],np.append(x,phases[phases>x[-1]]))
	y=np.append(waves[waves<y[0]],np.append(y,waves[waves>y[-1]]))
	func=interp2d(x,y,z)
	finalArr=func(phases,waves)
	finalSpec=Table(finalArr,names=['p:'+str(p) for p in phases])
	tempTable=Table()
	tempTable['wave']=waves
	finalSpec=hstack([tempTable,finalSpec])
	#ascii.write(finalSpec,os.path.join(__dir__,'data','default','type'+snType,snType+'_aveSpec.dat'),overwrite=True)
	#fig=plt.figure()
	#ax=fig.gca()
	#ax.plot(finalSpec['wave'],finalSpec['0'])
	#plt.show()
	return(finalSpec)

def _readSpec(spec):
	with open(spec) as f:
		for line in f:
			line=line.split()

			num=int(line[4])+2
			data=np.loadtxt(spec,skiprows=num)
			time=data[0][1:]
			if len(np.unique(time))!=len(time):
				print(os.path.basename(spec)+" has duplicate columns, skipping...")
				sys.exit()
			data=data[1:][:]
			allSpec=Table(data,names=['wave']+['p:'+str(x) for x in time])
			break
	return(allSpec)

def _chisquared(obs,actual):
	return(np.sum(((obs-actual)**2)))

def _specChi(const,sedFlux,tempFlux):
	return(_chisquared(sedFlux,tempFlux*const))

def _findBest(snType,wave,flux):
	flux=flux[wave>4200][wave[wave>4200]<6500]
	wave=wave[wave>4200][wave[wave>4200]<6500]

	mySpline=spl(wave,flux,k=5,ext=1)
	tempSpline=mySpline(wave)
	splineRemoved=np.log10(flux/tempSpline)
	sne,types=np.loadtxt(os.path.join(__dir__,'data','spectra','snTypes.ref'),dtype='str',unpack=True)
	mySpecs=sne[types==snType]
	bestChi=np.inf
	bestSpec=None
	bestConst=None
	for spec in mySpecs:
		data=_readSpec(os.path.join(__dir__,'data','spectra',spec+'.lnw'))
		x=np.array([float(col[2:]) for col in data.colnames if col !='wave'])


		if len(x)<2:
			continue
		y=data['wave']

		data.remove_column('wave')
		func=interp2d(x,y,np.array([list(r) for r in np.array(data)]))
		tempFlux=np.transpose(func(0,wave))[0]
		res=minimize(_specChi,np.array([1]),args=(splineRemoved,tempFlux))
		if res.fun<bestChi:
			bestChi=res.fun
			bestSpec=spec
			bestConst=res.x
	return(_readSpec(os.path.join(__dir__,'data','spectra',bestSpec+'.lnw')),bestConst)


def _addCCSpec(snType,sedFile,oldWave,specList=None,specName=None):
	phase,newWave,flux=sncosmo.read_griddata_ascii(sedFile)
	wave=np.append(newWave[newWave<oldWave[0]],newWave[newWave>oldWave[-1]])
	if len(wave)==0:
		return(flux)
	#flux=flux[:][newWave==wave]
	if specList:
		spec=createAveSpec(specList)
	elif specName:
		spec=os.path.join(__dir__,'data','spectra',specName)
		allSpec=_readSpec(spec)
	else:
		tempInd,tempVal=_find_nearest(phase,0)
		if snType=='II':
			snType='IIP'
		allSpec,bestConst=_findBest(snType,newWave,flux[tempInd])
		#spec=os.path.join(__dir__,'data','spectra',_defaultSpec[snType])
		#allSpec=_readSpec(spec)



	x=[float(x[2:]) for x in allSpec.colnames if x !='wave']
	y=allSpec['wave']
	allSpec.remove_column('wave')

	func=interp2d(x,y,np.array([list(r) for r in np.array(allSpec)]))
	tempFlux=np.transpose(func(phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]],newWave[newWave>4000][newWave[newWave>4000]<7500]))
	finalFlux=np.transpose(func(phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]],wave[wave>=y[0]][wave[wave>=y[0]]<=y[-1]]))
	if np.max(np.max(finalFlux))==0:
		return
	#if os.path.basename(sedFile)=='SDSS-015339.SED':
	#	print(finalFlux)
	#	print(wave[wave>=y[0]][wave[wave>=y[0]]<=y[-1]])
	#	sys.exit()

	splines=[]
	for p in phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]]:

		ind=np.where(phase==p)[0][0]
		mySpline=spl(newWave[newWave>4000][newWave[newWave>4000]<7500],flux[ind][newWave<7500][newWave[newWave<7500]>4000],k=5,ext=1)
		tempSpline=mySpline(newWave[newWave<7500][newWave[newWave<7500]>4000])
		splines.append(np.log10(flux[ind][newWave<7500][newWave[newWave<7500]>4000]/tempSpline))

	splines=np.array(splines)
	waves=wave[wave>=y[0]][wave[wave>=y[0]]<=y[-1]]
	for i in range(len(phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]])):
		#const1=math.fabs(np.max(splines[i])/np.max(tempFlux[i]))

		#const2=math.fabs(np.min(splines[i])/np.max(tempFlux[i]))
		const=minimize(_specChi,np.array([bestConst]),args=(splines[i],tempFlux[i])).x
		#const=np.nanmedian([math.fabs(k) for k in splines[i]/tempFlux[i]])
		finalFlux[i]*=const
		if i==4:
			final=const
		#if tempFlux[i][0]>splines[i][0]:
		#	constMinus=tempFlux[i][0]-splines[i][0]
		#	finalFlux[i]-=constMinus
		#else:
		#	constMinus=splines[i][0]-tempFlux[i][0]
		#	tempFlux[i]+=constMinus
		#splines[i]+=finalFlux[i]
		ind=np.where(phase==phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]][i])[0][0]
		for j in range(len(waves)):
			ind2=np.where(newWave==waves[j])[0][0]
			flux[ind][ind2]=max(0,flux[ind][ind2]+finalFlux[i][j])
	#fig=plt.figure()
	#ax=fig.gca()
	#ax.plot(newWave[newWave<7500][newWave[newWave<7500]>4000],tempFlux[4]*const)
	#ax.plot(newWave[newWave<7500][newWave[newWave<7500]>4000],splines[4])
	#plt.savefig('test_'+os.path.basename(sedFile[:-3])+'pdf',format='pdf',overwrite=True)
	#plt.close()

	#fig=plt.figure()
	#ax=fig.gca()
	#ax.plot(newWave[newWave<7500][newWave[newWave<7500]>4000],tempFlux[4]*final)
	#ax.plot(newWave[newWave<7500][newWave[newWave<7500]>4000],splines[4])
	#plt.savefig('test_'+os.path.basename(sedFile[:-3])+'pdf',format='pdf',overwrite=True)
	#plt.close()
	sncosmo.write_griddata_ascii(phase,newWave,flux,sedFile)
	return

def addCCSpec(snType,sed=None,waveRange=None,phase=None,wave=None,flux=None,specList=None,specName=None):
	if not sed:
		if not phase and not wave and not flux:
			print('Need SED file or phase,wave,flux arrays (found none).')
			sys.exit()
	else:
		phase,wave,flux=sncosmo.read_griddata_ascii(sed)




def addIaSpec():
	pass
#def main():
	#sne,types=np.loadtxt(os.path.join('data','spectra','snTypes.ref'),dtype='str',unpack=True)
	#mySne=sne[types=='IIP']
	#aveSpec=createAveSpec('II',[os.path.join('data','spectra',x+'.lnw') for x in mySne],np.arange(-50,150,1),np.arange(1000,20000,10))
#	#addSpec('Ia',None)
#	f=addSpec('Ib',np.arange(-5,21,1),np.arange(4000,17500,50),np.arange(4000,10000,50),np.zeros([len(np.arange(-5,21,1)),len(np.arange(4000,17500,50))]))


#if __name__=='__main__':
#	main()