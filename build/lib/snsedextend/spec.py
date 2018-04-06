import glob,os,sys,scipy,sncosmo
import numpy as np
from scipy.interpolate import interp1d,interp2d
from astropy.table import Table,hstack
from astropy.io import ascii
from scipy.interpolate import UnivariateSpline as spl

from .utils import __dir__,_defaultSpec
from helpers import _find_nearest

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
		spec=os.path.join(__dir__,'data','spectra',_defaultSpec[snType])
		allSpec=_readSpec(spec)



	x=[float(x[2:]) for x in allSpec.colnames if x !='wave']
	y=allSpec['wave']
	allSpec.remove_column('wave')


	func=interp2d(x,y,np.array([list(r) for r in np.array(allSpec)]))
	tempFlux=np.transpose(func(phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]],newWave[newWave>3500][newWave[newWave>3500]<9000]))
	finalFlux=np.transpose(func(phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]],wave[wave>=y[0]][wave[wave>=y[0]]<=y[-1]]))
	splines=[]
	for p in phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]]:

		ind=np.where(phase==p)[0][0]
		mySpline=spl(newWave[newWave>3500][newWave[newWave>3500]<9000],flux[ind][newWave<9000][newWave[newWave<9000]>3500],k=5,ext=1)
		tempSpline=mySpline(newWave[newWave<9000][newWave[newWave<9000]>3500])
		splines.append(np.log10(flux[ind][newWave<9000][newWave[newWave<9000]>3500]/tempSpline))

	splines=np.array(splines)
	waves=wave[wave>=y[0]][wave[wave>=y[0]]<=y[-1]]
	for i in range(len(phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]])):
		const1=np.median(tempFlux[i])/np.max(tempFlux[i])
		const2=np.median(splines[i])/np.max(splines[i])
		const=const1/const2
		finalFlux[i]*=const
		#if tempFlux[i][0]>splines[i][0]:
		#	constMinus=tempFlux[i][0]-splines[i][0]
		#	finalFlux[i]-=constMinus
		#else:
		#	constMinus=splines[i][0]-tempFlux[i][0]
		#	tempFlux[i]+=constMinus
		#splines[i]+=finalFlux[i]
		ind=np.where(phase==phase[phase>=x[0]][phase[phase>=x[0]]<=x[-1]][i])[0][0]
		for j in range(len(waves)):
			ind2=np.where(newWave==waves[i])[0][0]
			flux[ind][ind2]+=finalFlux[i][j]
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
#	f=addSpec('Ib',np.arange(-5,21,1),np.arange(4000,18000,50),np.arange(4000,10000,50),np.zeros([len(np.arange(-5,21,1)),len(np.arange(4000,18000,50))]))


#if __name__=='__main__':
#	main()