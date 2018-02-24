from __future__ import print_function

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys,warnings,os
import matplotlib.cm as cm

warnings.simplefilter('ignore')

import numpy as np
from scipy.interpolate import RectBivariateSpline
import sncosmo,sys,math
import matplotlib.pyplot as plt

class ComboSource(sncosmo.Source):

	_param_names = ['amplitude']
	param_names_latex = ['A']   # used in plotting display

	def __init__(self, phase, wave, flux1, name=None, version=None):
		self.name = name
		self.version = version
		self._phase = phase
		self._wave = wave
		
		self._model_flux1 = RectBivariateSpline(phase, wave, flux1, kx=3, ky=3)
		self._parameters = np.array([1.0])  # initial parameters

	def _flux(self, phase, wave):
		amplitude = self._parameters
		return amplitude  * self._model_flux1(phase, wave)

from sncosmo.builtins import DATADIR
from matplotlib import pyplot as plt
import glob,os

files=glob.glob(os.path.join('/Users','jpierel','rodney','snsedextend','SEDs','typeII','*.SED'))#[0:2]

wave = np.linspace(2000.0, 20000.0, 500)
days=0
w=1.0

for f in ['j','h','ks']:
		wave,trans=np.loadtxt(os.path.join('/Users','jpierel','rodney','snsedextend','filters',f+'.dat'),unpack=True)
		wave*=10000
		sncosmo.registry.register(sncosmo.Bandpass(wave,trans,name='paritel::'+f),force=True)

	
def uv_weighted_avg(values):
	"""Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.

	samplemean : include the correction factor for the "sample"
	weighted mean (as opposed to the mean of the full population): (see
	http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)

	"""
	samplemean=True
	weights=1/(uvdata['error'][np.where(uvdata['mag']==values[0])[0][0]:np.where(uvdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (values[0])
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (average)

def getUVErrors(values):

	samplemean=True
	weights=1/(uvdata['error'][np.where(uvdata['mag']==values[0])[0][0]:np.where(uvdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (1/np.sqrt(weights[0]))
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (np.sqrt(variance)/np.sqrt(samplecorrection))


def ir_weighted_avg(values):
	"""Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.

	samplemean : include the correction factor for the "sample"
	weighted mean (as opposed to the mean of the full population): (see
	http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)

	"""
	samplemean=True
	weights=1/(irdata['error'][np.where(irdata['mag']==values[0])[0][0]:np.where(irdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (np.array(values[0]))
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (average)

def getIRErrors(values):
	samplemean=True 
	weights=1/(irdata['error'][np.where(irdata['mag']==values[0])[0][0]:np.where(irdata['mag']==values[-1])[0][0]+1])**2
	N = len(np.nonzero(weights)[0])
	if N == 1:
		return (np.array(1/np.sqrt(weights[0])))
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	if samplemean:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection = 1
	return (np.sqrt(variance)/np.sqrt(samplecorrection))


def uvweightedMean(x):
	top=0
	bot=0
	errors=uvdata['error'][np.where(uvdata['mag']==x[0])[0][0]:np.where(uvdata['mag']==x[-1])[0][0]+1]
	
	for i in range(len(x)):
		top=top+(x[i]/(errors[i])**2)
		bot=bot+(1/(errors[i])**2)
	return(np.array(top/bot))

def irweightedMean(x):
	top=0
	bot=0
	errors=irdata['error'][np.where(irdata['mag']==x[0])[0][0]:np.where(irdata['mag']==x[-1])[0][0]+1]
	for i in range(len(x)):
		top=top+(x[i]/(errors[i])**2)
		bot=bot+(1/(errors[i])**2)
	return(np.array(top/bot))

def getErrors(x):
	err=0
	N = len(np.nonzero(x)[0])
	for i in range(len(x)):
		err=err+(1/(x[i])**2)
	if N>1:
		samplecorrection = np.float((N-1))/np.float(N)
	else:
		samplecorrection=1
	return((1/(err)**.5)/np.sqrt(samplecorrection))
		#rr=
	#return(np.mean(x))

t='II'
uvcol={
	'II':['V-r'],
}

UVcolorOrder=uvcol[t]
colors=ascii.read(os.path.join('type'+t,'tables','all'+t+'Colors.dat'))

plotColors=cm.nipy_spectral(np.linspace(0,1,15))
sn=np.unique(colors['SN'])
sne={sn[i]:plotColors[i] for i in range(len(sn))}

allUVColors=[]

allUVerr=[]

uvcolors=[]

uvtime=[]

uvnames=[]

for row in colors:
	for col in UVcolorOrder:
		if row[col]:
			allUVColors.append(row[col])
			allUVerr.append(row[col[0]+col[-1]+'_err'])
			uvcolors.append(col)
			uvtime.append(row['time'])
			uvnames.append(row['SN'])
			break

uvdata=Table([uvtime,uvcolors,allUVColors,allUVerr,uvnames],names=['time','color','mag','error','SN'],masked=True)
bin=np.array(np.trunc(uvdata['time']/.001))
uvgrouped=uvdata.group_by(bin)

uvtime=uvgrouped['time'].groups.aggregate(np.mean)
#uvmag=uvgrouped['mag'].groups.aggregate(uv_weighted_avg)
#uvmagerr=uvgrouped['mag'].groups.aggregate(getUVErrors)
uvmag=uvgrouped['mag'].groups.aggregate(uvweightedMean)
uvmagerr=uvgrouped['error'].groups.aggregate(getErrors)
uvnames=uvgrouped['SN']
plotColors=np.array([sne[uvnames[i]] for i in range(len(uvnames))])

temp=Table([np.array(uvtime),np.array(uvmag),np.array(uvmagerr)],names=('time','mag','magerr'))
ascii.write(temp,os.path.join('type'+t,'tables','uv.dat'))
fig=plt.figure()
ax = fig.add_subplot(111,frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel("Time (Since Peak)",fontsize=18)
plt.ylabel("Color Magnitude (Vega)",fontsize=18)
ax.set_title('Type II V-r Color',fontsize=18)
axes=[]
plots=[]
i=1
for filename in files:
	phase1, wave1, flux1 = sncosmo.read_griddata_ascii(filename)
	source = ComboSource(phase1, wave1, flux1, name='extrapolated')
	source.set(amplitude=w)
	model = sncosmo.Model(source=source)
	model.set(z=0.0)
	thisColor=[]
	for time in uvtime:
		thisColor.append(model.bandmag('bessellv','vega',time)-model.bandmag('sdss::r','vega',time))
	axes.append(fig.add_subplot(math.ceil(float(len(files))/4),4,i))
	axes[i-1].scatter(uvtime, thisColor,marker='*')
	axes[i-1].text(.45*np.max(wave),.8*np.max(source.flux(float(days),wave)),filename,fontsize=6)
	axes[i-1].tick_params(labelcolor='black', top='off', right='off')
	#axes[i-1].set_xlabel('Color Magnitude (Vega)')
	#axes[i-1].set_ylabel('Time (Since Peak)')
	for j in range(len(plotColors)):
		axes[i-1].errorbar(np.array(uvtime)[j],np.array(uvmag)[j],yerr=np.array(uvmagerr)[j],fmt='x',color=plotColors[j],label=str(np.array(uvnames[j])))
	axes[i-1].invert_yaxis()
	axes[i-1].text(.45*np.max(uvtime),.8*np.min(uvmag),os.path.basename(filename)[:-3],fontsize=6)
	i+=1
#ax.legend(plots,labels,fontsize=8,numpoints=1,loc=3)
#plt.title('V-r Average Color, Type '+t+' (Bin Size=2 Days)',size=15)
#fig.text(0.5, 0.02, 'Time (Since Peak)', ha='center',size=20)
#fig.text(0.02, 0.5, 'Color Magnitude (Vega)', va='center', rotation='vertical',size=20)
fig.savefig(os.path.join("type"+t,'plots','UVAverageColor.pdf'),format='pdf')
#plt.show()

