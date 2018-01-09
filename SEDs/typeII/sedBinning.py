from __future__ import print_function

import numpy as np
from scipy.interpolate import RectBivariateSpline
import sncosmo,sys,os
import matplotlib.pyplot as plt
from sncosmo.builtins import DATADIR
from matplotlib import pyplot as plt
import glob,math

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

def createSncosmoSED(filename):
    #phase1, wave1, flux1 = _getsed(filename) #here the flux should be in ergs/s or ergs/A not sure
    #flux1=np.array(flux1[14])/(sncosmo.constants.HC_ERG_AA/np.array(wave1[14])) #now should be in photons
    phase1,wave1,flux1=sncosmo.read_griddata_ascii(filename)
    
    source = ComboSource(phase1, wave1, flux1, name=filename[:-4])

    return(source)

files=glob.glob('*.SED')
fig=plt.figure()
ax=fig.gca()
med=[]
for filename in files:
	sed=createSncosmoSED(filename)
	model=sncosmo.Model(sed)
	time=[0]
	colors=[model.color('bessellv','sdss::r','vega',x) for x in time]
	#time=np.linspace(-20,50,100)
	med.append(colors[0])
med=np.median(med)
print(med)
for filename in files:
	sed=createSncosmoSED(filename)
	model=sncosmo.Model(sed)
	time=[0]
	colors=[model.color('bessellv','sdss::r','vega',x) for x in time]
	if colors[0] > med+.06:
		ax.scatter(time,colors,color='r')
	elif colors[0]<med-.06:
		ax.scatter(time,colors,color='b')
	else:
		ax.scatter(time,colors,color='k')
ax.fill_between([-5,5],[med+.06,med+.06],med-.06,
        where=[True,True],color='b', alpha=.3,edgecolor='k')
#ax.plot([0,30],[med,med],l)
plt.show()







