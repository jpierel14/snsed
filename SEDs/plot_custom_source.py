"""
===========================
Creating a new Source class
===========================

Extending sncosmo with a custom type of Source.

A ``Source`` is something that specifies a spectral timeseries as
a function of an arbitrary number of parameters. For example, the SALT2
model has three parameters (``x0``, ``x1`` and ``c``) that determine a
unique spectrum as a function of phase. The ``SALT2Source`` class implements
the behavior of the model: how the spectral time series depends on those
parameters.

If you have a spectral timeseries model that follows the behavior of one of
the existing classes, such as ``TimeSeriesSource``, great! There's no need to
write a custom class. However, suppose you want to implement a model that
has some new parameterization. In this case, you need a new class that
implements the behavior.

In this example, we implement a new type of source model. Our model is a linear
combination of two spectral time series, with a parameter ``w`` that
determines the relative weight of the models.
"""

from __future__ import print_function

import numpy as np
from scipy.interpolate import RectBivariateSpline
import sncosmo,sys
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


########################################################################
# ... and that's all that we need to define!: A couple class attributes
# (``_param_names`` and ``param_names_latex``, an ``__init__`` method,
# and a ``_flux`` method. The ``_flux`` method is guaranteed to be passed
# numpy arrays for phase and wavelength.
#
# We can now initialize an instance of this source from two spectral time
# series:

#Just as an example, we'll use some undocumented functionality in
# sncosmo to download the Nugent Ia and 2p templates. Don't rely on this
# the `DATADIR` object, or these paths in your code though, as these are
# subject to change between version of sncosmo!
from sncosmo.builtins import DATADIR
#filename='SNLS-04D1la.SED'
#filename='SDSS-018892.SED'
#['SDSS-013195.SED','SDSS-014475.SED','SDSS-015475.SED','SDSS-017548.SED']
filename='SDSS-000020.SED'
filename2='../overflow_snsedextend/SDSS-003818.SED'
#filename='../Nugent+Scolnic_IIL.SED'
phase1, wave1, flux1 = sncosmo.read_griddata_ascii(filename)
phase2, wave2, flux2 = sncosmo.read_griddata_ascii(filename2)

# In our __init__ method we defined above, the two fluxes need to be on
# the same grid, so interpolate the second onto the first:

source = ComboSource(phase1, wave1, flux1, name='extrapolated')
source2=ComboSource(phase2, wave2, flux2, name='original')

##########################################################################
# We can get a summary of the Source we created:


##########################################################################
# Get a spectrum at phase 10 for different parameters:

from matplotlib import pyplot as plt

wave = np.linspace(2000.0, 20000.0, 500)
w=1.0
source.set(amplitude=w)
ax=plt.gca()
plt.plot(wave, source.flux(0., wave), label='w={:3.1f}'.format(w))
ax.set_title('SED: SDSS-013195 (Type Ic), UV Extrapolation, Days from Peak=10')
ax.set_xlabel('Wavelength (Angstrom)')
ax.set_ylabel('Flux')
plt.show()
sys.exit()
##########################################################################
# The w=0 spectrum is that of the Ia model, the w=1 spectrum is that of
# the IIp model, while intermediate spectra are weighted combinations.
#
# We can even fit the model to some data!
model = sncosmo.Model(source=source)
for f in ['j','h','ks']:
    wave,trans=np.loadtxt(os.path.join('/Users','jpierel','rodney','snsedextend','filters',f+'.dat'),unpack=True)
    wave*=10000
    sncosmo.registry.register(sncosmo.Bandpass(wave,trans,name='paritel::'+f),force=True)
model.set(z=0.021308074)
print(model.bandmag('bessellux','vega',0)-model.bandmag('bessellv','vega',0))
'''
wave,trans=np.loadtxt('snsedextend/kband/tophatK.dat',unpack=True)
band=sncosmo.Bandpass(wave,trans,name='tophatk')
sncosmo.registry.register(band)

model = sncosmo.Model(source=source)
model2= sncosmo.Model(source=source2)
Hvals=[]
Jvals=[]
Kvals=[]
for i in np.arange(0,1.01,.01):
    model.set(z=i)
    #Hvals.append(model.bandmag('bessellv','ab',0)-model.bandmag('f160w','ab',0))
    #Jvals.append(model.bandmag('bessellv','ab',0)-model.bandmag('cspjs','ab',0))
    Kvals.append(model.bandmag('bessellv','ab',0)-model.bandmag('tophatk','ab',0))
tgrid=np.arange(0,1.01,.01)
#plt.plot(tgrid,Hvals,'-',color='red')
fig=plt.figure()


ax=fig.add_subplot(111)

ax.invert_yaxis()
ax.set_title('V-K=0.95, Filter=Tophat')
ax.set_xlabel('Redshift')
ax.set_ylabel('Peak Magnitude (AB)')
ax.plot(tgrid,Kvals,'-',color='red')
plt.savefig('/Users/jpierel/Desktop/V-K')
'''
'''
tgrid=np.linspace(model.mintime(),model.maxtime()+1,int(model.maxtime()-model.mintime())+1)
mflux=model.bandflux('bessellv',tgrid)
plt.plot(tgrid,mflux)
plt.show()
'''
#print(model.bandmag('bessellv','ab',0)-model.bandmag('f160w','ab',0))

'''
data = sncosmo.read_lc('refsdalS1_psfphot.dat')
result, fitted_model = sncosmo.fit_lc(data, model,['t0', 'amplitude'])
result2, fitted_model2 = sncosmo.fit_lc(data, model2,['t0', 'amplitude'])

sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
plt.show()
'''
##########################################################################
# The fact that the fitted value of w is closer to 0 than 1 indicates that
# the light curve looks more like the Ia template than the IIp template.
# This is generally what we expected since the example data here was
# generated from a Ia template (although not the Nugent template!).
