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

from scipy import interpolate as scint
from numpy import arange,loadtxt,unique,where
from scipy.integrate import simps

def _extrap_helper(wave1,trans1,interpFunc1,wave2,interpFunc2):
        global idx,idx2
        idx,val=_find_nearest(w,int(math.ceil(wave1[0]/wavestep))*wavestep)
        idx2,val2=_find_nearest(w,int(math.floor(wave1[-1]/wavestep))*wavestep)
        interp=interpFunc1(arange(val,val2+1,wavestep))
        area=simps(f[idx:idx2+1]*interp,w[idx:idx2+1])#area in the band we have for color calculation
        val=int(math.ceil(wave2[0]/wavestep))*wavestep
        val2=int(math.floor(wave2[-1]/wavestep))*wavestep
        wave=arange(val,val2+1,wavestep)
        trans=interpFunc2(wave)
        print(area)# / sncosmo.constants.HC_ERG_AA)
        
        return(wave,trans,area)
def _find_nearest(array,value):
    """
    Helper function to find the nearest value in an array
    :param array: The array to search
    :param value: The value you want to match
    """

    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1,array[idx-1]
    else:
        return idx,array[idx]

def _getsed(sedfile) : 
    """
    Helper function that reads an SED file into a useable format.
    :param sedfile: SED filename to read (string)
    """

    d,w,f = loadtxt( sedfile, unpack=True,comments='#') 
    
    days = unique( d ) 

    dlist = [ d[ where( d == day ) ] for day in days ]
    wlist = [ w[ where( d == day ) ] for day in days ]
    flist = [ f[ where( d == day ) ] for day in days ]

    return( dlist, wlist, flist )

def createSncosmoSED(filename):
    #phase1, wave1, flux1 = _getsed(filename) #here the flux should be in ergs/s or ergs/A not sure
    #flux1=np.array(flux1[14])/(sncosmo.constants.HC_ERG_AA/np.array(wave1[14])) #now should be in photons
    '''
    global w
    w=np.array(wave1[14])
    global wavestep
    wavestep=wave1[14][1]-wave1[14][0]
    global f
    f=np.array(flux1[14])
    
    bWave=sncosmo.get_bandpass('bessellv').wave
    bTrans=sncosmo.get_bandpass('bessellv').trans
    bInterpFunc=scint.interp1d(bWave,bTrans)
    rWave=sncosmo.get_bandpass('sdss::r').wave
    rTrans=sncosmo.get_bandpass('sdss::r').trans
    rInterpFunc=scint.interp1d(rWave,rTrans)
    rWave,rTrans,bArea=_extrap_helper(bWave,bTrans,bInterpFunc,rWave,rInterpFunc)
    
    
    #flux1=np.array(flux1)/(sncosmo.constants.HC_ERG_AA/np.array(wave1))
    '''
    phase1,wave1,flux1=sncosmo.read_griddata_ascii(filename)
    
    source = ComboSource(phase1, wave1, flux1, name=filename[:-4])

    return(source)

def plotSncosmoSED(sed,phase=0.0,minwave=1500.0,maxwave=20000.0,wavestep=500,amplitude=1.0,showfig=True,savefig=False,outFilename="mySED.pdf",fontsize=18):
    
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm
    wave=np.linspace(minwave,maxwave,wavestep)
    '''
    cmap = ListedColormap(['k', 'b', 'k','gold','darkorange','r','k'])
    norm = BoundaryNorm([0, sncosmo.get_bandpass('bessellux').wave[0], sncosmo.get_bandpass('bessellux').wave[-1],
        sncosmo.get_bandpass('paritel::j').wave[0], sncosmo.get_bandpass('paritel::j').wave[-1],
        sncosmo.get_bandpass('paritel::h').wave[-1],sncosmo.get_bandpass('paritel::ks').wave[-1],maxwave], cmap.N)

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be numlines x points per line x 2 (x and y)
    points = np.array([wave, sed.flux(float(phase),wave)]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create the line collection object, setting the colormapping parameters.
    # Have to set the actual values used for colormapping separately.
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(wave)
    lc.set_linewidth(3)
    '''
    fig = plt.figure()
    ax=fig.gca()
    #ax.add_collection(lc)
    sed.set(amplitude=amplitude)
    
    ax.plot(wave, sed.flux(float(phase), wave),'k',linewidth=3)
    
    ax.fill_between(wave,sed.flux(float(phase), wave),
        where=[x>sncosmo.get_bandpass('bessellux').wave[0] and x < sncosmo.get_bandpass('bessellux').wave[-1] for x in wave],
        color='b', alpha=.3,edgecolor='k')
    ax.fill_between(wave,sed.flux(float(phase), wave),
        where=[x>sncosmo.get_bandpass('paritel::j').wave[0] and x < sncosmo.get_bandpass('paritel::j').wave[-1] for x in wave],
        color='gold', alpha=.3,edgecolor='k')
    ax.fill_between(wave,sed.flux(float(phase), wave),
        where=[x>sncosmo.get_bandpass('paritel::j').wave[-1] and x < sncosmo.get_bandpass('paritel::h').wave[-1] for x in wave],
        color='coral', alpha=.3)
    ax.fill_between(wave,sed.flux(float(phase), wave),
        where=[x>sncosmo.get_bandpass('paritel::h').wave[-1] and x < sncosmo.get_bandpass('paritel::ks').wave[-1] for x in wave],
        color='r', alpha=.3,edgecolor='k')
    
    ax.plot(np.append([0],wave),np.zeros(len(wave)+1),'--',color='black')
    #idx,val=_find_nearest(w,int(math.ceil(sncosmo.get_bandpass('paritel::j').wave[0]/wavestep))*wavestep)
    #idx2,val2=_find_nearest(w,int(math.ceil(sncosmo.get_bandpass('paritel::j').wave[-1]/wavestep))*wavestep)
    #ax.fill(w[idx:idx2+1],f[idx:idx2+1])
    if phase != 0:
        ax.set_title('Type II SED Extrapolation, %i Days from Peak'%phase,fontsize=18)
    else:
        ax.set_title('Type II SED Extrapolation, at Peak',fontsize=18)
    plt.xlabel("Wavelength (Angstrom)",fontsize=fontsize)
    plt.ylabel("Flux",fontsize=fontsize)
    if savefig:
        plt.savefig(os.path.basename(sed.name+".pdf"),format='pdf')
    if showfig:
        plt.show()
    return(fig)
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

#filename='SNLS-04D1la.SED'
#filename='SDSS-018892.SED'
#['SDSS-013195.SED','SDSS-014475.SED','SDSS-015475.SED','SDSS-017548.SED']
def main():
    files=glob.glob('*.SED')
    #filename='SDSS-000020.SED'
    #filename='../Nugent+Scolnic_IIL.SED'
    wave = np.linspace(1500.0, 9000.0, 500)
    days=0
    w=1.0
    fig = plt.figure()
    axes=[]
    ax = fig.add_subplot(111,frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("Wavelength (Angstrom)",fontsize=18)
    plt.ylabel("Flux",fontsize=18)
    ax.set_title('Type II SED Extrapolation, Days from Peak=%i'%days,fontsize=18)
    i=1
    for filename in files:
        print(filename)
        phase1, wave1, flux1 = sncosmo.read_griddata_ascii(filename)

        # In our __init__ method we defined above, the two fluxes need to be on
        # the same grid, so interpolate the second onto the first:

        source = ComboSource(phase1, wave1, flux1, name='extrapolated')

        ##########################################################################
        # We can get a summary of the Source we created:


        ##########################################################################
        # Get a spectrum at phase 10 for different parameters:
        source.set(amplitude=w)
        axes.append(fig.add_subplot(math.ceil(float(len(files))/3),3,i))
        axes[i-1].plot(wave, source.flux(float(days), wave))
        axes[i-1].text(.45*np.max(wave),.8*np.max(source.flux(float(days),wave)),filename,fontsize=6)
        #axes[i-1].locator_params(axis='y', nticks=4)
        #axes[i-1].locator_params(axis='x', nticks=6)
        axes[i-1].tick_params(labelcolor='none', top='off', right='off')

        #axes[i-1].set_xticks(np.linspace(axes[i-1].get_xticks()[0],axes[i-1].get_xticks()[-1],len(axes[i-1].get_xticks())/2))
        #axes[i-1].set_yticks(np.linspace(axes[i-1].get_yticks()[0],axes[i-1].get_yticks()[-1],len(axes[i-1].get_yticks())/2))

        i+=1

        
        
        #ax.plot(wave, source.flux(0., wave), label='w={:3.1f}'.format(w))

    plt.savefig("typeII2.pdf",format='pdf')
    #plt.show()
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
if __name__ == "__main__": main()