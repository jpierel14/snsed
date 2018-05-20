import warnings
import numpy as np
import astropy.units as u
from astropy.modeling.blackbody import FLAM
from astropy.modeling.models import BlackBody1D
from scipy.optimize import minimize

from scipy.signal import savgol_filter as smooth


warnings.simplefilter('ignore')


def _chisq(obs,actual):
    return(np.sum((obs-actual)**2))

def _bbChi(stuff,wave,flux):
    temp,const=stuff
    bb=BlackBody1D(temperature=temp*u.K)
    bbFlux=[x.value for x in bb(wave).to(FLAM,u.spectral_density(wave))]

    bbFlux=[x*const for x in bbFlux]
    return _chisq(bbFlux,flux)


def fluxFromBB(temp,const,wave):
    bb=BlackBody1D(temperature=temp*u.K)
    bbFlux=[x.value for x in bb(wave*u.AA).to(FLAM,u.spectral_density(wave*u.AA))]
    bbFlux=np.array([x*const for x in bbFlux])
    return(bbFlux)

def getBB(phase,wave,flux):

    #fluxes=[]

    newWave=np.arange(5000,30000,50)
    allFlux=[]
    for p in range(len(phase)):
        if phase[p]>3 or phase[p]<-3:
            allFlux.append(np.zeros(len(newWave)))
        else:
            flux2=flux[p][wave<9000][wave[wave<9000]>4200]
            wave2=wave[wave<9000][wave[wave<9000]>4200]*u.AA
            res=minimize(_bbChi,np.array([6000,1]),args=(wave2,smooth(flux2,len(flux2)/3,2)),bounds=((0,None),(0,None)))

            temp,const=res.x
            bbFlux=fluxFromBB(temp,const,newWave)
            allFlux.append(bbFlux)
    return(newWave,np.array(allFlux))



