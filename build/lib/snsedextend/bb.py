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

def getBB(phase,wave,flux,name,startWave=4200,endWave=9000,winConst=6,constant=1,ax=None):

    newWave=np.arange(5000,30000,50)
    allFlux=[]
    for p in range(len(phase)):
        if phase[p]>3 or phase[p]<-3:
            allFlux.append(np.zeros(len(newWave)))
        else:
            flux2=flux[p][wave<endWave][wave[wave<endWave]>startWave]*constant
            wave2=wave[wave<endWave][wave[wave<endWave]>startWave]*u.AA
            win=int(len(flux2)/winConst) if int(len(flux2)/winConst)%2==1 else int(len(flux2)/winConst)+1
            res=minimize(_bbChi,np.array([6000,1]),args=(wave2,smooth(flux2,win,2)),bounds=((1000,20000),(1,None)))

            temp,const=res.x
            const/=constant
            bbFlux=fluxFromBB(temp,const,newWave)
            allFlux.append(bbFlux)
    return(newWave,np.array(allFlux))



