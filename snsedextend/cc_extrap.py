#! /usr/bin/env python
#S.rodney & J.R. Pierel
# 2017.03.31

import os,sys,glob,sncosmo,snsedextend
from numpy import *
from pylab import * 
from scipy import interpolate as scint
from scipy import stats
from astropy.table import Table
from scipy.interpolate import RectBivariateSpline

from .utils import *
from .helpers import *
from .BIC import BICrun

__all__=['extendCC','extendIa','fitColorCurve']

#Automatically sets environment to your SNDATA_ROOT file (Assumes you've set this environment variable,otherwise will use the builtin version)
try:
    sndataroot = os.environ['SNDATA_ROOT']
except:
    os.environ['SNDATA_ROOT']=os.path.join(os.path.dirname(snsedextend),'SNDATA_ROOT')
    sndataroot=os.environ['SNDATA_ROOT']

class snSource(sncosmo.Source):

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


def _getZP(band,zpsys):
    """
    Helper function that calculates the zero-point of a given band
    :param band: filename of band transmission function to get zero-point of 
    """
    ab=sncosmo.get_magsystem(zpsys)
    return(ab.band_flux_to_mag(1,band))


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

def _getFilter(band):
    w,t=loadtxt(_filters[band],unpack=True,comments='#')
    return(w,t)

def plotsed( sedfile,day, normalize=False,MINWAVE=None,MAXWAVE=None,**kwarg): 
    """
    Simple plotting function to visualize the SED file.
    :param sedfile: SED filename to read (string)
    :param day: The day you want to plot, or 'all'
    :param normalize: Boolean to normalize the data
    :param MINWAVE: Allows user to plot in a specific window, left edge
    :param MAXWAVE: Allows user to plot in a specific window, right edge
    """

    dlist,wlist,flist = _getsed( sedfile ) 

    for i in range( len(wlist) ) : 
        if MINWAVE:
            idx1,val= _find_nearest(wlist[i],MINWAVE)
        else:
            idx1=None
        if MAXWAVE:
            idx2,val=_find_nearest(wlist[i],MAXWAVE)
        else:
            idx2=None

        thisday = dlist[i][0]

        if day!='all' : 
            if abs(thisday-day)>0.6 : continue
        if normalize : 
            plot( wlist[i], flist[i]/flist[i].max()+thisday, **kwarg )
            show()
        else : 
            plot( wlist[i][idx1:idx2], flist[i][idx1:idx2], label=str(thisday), **kwarg )
            show()

def _extrapolate_uv(band,rArea,color,xTrans,xWave,f,w,niter,log,index,bandDict,zpsys,accuracy=1e-9):
    """
    Algorithm for extrapolation of SED into the UV
    """
    wavestep=w[1]-w[0]
    idx,val=_find_nearest(w,xWave[-1])
    idx2,val2=_find_nearest(w,xWave[0])
    if val-xWave[-1]>wavestep: #then need to extrapolate between end of SED and left edge of band
        Nstep = len( arange(xWave[-1], val, wavestep ) )
        wextRed =  sorted( [ xWave[-1] + (j+1)*wavestep for j in range(Nstep) ] )
        fextRed = array( [ max(0,f[idx]) for wave in wextRed ] )#just a flat slope from the end of the SED to the edge of the band
        w2 = append(wextRed,w[idx:])
        f2 = append(fextRed,f[idx:])

    else: #if the SED extends to the band, then cut it off at the band for new extrapolation
        w2=w[idx:]
        f2=f[idx:]
    if idx > idx2:
        w1=w[:idx2]
        f1=f[:idx2]
    else:
        w1=[]
        f1=[]
    w=w2
    f=f2
    ms=sncosmo.get_magsystem(zpsys)
    #print(rArea)
    bArea= ms.band_mag_to_flux(rArea+color,bandDict[band])*sncosmo.constants.HC_ERG_AA #area in the band we are extrapolating into (mag to flux)
    #bArea=rArea*color
    x1=w[0]
    x2=xWave[0]
    interpFunc=scint.interp1d(append(x2,x1),append(0,f[0]))
    #area=simps(xTrans*interpFunc(xWave),dx=wavestep) #this is the area if the slope is such that the flux is 0 at the right edge of the band
    area=np.sum(interpFunc(xWave)*xTrans*xWave*gradient(xWave))
    i=1
    extra=False
    #print(color,bArea,area)
    if area>bArea: #then we need flux=0, but slope of line steeper
        i=-1
        while area>bArea:
            i+=1
            x2=xWave[i]
            interpFunc=scint.interp1d(append(x2,w[0]),append(0,f[0]))
            idx=list(xWave).index(x2)
            #print(len(xTrans[0:idx+1]),len(arange(w[-1],x2+wavestep/10,wavestep)))
            area=np.sum(interpFunc(arange(x2,w[0]+wavestep/10,wavestep))*xTrans[idx:]*arange(x2,w[0]+wavestep/10,wavestep)*gradient(arange(x2,w[0]+wavestep/10,wavestep)))
        i=1
        x3=x2-wavestep
        idx2=list(xWave).index(x2)
        idx3=idx2-1
        while(abs(area-bArea)/(area+bArea)>accuracy and i<=niter): #this loop runs if necessary accuracy wasn't reached because of wavelength binning, it changes the starting flux slightly
            if i==1:
                y3=f[0]/2
                while area<bArea:
                    y3*=2
                    interpFunc1=scint.interp1d(append(x2,w[0]),append(y3,f[0]))
                    interpFunc2=scint.interp1d(append(x3,x2),append(0,y3))
                    #area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
                    area1=np.sum(interpFunc1(arange(x2,w[0]+wavestep/10,wavestep))*xTrans[idx:]*arange(x2,w[0]+wavestep/10,wavestep)*gradient(arange(x2,w[0]+wavestep/10,wavestep)))
                    #print(len(arange(x2,x3+wavestep/10,wavestep)))
                    area2=np.sum(interpFunc2(arange(x3,x2+wavestep/10,wavestep))*xTrans[idx3:idx2+1]*arange(x3,x2+wavestep/10,wavestep)*gradient(arange(x3,x2+wavestep/10,wavestep)))
                    area=area1+area2
                y1=0
                y2=y3
            if area>bArea:
                y3=y2
            else:
                y1=y2
            y2=(y3+y1)/2
            interpFunc1=scint.interp1d(append(x2,w[0]),append(y2,f[0]))
            interpFunc2=scint.interp1d(append(x3,x2),append(0,y2))
            #area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
            area1=np.sum(interpFunc1(arange(x2,w[0]+wavestep/10,wavestep))*xTrans[idx:]*arange(x2,w[0]+wavestep/10,wavestep)*gradient(arange(x2,w[0]+wavestep/10,wavestep)))
            #print(len(arange(x2,x3+wavestep/10,wavestep)))
            area2=np.sum(interpFunc2(arange(x3,x2+wavestep/10,wavestep))*xTrans[idx3:idx2+1]*arange(x3,x2+wavestep/10,wavestep)*gradient(arange(x3,x2+wavestep/10,wavestep)))
            area=area1+area2
            i+=1
        y1=f[0]
        extra=True

    elif area<bArea:#then the flux at the left edge of the band must be greater than 0
        y1=0
        y3=max(f)#initial upper bound is max of SED
        y2=y3
        x2=xWave[0]
        tried=False
        while(abs(area-bArea)/(area+bArea)>accuracy and i<=niter):
            if area>bArea:
                y3=y2
            else:
                if not tried:
                    while area < bArea:#this will increase the upper bound until the area is greater than the necessary area
                        y3*=2
                        interpFunc=scint.interp1d(append(x2,w[0]),append(y3,f[0]))
                        #area=simps(xTrans*interpFunc(xWave),dx=wavestep)
                        area=np.sum(interpFunc(xWave)*xTrans*xWave*gradient(xWave))
                    tried=True
                else:
                    y1=y2
            y2=(y3+y1)/2

            interpFunc=scint.interp1d(append(x2,w[0]),append(y2,f[0]))
            #area=simps(xTrans*interpFunc(xWave),dx=wavestep)
            area=np.sum(interpFunc(xWave)*xTrans*xWave*gradient(xWave))
            i+=1
        y1=f[0]
    else:
        y1=f[0]
        y2=0
    if abs(area-bArea)>.001*bArea:
        log.write('WARNING: The parameters you chose led to an integration in the %s-Band of %e instead of %e for index %i \n'%(band,vArea/area,color,index))
    (a,b,rval,pval,stderr)=stats.linregress(append(x2,w[0]),append(y2,y1))
    if extra:

        Nstep3 = len( arange( x2,w[0],  wavestep ) )
        wextRed3 =  sorted( [ x2 + (j+1)*wavestep for j in range(Nstep3) ])
        fextRed3 = array( [ max( 0, a * wave + b ) for wave in wextRed3 ] )
        if x3>xWave[0]:
            Nstep1=len(arange(  xWave[0],x3,  wavestep ) )
            wextRed1 =  sorted( [ xWave[0] + (j+1)*wavestep for j in range(Nstep1) ])
            fextRed1 = array( [ 0 for wave in wextRed1 ] )
        else:
            wextRed1=[]
            fextRed1=[]

        (a,b,rval,pval,stderr)=stats.linregress(append(x3,x2),append(0,y2))
        Nstep2=len(arange( x3, x2,  wavestep ) )
        wextRed2 =  sorted( [ x3 + (j+1)*wavestep for j in range(Nstep2) ])
        fextRed2 = array( [ max( 0, a * wave + b ) for wave in wextRed2 ] )



        wextRed=np.append(wextRed1,np.append(wextRed2,wextRed3))
        fextRed=np.append(fextRed1,np.append(fextRed2,fextRed3))
    else:
        Nstep = len( arange( xWave[0],w[0],  wavestep ) )
        wextRed =  sorted( [ xWave[0] + (j)*wavestep for j in range(Nstep) ] )
        fextRed = array( [ max( 0, a * wave + b ) if wave >= xWave[0] else max(0,a*xWave[0]+b) for wave in wextRed ] )


    w = append(w1, append(wextRed,w))
    f = append(f1, append(fextRed,f))

    return(w,f)

def _createSED(filename,rescale=False):
    phase,wave,flux=sncosmo.read_griddata_ascii(filename)
    if rescale:
        for i in range( len(phase) ) :
            flux[i]=flux[i]/sncosmo.constants.HC_ERG_AA
    return(snSource(phase,wave,flux,name='extrapolated'))

def _getPartialBandMag(model,phase,wave,band,zpsys,wavestep):
    interpFunc=scint.interp1d(band.wave,band.trans)
    waves=arange(band.wave[0],wave[-1]+1,wavestep)
    transInterp=interpFunc(waves)
    sncosmo.registry.register(sncosmo.Bandpass(waves,transInterp,name='tempBand'),force=True)

    return(model.bandmag('tempBand',zpsys,phase)) #magnitude in the overlap of the band


def _checkBandOverlap(band1,band2):
    if (band1.wave[0]<band2.wave[0] and band1.wave[-1]>band2.wave[0]) or (band2.wave[0]<band1.wave[0] and band2.wave[-1]>band1.wave[0]):
        return True
    return False

def _extrapolate_ir(band,vArea,color,xTrans,xWave,d,f,w,niter,log,index,doneIR,bandDict,zpsys,model,bandsDone,accuracy=1e-12):
    """
    Algorithm for extrapolation of SED into the IR
    """
    wavestep=w[1]-w[0]
    idx,val=_find_nearest(w,xWave[0]) #closest wavelength in the sed to the start of the band
    idx2,val2=_find_nearest(w,xWave[-1]) #closest wavelength in the sed to the end of the band
    overlapMag=0
    temp1=xTrans[xWave<=w[-1]]
    temp=xWave[xWave<=w[-1]]
    if xWave[0]-val>wavestep: #then need to extrapolate between end of SED and left edge of band
        Nstep = len( arange( val, xWave[0],  wavestep ) )
        wextRed =  sorted( [ val + (j+1)*wavestep for j in range(Nstep) ] )
        fextRed = array( [ max(0,f[idx]) for wave in wextRed ] )#just a flat slope from the end of the SED to the edge of the band
        w2 = append(w[:idx+1], wextRed)
        f2 = append(f[:idx+1], fextRed)
    else: #if the SED extends to the band:
        if not any([_checkBandOverlap(bandDict[x],bandDict[band]) for x in bandsDone]): #  then cut it off at the band for new extrapolation
            w2=w[:idx+1]
            f2=f[:idx+1]
        else:#calculate the area for the overlap between bands, then subtract that from what is necessary and begin extrapolation at the right edge of the last extrapolation
            w2=w
            f2=f
            bandIdx,bandVal=_find_nearest(xWave,w[-1])#int(math.ceil(wave1[0]/wavestep))*wavestep
            tempIdx,tempVal=_find_nearest(w,xWave[0])
            if xWave[0]<w[-1]:
                if bandVal<w[-1]:
                    xTrans=xTrans[bandIdx+1:]
                    xWave=xWave[bandIdx+1:]
                else:
                    xTrans=xTrans[bandIdx:]
                    xWave=xWave[bandIdx:]


    w=w2
    f=f2
    ms=sncosmo.get_magsystem(zpsys)

    try:
        overlapArea=np.sum(f[tempIdx:]*temp*temp1*gradient(temp))
    except:
        overlapArea=0

    bArea= ms.band_mag_to_flux(vArea-color,bandDict[band])*sncosmo.constants.HC_ERG_AA-overlapArea #area in the band we are extrapolating into (mag to flux), everything should be in ergs at this point
    #bArea=vArea/color
    x2=xWave[-1]
    x1=w[-1]
    interpFunc=scint.interp1d(append(x1,x2),append(f[-1],0))
    #interpFunc2=scint.interp1d(bandDict[band].wave,bandDict[band].dwave)
    #area=simps(xTrans*interpFunc(xWave),dx=wavestep) #this is the area if the slope is such that the flux is 0 at the right edge of the band
    area=np.sum(interpFunc(xWave)*xTrans*xWave*gradient(xWave))
    i=1
    #print(d[0],color)
    extra=False
    #print(bArea,area)
    if area>bArea: #then we need flux=0, but slope of line steeper

        #last=0
        i=0
        while area>bArea:
            i-=1
            x2=xWave[i]
            interpFunc=scint.interp1d(append(w[-1],x2),append(f[-1],0))
            idx=list(xWave).index(x2)
            #print(len(xTrans[0:idx+1]),len(arange(w[-1],x2+wavestep/10,wavestep)))
            area=np.sum(interpFunc(arange(w[-1],x2+wavestep/10,wavestep))*xTrans[:idx+1]*arange(w[-1],x2+wavestep/10,wavestep)*gradient(arange(w[-1],x2+wavestep/10,wavestep)))
        i=1
        x3=x2+wavestep
        idx2=list(xWave).index(x2)
        idx3=idx2+1
        while(abs(area-bArea)/(area+bArea)>accuracy and i<=niter): #this loop runs if necessary accuracy wasn't reached because of wavelength binning, it changes the starting flux slightly
            if i==1:
                y3=f[-1]/2
                while area<bArea:
                    y3*=2
                    interpFunc1=scint.interp1d(append(w[-1],x2),append(f[-1],y3))
                    interpFunc2=scint.interp1d(append(x2,x3),append(y3,0))
                    #area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
                    area1=np.sum(interpFunc1(arange(w[-1],x2+wavestep/10,wavestep))*xTrans[:idx2+1]*arange(w[-1],x2+wavestep/10,wavestep)*gradient(arange(w[-1],x2+wavestep/10,wavestep)))
                    #print(len(arange(x2,x3+wavestep/10,wavestep)))
                    area2=np.sum(interpFunc2(arange(x2,x3+wavestep/10,wavestep))*xTrans[idx2:idx3+1]*arange(x2,x3+wavestep/10,wavestep)*gradient(arange(x2,x3+wavestep/10,wavestep)))
                    area=area1+area2
                y1=0
                y2=y3
            if area>bArea:
                y3=y2
            else:
                y1=y2
            y2=(y3+y1)/2
            interpFunc1=scint.interp1d(append(w[-1],x2),append(f[-1],y2))
            interpFunc2=scint.interp1d(append(x2,x3),append(y3,0))
            #area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
            area1=np.sum(interpFunc1(arange(w[-1],x2+wavestep/10,wavestep))*xTrans[0:idx2+1]*arange(w[-1],x2+wavestep/10,wavestep)*gradient(arange(w[-1],x2+wavestep/10,wavestep)))
            area2=np.sum(interpFunc2(arange(x2,x3+wavestep/10,wavestep))*xTrans[idx2:idx3+1]*arange(x2,x3+wavestep/10,wavestep)*gradient(arange(x2,x3+wavestep/10,wavestep)))
            area=area1+area2
            i+=1
        y1=f[-1]
        extra=True

    elif area<bArea:#then the flux at the right edge of the band must be greater than 0
        y1=0
        y3=max(f)#initial upper bound is max of SED
        y2=y3
        x2=xWave[-1]
        tried=False
        while(abs(area-bArea)/(area+bArea)>accuracy and i<=niter):
            if area>bArea:
                y3=y2
            else:
                if not tried:
                    while area < bArea:#this will increase the upper bound until the area is greater than the necessary area
                        y3*=2
                        interpFunc=scint.interp1d(append(w[-1],x2),append(f[-1],y3))
                        #area=simps(xTrans*interpFunc(xWave),dx=wavestep)
                        area=np.sum(interpFunc(xWave)*xTrans*xWave*gradient(xWave))
                    tried=True
                else:
                    y1=y2
            y2=(y3+y1)/2

            interpFunc=scint.interp1d(append(w[-1],x2),append(f[-1],y2))
            #area=simps(xTrans*interpFunc(xWave),dx=wavestep)
            area=np.sum(interpFunc(xWave)*xTrans*xWave*gradient(xWave))
            i+=1
        y1=f[-1]
    else:
        y1=f[-1]
        y2=0
    #print(ms.band_flux_to_mag((area+overlapArea)/sncosmo.constants.HC_ERG_AA,bandDict[band]))
    if abs(area-bArea)>.001*bArea:
        log.write('WARNING: The parameters you chose led to an integration in the %s-Band of %e instead of %e for index %i \n'%(band,vArea-ms.band_flux_to_mag((area+overlapArea)/sncosmo.constants.HC_ERG_AA,bandDict[band]),color,index))
    (a,b,rval,pval,stderr)=stats.linregress(append(w[-1],x2),append(y1,y2))


    if extra:
        Nstep1 = len( arange( w[-1], x2,  wavestep ) )
        wextRed1 =  sorted( [ w[-1] + (j+1)*wavestep for j in range(Nstep1) ])
        fextRed1 = array( [ max( 0, a * wave + b ) for wave in wextRed1 ] )
        (a,b,rval,pval,stderr)=stats.linregress(append(x2,x3),append(y2,0))
        Nstep2=len(arange( x2, x3,  wavestep ) )
        wextRed2 =  sorted( [ x2 + (j+1)*wavestep for j in range(Nstep2) ])
        fextRed2 = array( [ max( 0, a * wave + b ) for wave in wextRed2 ] )
        if x3<xWave[-1]:
            Nstep3=len(arange( x3, xWave[-1],  wavestep ) )
            wextRed3 =  sorted( [ x3 + (j+1)*wavestep for j in range(Nstep3) ])
            fextRed3 = array( [ 0 for wave in wextRed3 ] )
        else:
            wextRed3=[]
            fextRed3=[]
        wextRed=np.append(wextRed1,np.append(wextRed2,wextRed3))
        fextRed=np.append(fextRed1,np.append(fextRed2,fextRed3))
    else:
        Nstep = len( arange( w[-1], xWave[-1],  wavestep ) )
        wextRed =  sorted( [ w[-1] + (j+1)*wavestep for j in range(Nstep) ] )

        fextRed = array( [ max( 0, a * wave + b ) if wave <= xWave[-1] else max(0,a*xWave[-1]+b) for wave in wextRed ] )

    w = append(w, wextRed)
    f = append(f, fextRed)
    return(w,f)



def _extrapolatesed(sedfile, newsedfile,color,table,time,modColor, bands,zpsys,bandsDone,niter=50):
    """ Interpolate the given transmission function to the wavestep of the SED, then
    get the area in the V band for color calculation and run the extrapolation algorithm.
    :param sedfile: SED file being extrapolated
    :param newsedfile: filename of output SED file
    :param fVH 
    """
    dlist,wlist,flist = _getsed( sedfile ) #first time this is read in, should be in ergs/s/cm^2/AA
    '''
    tempTime=[]
    tempColor=[]
    extrap=True
    
    if dlist[0][0]<table[_get_default_prop_name('time')][0]:
        tempTime=append(math.floor(dlist[0][0]),time)
        tempColor=append(table[color][0],modColor)
    else:
        extrap=False
    if dlist[-1][0]>time[-1]:
        if np.any(tempTime):
            tempTime=append(tempTime,math.ceil(dlist[-1][0]))
            tempColor=append(tempColor,modColor[-1])
        else:
            tempTime=append(tempTime,append(time,math.ceil(dlist[-1][0])))
            tempColor=append(tempColor,append(table[color][-1],modColor))
    elif dlist[-1][0]<=time[-1] and not extrap:
        extrap=False
    if not extrap:
        tempTime=time
        tempColor=modColor
    '''
    i=0

    while dlist[i][0]<time[0]:
        i+=1
    dlist=dlist[i:]
    wlist=wlist[i:]
    flist=flist[i:]
    i=-1
    while dlist[i][0]>time[-1]:
        i-=1
    if i!=-1:
        dlist=dlist[:i+1]
        wlist=wlist[:i+1]
        flist=flist[:i+1]


    blue=color[0]
    red=color[-1]
    bWave=bands[blue].wave
    bTrans=bands[blue].trans
    rWave=bands[red].wave
    rTrans=bands[red].trans
    bInterpFunc=scint.interp1d(bWave,bTrans)
    rInterpFunc=scint.interp1d(rWave,rTrans)
    cInterpFunc=scint.interp1d(time,modColor)
    tempTime=[x[0] for x in dlist]
    colorData=cInterpFunc(tempTime)

    sed=_createSED(sedfile,rescale=False) #now the original sed is in photons/s/cm^2
    model=sncosmo.Model(sed)

    fout = open( newsedfile, 'w' )
    log= open('./error.log','w')

    def _extrap_helper(wave1,interpFunc1,wave2,interpFunc2,known,currentPhase):
        #idx,val=_find_nearest(w,int(math.ceil(wave1[0]/wavestep))*wavestep)
        #idx2,val2=_find_nearest(w,int(math.floor(wave1[-1]/wavestep))*wavestep)
        #interp=interpFunc1(arange(val,val2+1,wavestep))
        #area=simps(f[idx:idx2+1]*interp,w[idx:idx2+1])#area in the band we have for color calculation
        area=model.bandmag(bands[known],zpsys,currentPhase)

        #area=bandMag(sedfile,currentPhase,bands[known],zpsys)#magnitude in the band we have for color calculation
        val=int(math.ceil(wave2[0]/wavestep))*wavestep
        val2=int(math.floor(wave2[-1]/wavestep))*wavestep
        wave=arange(val,val2+1,wavestep)
        trans=interpFunc2(wave)
        return(wave,trans,area)
    UV=False
    IR=False
    for i in range( len(dlist) ) :
        d,w,f = dlist[i],wlist[i],flist[i]
        #if not bandsDone:
            #f=f/sncosmo.constants.HC_ERG_AA #flux is now in photons/s/cm^2

        wavestep = w[1] - w[0]
        if bands[blue].wave_eff<=_UVrightBound:
            bWave,bTrans,rArea=_extrap_helper(rWave,rInterpFunc,bWave,bInterpFunc,red,d[0])
            wnew,fnew=_extrapolate_uv(blue,rArea,colorData[i],bTrans,bWave,f,w,niter,log,i,bands,zpsys)
            UV=True
        elif bands[red].wave_eff>=_IRleftBound:
            rWave,rTrans,bArea=_extrap_helper(bWave,bInterpFunc,rWave,rInterpFunc,blue,d[0])
            wnew,fnew=_extrapolate_ir(red,bArea,colorData[i],rTrans,rWave,d,f,w,niter,log,i,IR,bands,zpsys,model,bandsDone)
            IR=True
        else:
            raise RuntimeError("You supplied a color that does not support extrapolation to the IR or UV!")
        for j in range( len( wnew ) ) :
            fout.write("%5.1f  %10i  %12.7e \n"%( d[0], wnew[j], fnew[j] ))
    fout.close() 
    log.close()
    return( UV,IR )

def _boundUVsed(sedfile):
    dlist,wlist,flist = _getsed(sedfile)
    fout = open( sedfile, 'w' )
    for i in range(len(dlist)):
        d,w,f = dlist[i],wlist[i],flist[i]
        if w[0]>_UVleftBound:
            wavestep=w[1]-w[0]
            (a,b,rval,pval,stderr)=stats.linregress(append(_UVleftBound,w[0]),append(0,f[0]))
            Nstep = len( arange( _UVleftBound, w[0],  wavestep ) )
            wext =  sorted( [ _UVleftBound + (j)*wavestep for j in range(Nstep) ] )
            fext = array( [ max( 0, a * wave + b ) for wave in wext ] )
            wnew=append(wext,w)
            fnew=append(fext,f)
        else:
            wnew=w
            fnew=f
        for j in range( len( wnew ) ) :
            fout.write("%5.1f  %10i  %12.7e \n"%( d[0], wnew[j], fnew[j] ))
    fout.close()

def _boundIRsed(sedfile):
    dlist,wlist,flist = _getsed(sedfile)
    fout = open( sedfile, 'w' )
    for i in range(len(dlist)):
        d,w,f = dlist[i],wlist[i],flist[i]
        if w[-1]<_IRrightBound:
            wavestep=w[1]-w[0]
            (a,b,rval,pval,stderr)=stats.linregress(append(w[-1],_IRrightBound),append(f[-1],0))
            Nstep = len( arange( w[-1], _IRrightBound,  wavestep ) )
            wext =  sorted( [ w[-1] + (j+1)*wavestep for j in range(Nstep) ] )
            fext = array( [ max( 0, a * wave + b ) for wave in wext ] )
            wnew=append(w,wext)
            fnew=append(f,fext)
        else:
            wnew=w
            fnew=f
        for j in range( len( wnew ) ) :
            fout.write("%5.1f  %10i  %12.7e \n"%( d[0], wnew[j], fnew[j] ))
    fout.close()



def bandMag(filename,currentPhase,band,zpsys,rescale=False):
    sed=_createSED(filename,rescale=rescale)
    model=sncosmo.Model(sed)
    return(model.bandmag(band,zpsys,currentPhase))




def getExtremes(time,curve,original,color):
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

    error=getErrors(np.array(original[color[0]+color[-1]+'_err'])[np.array(original[color[0]+color[-1]+'_err'])!=0])

    blue=curve-9*error
    red=curve+9*error
    plt.plot(time,curve,color='k')
    plt.plot(time,blue,color='b')
    plt.plot(time,red,color='r')
    plt.show()
    sys.exit()

def extendIa():
    pass

def fitColorCurve(table,confidence=50,type='II',verbose=True):

    colors=[x for x in table.colnames if len(x)==3 and x[1]=='-']
    result=dict([])
    print("Running BIC for color:")
    for color in colors:
        print('     '+color)
        tempTable=table[~table[color].mask]
        tempTable=Table([tempTable['time'],tempTable[color],tempTable[color[0]+color[-1]+'_err']],names=('time','mag','magerr'))
        temp=BICrun(tempTable,type)
        result[color]=Table([temp['x'],temp[str(confidence*10)]],names=('time',color))
    return(result)

def extendCC(colorTable,colorCurveDict,outFileLoc='.',bandDict=_filters,colors=None,zpsys='AB',sedlist=None,showplots=None,verbose=True):
    """
    Function called in main which runs the extrapolation algorithm.
    :param sedlist: list of files to analyze (or none if all in current directory)
    :param fVJ: input V-J
    :param fVH: input V-H
    :param fVK: input V-K
    :param showplots: None if you don't want plotting, otherwise it's 'all' if you want to plot all days or a single day (i.e. 3)
    """
    #import pickle
    colorTable=_standardize(colorTable)
    #print('Running Bayesian Information Criterion...')
    #temp=snsedextend.BIC.run()
    #sys.exit()
    #with open('tempBIC.txt','wb') as handle:
    #    pickle.dump(temp,handle)
    ###### for debugging#####
    #with open('tempBIC.txt','rb') as handle:
    #    temp=pickle.loads(handle.read())


    #bandDict=dict((k.upper(),v) for k,v in bandDict.iteritems())
    if not isinstance(colorTable,Table):
        raise RuntimeError('Colors argument must be an astropy table.')

    for band in bandDict:
        if not isinstance(bandDict[band],sncosmo.Bandpass):
            bandDict=_bandCheck(bandDict,band)
    if not colors:
        print("No extrapolation colors defined, assuming: ",colorCurveDict.keys())
        colors=[x for x in colorCurveDict.keys() if 'U-' in x]+[x for x in colorCurveDict.keys() if '-J' in x]+[x for x in colorCurveDict.keys() if '-H' in x]+[x for x in colorCurveDict.keys() if '-K' in x]
    else:
        tempColors=[x for x in colors if 'U-' in x]+[x for x in colors if '-J' in x]+[x for x in colors if '-H' in x]+[x for x in colors if '-K' in x]
        if len(tempColors)>4:
            raise RuntimeError("Only can do 4 colors!")


    bands=append([col[0] for col in colors],[col[-1] for col in colors])
    for band in _filters.keys():
        if band not in bandDict.keys() and band in bands:
            bandDict[band]=sncosmo.get_bandpass(_filters[band])
    if not sedlist:
        sedlist = glob.glob(os.path.join(sndataroot,"snsed","NON1A","*.SED"))
    else:
        sedlist=[sedlist] if not isinstance(sedlist,(tuple,list)) else sedlist
        sedlist=[os.path.join(sndataroot,"snsed","NON1A",sed) for sed in sedlist]
    for sedfile in sedlist:
        newsedfile=os.path.join(outFileLoc,os.path.basename(sedfile))
        if verbose:
            print("EXTRAPOLATING %s"%newsedfile)
        once=False
        boundUV=False
        boundIR=False
        bandsDone=[]

        for color in colors:
            if verbose:
                print('     Running %s extrapolation.'%color)
            if color[0] not in bandDict.keys():
                raise RuntimeError('Band "%s" defined in color "%s", but band is not defined.')
            if color[-1] not in bandDict.keys():
                raise RuntimeError('Band "%s" defined in color "%s", but band is not defined.')


            #getExtremes(colorCurveDict[color]['time'],colorCurveDict[color][color],colorTable,color)

            tempMask=colorTable[color].mask

            colorTable[color].mask=tempMask
            if once:
                sedfile=newsedfile
            tempTable=colorTable[~colorTable[color].mask]
            UV,IR=_extrapolatesed(sedfile, newsedfile,color,tempTable,colorCurveDict[color]['time'],colorCurveDict[color][color], bandDict,zpsys,bandsDone,niter=50)
            if UV:
                boundUV=True
                bandsDone.append(color[0])
            if IR:
                boundIR=True
                bandsDone.append(color[-1])
            once=True


        if showplots:
            plotsed(newsedfile,day=showplots)
        if boundUV:
            _boundUVsed(newsedfile)
        if boundIR:
            _boundIRsed(newsedfile)
        if verbose:
            print("     Done with %s.\a\a\a"%os.path.basename(sedfile))



