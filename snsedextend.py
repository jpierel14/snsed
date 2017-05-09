#! /usr/bin/env python
#S.rodney & J.R. Pierel
# 2017.03.31

import os,sys,getopt,warnings,math,glob,shutil
from numpy import *
from pylab import * 
from scipy import interpolate as scint
from scipy import stats
from scipy.integrate import simps
from copy import deepcopy
#Automatically sets environment to your SNDATA_ROOT file (Assumes you've set this environment variable,otherwise will use the builtin version)
try:
    sndataroot = os.environ['SNDATA_ROOT']
except:
    os.environ['SNDATA_ROOT']=os.path.join(os.path.dirname(snsedextend),'SNDATA_ROOT')
    sndataroot=os.environ['SNDATA_ROOT']

#default transmission files, user can define their own 
_filters={
    'V':'bessellv',
    'J':'tophatJ',
    'H':'tophatH',
    'K':'tophatK'
}

#dictionary of zero-points (calculated later based on transmission files)
_zp={
    'AB':{
        'V':None,
        'J':None,
        'H':None,
        'K':None
        },
    'Vega':{
        'V':-0.02,
        'J':-0.91,
        'H':-1.39,
        'K':-1.85   
        }
}

#just a helpter dictionary
_colToFilter={
    'V-J':'J',
    'V-H':'H',
    'V-K':'K'
}

def _findFile(filename):
    """Helper function that does a quick recurisive directory seach for a transmission file.
    :param filename: Transmission file to find
    """
    for root, subFolders, files in os.walk('.'):
        if filename in files:
            return(os.path.abspath(os.path.join(root,filename)))

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

def _getZP(band):
    """
    Helper function that calculates the zero-point of a given band
    :param band: filename of band transmission function to get zero-point of 
    """
    c=299792458
    wave,trans=_getFilter(band)
    fInterp=scint.interp1d(wave,trans)
    if wave[1]-wave[0]>10:
        x=np.arange(wave[0],wave[-1]+1,10)
    else:
        x=wave
    y=fInterp(x)
    xnu=c/(x*1E-10)

    arr=xnu.argsort()
    y2=y[arr[0::]]
    x2=np.sort(xnu)
    F=simps(y2*3.63E-20,x2)
    return(2.5*log10(F))

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


def extrapolate_band(band,vArea,color,xTrans,xWave,f,w,niter,log,index,accuracy=1e-9):
    """
    Algorithm for extrapolation of SED into the IR
    """
    wavestep=w[1]-w[0]
    idx,val=_find_nearest(w,xWave[0])
    if xWave[0]-val>wavestep: #then need to extrapolate between end of SED and left edge of band
        Nredstep = len( arange( val, xWave[0],  wavestep ) )
        wextRed =  sorted( [ val + (j+1)*wavestep for j in range(Nredstep) ] )
        fextRed = array( [ max(0,f[idx]) for wave in wextRed ] )#just a flat slope from the end of the SED to the edge of the band
        w = append(w[:idx+1], wextRed)
        f = append(f[:idx+1], fextRed)

    else: #if the SED extends to the band, then cut it off at the band for new extrapolation
        w=w[:idx+1]
        f=f[:idx+1]

    bArea=vArea/color
    x2=xWave[-1]
    x1=w[-1]
    interpFunc=scint.interp1d(append(x1,x2),append(f[-1],0))
    area=simps(xTrans*interpFunc(xWave),dx=wavestep) #this is the area if the slope is such that the flux is 0 at the right edge of the band
    i=1
    if area>bArea: #then we need flux=0, but slope of line steeper
        last=0
        while(abs(area-bArea)/(area+bArea)>accuracy and i<=niter):
            if area>bArea:
                x3=x2
            else:
                x1=x2
            x2=(x3+x1)/2
            idx,x2=_find_nearest(xWave,x2)
            interpFunc=scint.interp1d(append(w[-1],x2),append(f[-1],0))
            area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
            i+=1
            if area==last: #break out of loop once we can't get more accurate because of wavelength binning
                break
            last=area

        i=1
        while(abs(area-bArea)/(area+bArea)>accuracy and i<=niter): #this loop runs if necessary accuracy wasn't reached because of wavelength binning, it changes the starting flux slightly
            if i==1:
                y3=f[-1]
                while area<bArea:
                    y3*=2
                    interpFunc=scint.interp1d(append(w[-1],x2),append(y3,0))
                    area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
                y1=f[-1]
                y2=y3
            if area>bArea:
                y3=y2
            else:
                y1=y2
            y2=(y3+y1)/2
            interpFunc=scint.interp1d(append(w[-1],x2),append(y2,0))
            area=simps(xTrans[0:idx+1]*interpFunc(arange(w[-1],x2+wavestep/10,wavestep)),dx=wavestep)
            i+=1
        y1=y2
        y2=0
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
                        area=simps(xTrans*interpFunc(xWave),dx=wavestep)
                    tried=True
                else:
                    y1=y2
            y2=(y3+y1)/2

            interpFunc=scint.interp1d(append(w[-1],x2),append(f[-1],y2))
            area=simps(xTrans*interpFunc(xWave),dx=wavestep)
            i+=1
        y1=f[-1]
    else:
        y1=f[-1]
        y2=0
    if abs(area-bArea)>.001*bArea:
        log.write('WARNING: The parameters you chose led to an integration in the %s-Band of %e instead of %e for index %i \n'%(band,vArea/area,color,index))
    (a,b,rval,pval,stderr)=stats.linregress(append(w[-1],x2),append(y1,y2))
    Nredstep = len( arange( w[-1], xWave[-1],  wavestep ) )
    wextRed =  sorted( [ w[-1] + (j+1)*wavestep for j in range(Nredstep) ] )
    fextRed = array( [ max( 0, a * wave + b ) if wave <= xWave[-1] else max(0,a*xWave[-1]+b) for wave in wextRed ] )
    w = append(w, wextRed)
    f = append(f, fextRed)
    return(w,f)
        



def extrapolatesed_linear(sedfile, newsedfile, color, niter=50):
    """ Interpolate the given transmission function to the wavestep of the SED, then
    get the area in the V band for color calculation and run the extrapolation algorithm.
    :param sedfile: SED file being extrapolated
    :param newsedfile: filename of output SED file
    :param fVH 
    """
    colors=[fVJ,fVH,fVK]
    dlist,wlist,flist = _getsed( sedfile ) 
    dlistnew, wlistnew, flistnew = [],[],[]

    vWave,vTrans=_getFilter('V')
    vInterpFunc=scint.interp1d(vWave,vTrans)
    fout = open( newsedfile, 'w' )
    log= open('./error.log','w')
    for i in range( len(dlist) ) : 
        d,w,f = dlist[i],wlist[i],flist[i]

        wavestep = w[1] - w[0]

        idx,val=_find_nearest(w,int(math.ceil(vWave[0]/10))*10)
        idx2,val2=_find_nearest(w,int(math.floor(vWave[-1]/10))*10)
        interp=vInterpFunc(arange(val,val2+1,wavestep))
        vArea=simps(f[idx:idx2+1]*interp,w[idx:idx2+1])#area in the v-band for color calculation
        fnew=deepcopy(f)
        wnew=deepcopy(w)
        for col in colors:
            if col:
                tempName=[ k for k,v in locals().iteritems() if v is col][0]
                tempWave,tempTrans=_getFilter(_colToFilter[tempName])
                interpFunc=scint.interp1d(tempWave,tempTrans)
                val=int(math.ceil(tempWave[0]/10))*10
                val2=int(math.floor(tempWave[-1]/10))*10
                tempWave=arange(val,val2+1,wavestep)
                tempTrans=interpFunc(tempWave)
                wnew,fnew=extrapolate_band(_colToFilter[tempName],vArea,col,tempTrans,tempWave,fnew,wnew,niter,log,i)
        for i in range( len( wnew ) ) :
            fout.write("%5.1f  %10i  %12.7e \n"%( d[0], wnew[i], fnew[i] ))
    fout.close() 
    log.close()
    return( newsedfile )

def extendNon1a(colorNames,colors,sedlist=None,showplots=None):
    """
    Function called in main which runs the extrapolation algorithm.
    :param sedlist: list of files to analyze (or none if all in current directory)
    :param fVJ: input V-J
    :param fVH: input V-H
    :param fVK: input V-K
    :param showplots: None if you don't want plotting, otherwise it's 'all' if you want to plot all days or a single day (i.e. 3)
    """

    if not sedlist:
        sedlist = glob.glob(os.path.join(sndataroot,"snsed","NON1A","*.SED"))
    else:
        sedlist=[os.path.join(sndataroot,"snsed","NON1A",sed) for sed in sedlist]
    for sedfile in sedlist : 
        newsedfile=os.path.basename(sedfile)
        print("EXTRAPOLATING %s"%newsedfile)
        if isinstance(colorNames,(tuple,list)) or isinstance(colors,(tuple,list)):
            if not isinstance(colors,(tuple,list)) or not isinstance(colorNames,(tuple,list)):
                raise RuntimeError('ColorNames and colors must both (or neither) be a list.')
            for color in colors:        
                extrapolatesed_linear(sedfile, newsedfile, color, niter=50 )
        if showplots:
            if showplots == 'all':
                plotsed(newsedfile,day=showplots)
            else:
                plotsed(newsedfile,day=showplots-9999)

        print("     Done with %s.\a\a\a"%os.path.basename(sedfile))

"""
I'm here. About to change it so that instead of defaulting to reading a transmission file, user will define a band name and sncosmo will find the band, or you can register
the band with sncosmo (I guess write a separate function to do this just to be comprehensive. Then I'm generalizing the code so that it doesn't just work for h,j,k but any
color and band combo you would want to choose.)
"""



def main():
    warnings.filterwarnings("ignore")
    for band in _filters:
        _filters[band]=_findFile(_filters[band])
    opts,args=getopt.getopt(sys.argv[1:],"i:p:v:j:h:k:",["vh=","vk=","jh=","jk=","vj=","vega"])
    mVJ=None
    mVH=None
    mVK=None
    mJH=None
    mJK=None
    showplots=None
    sedlist=None
    zpsys='AB'
    for opt,arg in opts:
        if opt == '-i':
            sedlist=arg.split(',')
        elif opt == '-p':
            showplots=arg
            if showplots != 'all':
                showplots=float(showplots)+9999 #this is just because I check if showplots exists later and if you choose zero it'll think that it doesn't exist
        elif opt == '-v':
            _filters[opt[-1].upper()]=_findFile(arg)
        elif opt == '-j':
            _filters[opt[-1].upper()]=_findFile(arg)
        elif opt == '-h':
            _filters[opt[-1].upper()]=_findFile(arg)
        elif opt == '-k':
            _filters[opt[-1].upper()]=_findFile(arg)
        elif opt == '--vh': #then v-h was given
            mVH=float(arg)
        elif opt=='--vk': #then v-k was given
            mVK=float(arg)
        elif opt == '--jh': #then j-h was given
            mJH=float(arg)
        elif opt == '--jk': #then j-k was given
            mJK=float(arg)
        elif opt == '--vj': #then v-j was given
            mVJ=float(arg)
        elif opt =='--vega':
            zpsys='Vega'
    if not mVH:
        if mVK and mJK and mJH:
            mVH=mVK-mJK+mJH
        elif mVJ and mJH:
            mVH=mVJ+mJH
    if not mVK: 
        if mVH and mJK and mJH:
            mVK=mVH-mJH+mJK
        elif mVJ and mJK:
            mVK=mVJ+mJK
    if not mVJ:
        if mVH and mJH:
            mVJ=mVH-mJH
        elif mVK and mJK:
            mVJ=mVK-mJK
    #calculate zero points from transmission files:
    _zp['AB']['V']=_getZP('V')
    #translate color to flux ratio
    if mVJ:
        _zp['AB']['J']=_getZP('J')
        fVJ=10**(-.4*(mVJ-(_zp[zpsys]['V']-_zp[zpsys]['J'])))
    else:
        fVJ=None
    if mVH:
        _zp['AB']['H']=_getZP('H')
        fVH=10**(-.4*(mVH-(_zp[zpsys]['V']-_zp[zpsys]['H'])))
    else:
        fVH=None
    if mVK:
        _zp['AB']['K']=_getZP('K')
        fVK=10**(-.4*(mVK-(_zp[zpsys]['V']-_zp[zpsys]['K'])))
    else: 
        fVK=None
    extendNon1a(sedlist,fVJ,fVH,fVK,showplots)


