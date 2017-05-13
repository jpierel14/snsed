#! /usr/bin/env python
#S.rodney & J.R. Pierel
# 2017.03.31

import os,sys,getopt,warnings,math,glob
from numpy import *
from pylab import * 
from scipy import interpolate as scint
from scipy import stats
from scipy.integrate import simps
import shutil
from copy import deepcopy
#Automatically sets environment to your SNDATA_ROOT file (Assumes you've set this environment variable)
sndataroot = os.environ['SNDATA_ROOT']

#User can define center of V,H,K bands (angstrom)
filters={
    'V':'vBand/bessellv.dat',
    'J':'jBand/tophatJ.dat',
    'H':'hBand/tophatH.dat',
    'K':'kBand/tophatK.dat'
}

#dictionary of zero-points
zp={
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

colToFilter={
    'fVJ':'J',
    'fVH':'H',
    'fVK':'K'
}


def getsed(sedfile) : 
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


def plotsed( sedfile,day, normalize=False,MINWAVE=None,MAXWAVE=None,**kwarg): 
    """
    Simple plotting function to visualize the SED file.
    :param sedfile: SED filename to read (string)
    :param day: The day you want to plot, or 'all'
    :param normalize: Boolean to normalize the data
    :param MINWAVE: Allows user to plot in a specific window, left edge
    :param MAXWAVE: Allows user to plot in a specific window, right edge
    """

    dlist,wlist,flist = getsed( sedfile ) 

    for i in range( len(wlist) ) : 
        if MINWAVE:
            idx1,val= find_nearest(wlist[i],MINWAVE)
        else:
            idx1=None
        if MAXWAVE:
            idx2,val=find_nearest(wlist[i],MAXWAVE)
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


def find_nearest(array,value):
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

def getFilter(band):
    w,t=loadtxt(filters[band],unpack=True,comments='#')
    return(w,t)

def extrapolate_band(band,vArea,color,xTrans,xWave,f,w,niter,log,index):
    wavestep=w[1]-w[0]
    idx,val=find_nearest(w,xWave[0])
    if xWave[0]-val>wavestep:
        Nredstep = len( arange( val, xWave[0],  wavestep ) )
        wextRed =  sorted( [ val + (j+1)*wavestep for j in range(Nredstep) ] )
        fextRed = array( [ max(0,f[idx]) for wave in wextRed ] )
        w = append(w[:idx+1], wextRed)
        f = append(f[:idx+1], fextRed)

    else:
        w=w[:idx+1]
        f=f[:idx+1]
    wN=w[-1]
    fN=f[-1]
    bArea=vArea/color
    x2=xWave[-1]
    x1=wN
    interpFunc=scint.interp1d(append(wN,x2),append(fN,0))
    area=simps(xTrans*interpFunc(xWave),dx=wavestep)
    temp=.5*(xWave[-1]-xWave[0])*fN
    #(m,b,rval,pval,stderr)=stats.linregress(append(wN,x2),append(fN,0))
    i=1
    if area>bArea:
        y2=fN
        y1=fN
        last=0
        while(abs(area-bArea)/(area+bArea)>.00000001 and i<niter):
            if area>bArea:
                x3=x2
                x2=(x3+x1)/2
                idx,x2=find_nearest(xWave,x2)
            else:
                x1=x2
                x2=(x3+x1)/2
                idx,x2=find_nearest(xWave,x2)
            interpFunc=scint.interp1d(append(wN,x2),append(fN,0))
            area=simps(xTrans[0:idx+1]*interpFunc(arange(wN,x2+1,wavestep)),dx=wavestep)
            i+=1
            if area==last:
                break
            last=area

        i=1
        while(abs(area-bArea)/(area+bArea)>.00000001 and i<niter):
            if i==1:
                y1=.5*fN
                y2=fN
                y3=1.5*fN
            if area>bArea:
                y3=y2
                y2=(y3+y1)/2
            else:
                y1=y2
                y2=(y3+y1)/2
            interpFunc=scint.interp1d(append(wN,x2),append(y2,0))
            area=simps(xTrans[0:idx+1]*interpFunc(arange(wN,x2+1,wavestep)),dx=wavestep)
            i+=1
        y1=y2
        y2=0
    elif area<bArea:
        y1=0
        y3=max(f)
        y2=y3
        x2=xWave[-1]
        tried=False
        while(abs(area-bArea)/(area+bArea)>.00000001 and i<niter):
            if area>bArea:
                y3=y2
                y2=(y3+y1)/2
            else:
                if not tried:
                    while area < bArea:
                        y3*=2
                        interpFunc=scint.interp1d(append(wN,x2),append(fN,y3))
                        area=simps(xTrans*interpFunc(xWave),dx=wavestep)
                    y2=(y3+y1)/2
                    tried=True
                else:
                    y1=y2
                    y2=(y3+y1)/2

            interpFunc=scint.interp1d(append(wN,x2),append(fN,y2))
            area=simps(xTrans*interpFunc(xWave),dx=wavestep)
            i+=1
        y1=fN
    else:
        y1=fN
        y2=0
    if abs(area-bArea)>.001*bArea:
        log.write('WARNING: The parameters you chose led to an integration in the %s-Band of %e instead of %e for day %i \n'%(band,vArea/area,color,index))
    (a,b,rval,pval,stderr)=stats.linregress(append(wN,x2),append(y1,y2))
    Nredstep = len( arange( wN, xWave[-1],  wavestep ) )
    wextRed =  sorted( [ wN + (j+1)*wavestep for j in range(Nredstep) ] )
    fextRed = array( [ max( 0, a * wave + b ) if wave <= xWave[-1] else max(0,a*xWave[-1]+b) for wave in wextRed ] )
    w = append(w, wextRed)
    f = append(f, fextRed)
    return(w,f)
        



def extrapolatesed_linear(sedfile, newsedfile, fVH,fVK,fVJ, niter=15):
    """ use a linear fit of the first/last Npt points on the SED
    to extrapolate to H band (if necessary), then calculates the slope
    necessary across H and K bands to end up with the user-defined values
    for V-H and V-K
    """
    colors=[fVJ,fVH,fVK]
    dlist,wlist,flist = getsed( sedfile ) 
    dlistnew, wlistnew, flistnew = [],[],[]

    vWave,vTrans=getFilter('V')
    vInterpFunc=scint.interp1d(vWave,vTrans)
    fout = open( newsedfile, 'w' )
    log= open('./error.log','w')
    for i in range( len(dlist) ) : 
        d,w,f = dlist[i],wlist[i],flist[i]

        wavestep = w[1] - w[0]

        idx,val=find_nearest(w,int(math.ceil(vWave[0]/10))*10)
        idx2,val2=find_nearest(w,int(math.floor(vWave[-1]/10))*10)
        interp=vInterpFunc(arange(val,val2+1,wavestep))
        vArea=simps(f[idx:idx2+1]*interp,w[idx:idx2+1])
        fnew=deepcopy(f)
        wnew=deepcopy(w)
        for col in colors:
            if col:
                tempName=[ k for k,v in locals().iteritems() if v is col][0]
                tempWave,tempTrans=getFilter(colToFilter[tempName])
                interpFunc=scint.interp1d(tempWave,tempTrans)
                val=int(math.ceil(tempWave[0]/10))*10
                val2=int(math.floor(tempWave[-1]/10))*10
                tempWave=arange(val,val2+1,wavestep)
                tempTrans=interpFunc(tempWave)
                wnew,fnew=extrapolate_band(colToFilter[tempName],vArea,col,tempTrans,tempWave,fnew,wnew,niter,log,i)
        for i in range( len( wnew ) ) :
            fout.write("%8.4f  %10i  %12.7e \n"%( d[0], wnew[i], fnew[i] ))
    fout.close() 
    log.close()
    return( newsedfile )

def extendNon1a(fVH,fVK,fVJ,sedlist,showplots):
    """
    Function called in main which runs the extrapolation algorithm.
    :param fVH: input V-H
    :param fVK: input V-K
    :param sedlist: list of files to analyze (or none if all in current directory)
    :param showplots: None if you don't want plotting, otherwise it's 'all' if you want to plot all days or a single day (i.e. 3)
    """
    import glob
    import shutil

    if not sedlist:
        sedlist = glob.glob(os.path.join(sndataroot,"snsed","NON1A","*.SED"))
    else:
        sedlist=[os.path.join(sndataroot,"snsed","NON1A",sed) for sed in sedlist]
    for sedfile in sedlist : 
        newsedfile=os.path.basename(sedfile)
        print("EXTRAPOLATING %s"%os.path.basename(sedfile))
        extrapolatesed_linear(sedfile, newsedfile, fVH,fVK,fVJ, niter=50 )
        if showplots:
            if showplots == 'all':
                plotsed(newsedfile,day=showplots)
            else:
                plotsed(newsedfile,day=showplots-9999)

        print("     Done with %s.\a\a\a"%os.path.basename(sedfile))


def getZP(band):
    c=299792458
    wave,trans=getFilter(band)
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

def findFile(filename):
    for root, subFolders, files in os.walk('.'):
        if filename in files:
            return(os.path.abspath(os.path.join(root,filename)))
def main():
    warnings.filterwarnings("ignore")
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
            filters[opt[-1].upper()]=findFile(arg)
        elif opt == '-j':
            filters[opt[-1].upper()]=findFile(arg)
        elif opt == '-h':
            filters[opt[-1].upper()]=findFile(arg)
        elif opt == '-k':
            filters[opt[-1].upper()]=findFile(arg)
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
    zp['AB']['V']=getZP('V')
    #translate color to flux ratio
    if mVJ:
        zp['AB']['J']=getZP('J')
        fVJ=10**(-.4*(mVJ-(zp[zpsys]['V']-zp[zpsys]['J'])))
    else:
        fVJ=None
    if mVH:
        zp['AB']['H']=getZP('H')
        fVH=10**(-.4*(mVH-(zp[zpsys]['V']-zp[zpsys]['H'])))
    else:
        fVH=None
    if mVK:
        zp['AB']['K']=getZP('K')
        fVK=10**(-.4*(mVK-(zp[zpsys]['V']-zp[zpsys]['K'])))
    else: 
        fVK=None
    extendNon1a(fVH,fVK,fVJ,sedlist,showplots)


if __name__ == '__main__':
    main()


