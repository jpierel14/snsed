#! /usr/bin/env python
#S.rodney & J.R. Pierel
# 2017.03.31

import os,sys,getopt,warnings
from numpy import *
from pylab import * 
from scipy import interpolate as scint
from scipy import stats
from scipy.integrate import simps
import shutil
#Automatically sets environment to your SNDATA_ROOT file (Assumes you've set this environment variable)
sndataroot = os.environ['SNDATA_ROOT']

#User can define center of V,H,K bands (angstrom)
JBAND=12355
HBAND=15414
KBAND=22000

#User can define width of V,H,K bands (angstrom)
jWidth=3000
hWidth=3390
kWidth=4000

#finds the left edge based on the two above variables
jLeftEdge=JBAND-jWidth/2
hLeftEdge=HBAND-hWidth/2
kLeftEdge=KBAND-kWidth/2

#dictionary of zero-points
zp={
    'AB':{
        'V':-13.75,
        'J':-14.16,
        'H':-14.51,
        'K':-15.11
        },
    'Vega':{
        'V':-13.77,
        'J':-15.07,
        'H':-15.9,
        'K':16.96   
        }
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


def extrapolate_band(band,vArea,xArea,fN,wN,width,rightEdge,nextEdge,wavestep,log,index):
    #geometrically matches area and then calculates necessary line slope
        bArea=vArea/xArea
        if fN*width > bArea: #negative slope
            if .5*fN*width>bArea:
                log.write('WARNING: Only part of %s band has positive flux for day %i \n'%(band,index+1))
                x=2*bArea/fN
                y=0

            else:
                x=width
                y=2*bArea/width-fN
            area=.5*(fN-y)*x+y*x
        else:
            x=width
            y=2*bArea/width-fN
            area=x*fN+.5*x*(y-fN)
        if abs(area-bArea)>.01*bArea:
            log.write('WARNING: The parameters you chose led to an integration in the %s-Band of %e instead of %e for day %i \n'%(band,vArea-area,xArea,index))
        (a,b,rval,pval,stderr)=stats.linregress(append(wN,wN+x),append(fN,y))

        Nredstep = len( arange( wN, nextEdge,  wavestep ) )
        wextRed =  sorted( [ wN + (j+1)*wavestep for j in range(Nredstep) ] )
        fextRed = array( [ max( 0, a * wave + b ) if wave <= rightEdge else max(0,a*rightEdge+b) for wave in wextRed ] )
        return(wextRed,fextRed)
        



def extrapolatesed_linear(sedfile, newsedfile, fVH,fVK,fVJ, Npt=4):
    """ use a linear fit of the first/last Npt points on the SED
    to extrapolate to H band (if necessary), then calculates the slope
    necessary across H and K bands to end up with the user-defined values
    for V-H and V-K
    """
    
    dlist,wlist,flist = getsed( sedfile ) 
    dlistnew, wlistnew, flistnew = [],[],[]

    vWave=np.loadtxt('vWave.dat')
    vTrans=np.loadtxt('vTrans.dat')
    vLeftEdge=vWave[0]
    vWidth=vWave[-1]-vWave[0]
    interpFunc=scint.interp1d(vWave,vTrans)
    fout = open( newsedfile, 'w' )
    log= open('./error.log','w')
    for i in range( len(dlist) ) : 
        d,w,f = dlist[i],wlist[i],flist[i]

        wavestep = w[1] - w[0]

        idx,val=find_nearest(w,vLeftEdge)
        idx2,val2=find_nearest(w,vLeftEdge+vWidth)
        nSteps = len( arange( val, vLeftEdge+vWidth,  wavestep ) )
        interp=interpFunc(np.arange(val,val2+1,10))
        
        vArea=simps(f[idx:idx+nSteps+1]*interp,w[idx:idx+nSteps+1])
        #if the data ends short of next band, extrapolates to next band
        if fVJ:
            idx,val=find_nearest(w,jLeftEdge)
            edge=jLeftEdge
        elif fVH:
            idx,val=find_nearest(w,hLeftEdge)
            edge=hLeftEdge
        elif fVK:
            idx,val=find_nearest(w,kLeftEdge)
            edge=kLeftEdge
        if abs(val-edge)>wavestep:
            '''
            # redward linear extrapolation from last N points
            wN = w[idx-Npt+1:idx+1]
            fN = f[idx-Npt+1:idx+1]
            (a,b,rval,pval,stderr)=stats.linregress(wN,fN)
            '''
            Nredstep = len( arange( val, edge,  wavestep ) )
            wextRed =  sorted( [ val + (j+1)*wavestep for j in range(Nredstep) ] )
            fextRed = array( [ max( 0, f[idx] ) for wave in wextRed ] )

            wnew = append(w[:idx+1], wextRed)
            fnew = append(f[:idx+1], fextRed)
        else:
            wnew=w[:idx+1]
            fnew=f[:idx+1]

        wN=wnew[-1]
        fN=fnew[-1]
        if fVJ:
            if fVH: 
                edge=hLeftEdge
            elif fVK:
                edge=kLeftEdge
            else:
                edge=max(np.max(w),jLeftEdge+jWidth)
            wextRed,fextRed=extrapolate_band('J',vArea,fVJ,fN,wN,jWidth,jLeftEdge+jWidth,edge,wavestep,log,i)
            wnew = append(wnew, wextRed)
            fnew = append(fnew, fextRed)

            wN=wnew[-1]
            fN=fnew[-1]
        if fVH:
            if fVK:
                edge=kLeftEdge
            else:
                edge=max(np.max(w),hLeftEdge+hWidth)
            wextRed,fextRed=extrapolate_band('H',vArea,fVH,fN,wN,hWidth,hLeftEdge+hWidth,edge,wavestep,log,i)
            wnew = append(wnew, wextRed)
            fnew = append(fnew, fextRed)

            wN=wnew[-1]
            fN=fnew[-1]
        if fVK:
            wextRed,fextRed=extrapolate_band('K',vArea,fVK,fN,wN,kWidth,kLeftEdge+kWidth,max(np.max(w),kLeftEdge+kWidth),wavestep,log,i)
            wnew = append(wnew, wextRed)
            fnew = append(fnew, fextRed)
        
        for i in range( len( wnew ) ) :
            fout.write("%5.1f  %10i  %12.7e \n"%( d[0], wnew[i], fnew[i] ))
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
        extrapolatesed_linear(sedfile, newsedfile, fVH,fVK,fVJ, Npt=4 )
        if showplots:
            if showplots == 'all':
                plotsed(newsedfile,day=showplots)
            else:
                plotsed(newsedfile,day=showplots-9999)

        print("     Done with %s.\a\a\a"%os.path.basename(sedfile))


def main():
    warnings.filterwarnings("ignore")
    opts,args=getopt.getopt(sys.argv[1:],"i:p:",["vh=","vk=","jh=","jk=","vj=","vega"])
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
    #translate color to flux ratio
    if mVJ:
        fVJ=10**(-.4*(mVJ-(zp[zpsys]['V']-zp[zpsys]['J'])))
    else:
        fVJ=None
    if mVH:
        fVH=10**(-.4*(mVH-(zp[zpsys]['V']-zp[zpsys]['H'])))
    else:
        fVH=None
    if mVK:
        fVK=10**(-.4*(mVK-(zp[zpsys]['V']-zp[zpsys]['K'])))
    else: 
        fVK=None
    extendNon1a(fVH,fVK,fVJ,sedlist,showplots)


if __name__ == '__main__':
    main()


