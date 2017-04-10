#! /usr/bin/env python
#S.rodney & J.R. Pierel
# 2017.03.31

import os,sys,getopt,warnings
from numpy import *
from pylab import * 

#Automatically sets environment to your SNDATA_ROOT file (Assumes you've set this environment variable)
sndataroot = os.environ['SNDATA_ROOT']

#User can define center of V,H,K bands (angstrom)
VBAND=5500
HBAND=15414
KBAND=22000

#User can define width of V,H,K bands (angstrom)
vWidth=2000
hWidth=3390
kWidth=4000

#finds the left edge based on the two above variables
vLeftEdge=VBAND-vWidth/2
hLeftEdge=HBAND-hWidth/2
kLeftEdge=KBAND-kWidth/2


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

def extrapolatesed_linear(sedfile, newsedfile, iVH,iVK, Npt=4):
    """ use a linear fit of the first/last Npt points on the SED
    to extrapolate to H band (if necessary), then calculates the slope
    necessary across H and K bands to end up with the user-defined values
    for V-H and V-K
    """

    from scipy import interpolate as scint
    from scipy import stats
    from scipy.integrate import simps
    import shutil
    
    dlist,wlist,flist = getsed( sedfile ) 
    dlistnew, wlistnew, flistnew = [],[],[]

    fout = open( newsedfile, 'w' )
    log= open('./error.log','w')
    for i in range( len(dlist) ) : 
        d,w,f = dlist[i],wlist[i],flist[i]

        wavestep = w[1] - w[0]

        idx,val=find_nearest(w,vLeftEdge)
        nSteps = len( arange( val, vLeftEdge+vWidth,  wavestep ) )
        vArea=simps(f[idx:idx+nSteps+1],w[idx:idx+nSteps+1])
        print(vArea)
        #if the data ends short of H band, extrapolates to H band
        idx,val=find_nearest(w,hLeftEdge)
        if abs(val-hLeftEdge)>wavestep:
            # redward linear extrapolation from last N points
            wN = w[idx-Npt+1:idx+1]
            fN = f[idx-Npt+1:idx+1]
            (a,b,rval,pval,stderr)=stats.linregress(wN,fN)
            Nredstep = len( arange( val, hLeftEdge,  wavestep ) )
            wextRed =  sorted( [ val + (j+1)*wavestep for j in range(Nredstep) ] )
            fextRed = array( [ max( 0, a * wave + b ) for wave in wextRed ] )

            wnew = append(w[:idx+1], wextRed)
            fnew = append(f[:idx+1], fextRed)
        else:
            wnew=w[:idx+1]
            fnew=f[:idx+1]

        wN=wnew[-1]
        fN=fnew[-1]

        #geometrically matches area and then calculates necessary line slope
        hArea=vArea-iVH
        if hArea<0:
            log.write('WARNING: Input of V-H led to a negative area in H Band for day %i. \n'%i)
        if fN > 2*hArea/hWidth:
            log.write('WARNING: Only part of H band has positive flux for day %i \n'%(i+1))
            x=2*hArea/fN
            y=0
        else:
            x=hWidth
            y=2*hArea/hWidth-fN
        area=x*y+.5*x*(fN-y)
        if abs(area-hArea)>.01*hArea:
            log.write('WARNING: The parameters you chose led to an integration in the H-Band of %e instead of %e for day %i \n'%(vArea-area,iVH,i))
        (a,b,rval,pval,stderr)=stats.linregress(append(wN,wN+x),append(fN,y))

        Nredstep = len( arange( wN, kLeftEdge,  wavestep ) )
        wextRed =  sorted( [ wN + (j+1)*wavestep for j in range(Nredstep) ] )
        fextRed = array( [ max( 0, a * wave + b ) for wave in wextRed ] )
        
        wnew = append(wnew, wextRed)
        fnew = append(fnew, fextRed)

        wN=wnew[-1]
        fN=fnew[-1]
        
        #geometrically matches area and then calculates necessary line slope
        kArea=vArea-iVK
        if hArea<0:
            log.write('WARNING: Input of V-K led to a negative area in K Band for day %i. \n'%i)
        if fN > 2*kArea/kWidth:
            log.write('WARNING: Only part of K band has positive flux for day %i \n'%(i))
            x=2*kArea/fN
            y=0
        else:
            x=kWidth
            y=2*kArea/kWidth-fN
        area=x*y+.5*x*(fN-y)
        if abs(area-kArea)>.01*kArea:
            log.write('WARNING: The parameters you chose led to an integration in the K-Band of %e instead of %e for day %i \n'%(vArea-area,iVK,i))

        (a,b,rval,pval,stderr)=stats.linregress(append(wN,wN+x),append(fN,y))

        Nredstep = len( arange( wN, kLeftEdge+kWidth,  wavestep ) )
        wextRed =  sorted( [ wN + (j+1)*wavestep for j in range(Nredstep) ] )
        fextRed = array( [ max( 0, a * wave + b ) for wave in wextRed ] )
        
        wnew = append(wnew, wextRed)
        fnew = append(fnew, fextRed)

        for i in range( len( wnew ) ) :
            fout.write("%5.1f  %10i  %12.7e \n"%( d[0], wnew[i], fnew[i] ))
    fout.close() 
    log.close()
    return( newsedfile )

def extendNon1a(iVH,iVK,sedlist,showplots):
    """
    Function called in main which runs the extrapolation algorithm.
    :param iVH: input V-H
    :param iVK: input V-K
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
        extrapolatesed_linear(sedfile, newsedfile, iVH,iVK, Npt=4 )
        if showplots:
            plotsed(newsedfile,day=showplots-1)

        print("     Done with %s.\a\a\a"%os.path.basename(sedfile))


def main():
    warnings.filterwarnings("ignore")
    opts,args=getopt.getopt(sys.argv[1:],"i:p:",["vh=","vk=","jh=","jk="])
    iVH=None
    iVK=None
    showplots=None
    sedlist=None
    for opt,arg in opts:
        if opt == '-i':
            sedlist=arg.split(',')
        elif opt == '-p':
            showplots=arg
            if showplots != 'all':
                showplots=float(showplots)+1 #this is just because I check if showplots exists later and if you choose zero it'll think that it doesn't exist
        elif opt == '--vh': #then v-h was given
            iVH=float(arg)
        elif opt=='--vk': #then v-k was given
            iVK=float(arg)
        elif opt == '--jh': #then j-h was given
            JH=float(arg)
        elif opt == '--jk': #then j-k was given
            JK=float(arg)
    if not iVH and iVK and JK and JH:
        iVH=iVK-JK+JH
    elif not iVH:
        raise RuntimeError("V-H not given.")
    if not iVK and iVH and JK and JH:
        iVK=iVH-JH+JK
    elif not iVK:
        raise RuntimeError("V-K not given.")
    iVH=10**(-.4*(iVH))
    iVK=10**(-.4*(iVK))
    extendNon1a(iVH,iVK,sedlist,showplots)


if __name__ == '__main__':
    main()


