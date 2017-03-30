#! /usr/bin/env python
#S.rodney
# 2011.05.04
"""
Extrapolate the Hsiao SED down to 300 angstroms
to allow the W filter to reach out to z=2.5 smoothly
in the k-correction tables
"""

import os,sys,getopt
from numpy import *
from pylab import * 

#sndataroot = os.environ['SNDATA_ROOT']
# sndataroot = os.environ['NON1A']

MINWAVE = 300     # min wavelength for extrapolation (Angstroms)
MAXWAVE = 18000   # max wavelength for extrapolation (Angstroms)

VBAND=5500
HBAND=15414
KBAND=22000

vWidth=3000
hWidth=17109-13719
kWidth=4000

vLeftEdge=VBAND-vWidth/2
hLeftEdge=HBAND-hWidth/2
kLeftEdge=KBAND-kWidth/2


def getsed(sedfile) : 
    d,w,f = loadtxt( sedfile, unpack=True,comments='#') 
    
    #d = d.astype(int)
    days = unique( d ) 

    dlist = [ d[ where( d == day ) ] for day in days ]
    wlist = [ w[ where( d == day ) ] for day in days ]
    flist = [ f[ where( d == day ) ] for day in days ]

    return( dlist, wlist, flist )

def plotsed( sedfile,day, normalize=False,MINWAVE=None,MAXWAVE=None,**kwarg): 
    dlist,wlist,flist = getsed( sedfile ) 
    #days = unique( dlist ) 
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
        #defaults = { 'label':str(thisday) } 
        #plotarg = dict( kwarg.items() + defaults.items() )
        if day!='all' : 
            if abs(thisday-day)>0.6 : continue
        if normalize : 
            plot( wlist[i], flist[i]/flist[i].max()+thisday, **kwarg )
            show()
        else : 
            #idx,val=find_nearest(wlist[i],hLeftEdge)
            #idx2,val=find_nearest(wlist[i],kLeftEdge+kWidth)
            plot( wlist[i][idx1:idx2], flist[i][idx1:idx2], label=str(thisday), **kwarg )
            show()
        # user_in=raw_input('%i : return to continue'%i)

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1,array[idx-1]
    else:
        return idx,array[idx]

def extrapolatesed_linear(sedfile, newsedfile, iVH,iVK, maxwave=MAXWAVE, Npt=4):
    """ use a linear fit of the first/last Npt  points on the SED
    to extrapolate  """

    from scipy import interpolate as scint
    from scipy import stats
    from scipy.integrate import simps
    import shutil
    
    dlist,wlist,flist = getsed( sedfile ) 
    dlistnew, wlistnew, flistnew = [],[],[]

    fout = open( newsedfile, 'w' )
    log= open('./error.log','w')
    error=False
    for i in range( len(dlist) ) : 
        d,w,f = dlist[i],wlist[i],flist[i]

        wavestep = w[1] - w[0]

        idx,val=find_nearest(w,vLeftEdge)
        nSteps = len( arange( val, vLeftEdge+vWidth,  wavestep ) )
        vArea=simps(f[idx:idx+nSteps+1],w[idx:idx+nSteps+1])
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
        
        if wN !=val or fN == 0:
            raise RuntimeError("Extrapolation from end of data to H band failed.")

        hArea=vArea-iVH
        if fN > 2*hArea/hWidth:
            log.write('WARNING: Only part of H band has positive flux for day %i \n'%(i+1))
            error=True
            x=2*hArea/fN
            y=0
        else:
            x=hWidth
            y=2*hArea/hWidth-fN
        area=x*y+.5*x*(fN-y)
        
        if abs(area-hArea)>.01*hArea:
            log.write('WARNING: The parameters you chose led to an integration in the H-Band of %e instead of %e for day %i \n'%(vArea-area,iVH,i+1))
            error=True
        (a,b,rval,pval,stderr)=stats.linregress(append(wN,wN+x),append(fN,y))

        Nredstep = len( arange( wN, kLeftEdge,  wavestep ) )
        wextRed =  sorted( [ wN + (j+1)*wavestep for j in range(Nredstep) ] )
        fextRed = array( [ max( 0, a * wave + b ) for wave in wextRed ] )
        
        wnew = append(wnew, wextRed)
        fnew = append(fnew, fextRed)

        wN=wnew[-1]
        fN=fnew[-1]

        if wN !=kLeftEdge or fN == 0:
            raise RuntimeError("Extrapolation from end of data to K band failed.")
    
        kArea=vArea-iVK
        if fN > 2*kArea/kWidth:
            log.write('WARNING: Only part of K band has positive flux for day %i \n'%(i+1))
            error=True
            x=2*kArea/fN
            y=0
        else:
            x=kWidth
            y=2*kArea/kWidth-fN
        area=x*y+.5*x*(fN-y)
        
        if abs(area-kArea)>.01*kArea:
            log.write('WARNING: The parameters you chose led to an integration in the K-Band of %e instead of %e for day %i \n'%(vArea-area,iVK,i+1))
            error=True

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
    import glob
    import shutil
    sedlist = glob.glob("./*.SED") if not sedlist else sedlist
    #sedlist=['sn91T-nugent.SED']
    for sedfile in sedlist : 
        #newsedfile =  'non1a/' + os.path.basename( sedfile )
        #newsedfile =  os.path.basename( sedfile )
        newsedfile=os.path.splitext(os.path.basename( sedfile ))[0]+'EXT'+os.path.splitext(os.path.basename( sedfile ))[1]

        print("EXTRAPOLATING %s"%sedfile)
        extrapolatesed_linear(sedfile, newsedfile, iVH,iVK, maxwave=MAXWAVE, Npt=4 )
        if showplots:
            plotsed(newsedfile,day=showplots-1)

        print("     Done with %s.\a\a\a"%sedfile)


def main():
    opts,args=getopt.getopt(sys.argv[1:],"i:p:",["vh=","vk=","hj=","jk="])
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
        elif opt == '--hj': #then h-j was given
            HJ=float(arg)
        elif opt == '--jk': #then j-k was given
            JK=float(arg)
    if not iVH and iVK and JK and HJ:
        iVH=iVK-JK-HJ
    elif not iVH:
        raise RuntimeError("V-H not given.")
    if not iVK and iVH and JK and HJ:
        iVK=iVH+HJ+JK
    elif not iVK:
        raise RuntimeError("V-K not given.")
    
    extendNon1a(iVH,iVK,sedlist,showplots)


if __name__ == '__main__':
    main()


