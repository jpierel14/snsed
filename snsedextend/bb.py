from __future__ import print_function
import sncosmo,glob,os,sys,math,snsedextend
import numpy as np
import astropy.units as u
from astropy.modeling.blackbody import FLAM
from astropy.modeling.models import BlackBody1D
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter as smooth
from astropy.io import ascii
from astropy.table import Table

from .utils import _filters
from .helpers import _find_nearest


def chisquared(obs,actual):
    return(np.sum((obs-actual)**2))

def bbChi(stuff,wave,flux):
    temp,const=stuff
    bb=BlackBody1D(temperature=temp*u.K)
    bbFlux=[x.value for x in bb(wave).to(FLAM,u.spectral_density(wave))]

    bbFlux=[x*const for x in bbFlux]
    return chisquared(bbFlux,flux)

def bbOut(table,name):
    ascii.write(table,name)

def fluxFromBB(temp,const,wave):
    bb=BlackBody1D(temperature=temp*u.K)
    bbFlux=[x.value for x in bb(wave*u.AA).to(FLAM,u.spectral_density(wave*u.AA))]
    bbFlux=np.array([x*const for x in bbFlux])
    return(bbFlux)

def getBB(phase,wave,flux,name):
    
    #fluxes=[]
    temperatures=[]
    constants=[]
    for p in range(len(phase)):
        flux2=flux[p][wave<9000][wave[wave<9000]>4200]
        wave2=wave[wave<9000][wave[wave<9000]>4200]*u.AA
        res=minimize(bbChi,np.array([6000,1]),args=(wave2,smooth(flux2,len(flux2)/3,2)),bounds=((0,None),(0,None)))

        temp,const=res.x
        bbFlux=fluxFromBB(temp,const,wave)

        temperatures.append(temp)
        constants.append(const)

    temperatures=np.array(temperatures)
    constants=np.array(constants)
    tempBB=Table([phase,temperatures,constants],names=['phase','temp','amplitude'])
    #bbOut(tempBB,name+'bb')
    return(tempBB)

def bbExtrap(bb,phase,wave,flux,name,typ,startWave=9000,endWave=55000):#
    extrapWave=np.arange(_find_nearest(wave,startWave)[1],endWave,wave[1]-wave[0])
    finalWave=np.append(wave[wave<extrapWave[0]],extrapWave)
    newFlux=[]
    for p in range(len(phase)):
        bbFlux=fluxFromBB(bb['temp'][bb['phase']==phase[p]],bb['amplitude'][bb['phase']==phase[p]],extrapWave)
        newFlux.append(np.append(flux[p][wave<extrapWave[0]],bbFlux))

    newFlux=np.array(newFlux)
    return(phase,finalWave,newFlux)

def standardize(data):
    bandDict=deepcopy(_filters)
    for band in bandDict:
        if not isinstance(bandDict[band],sncosmo.Bandpass):
            bandDict=snsed.colorCalc._bandCheck(bandDict,band)
    zpMag=sncosmo.get_magsystem('vega')
    temp=[bandDict[data['band'][l]].name for l in range(len(data))]
    data.remove_column('band')
    data['band']=temp
    #now make sure we have zero-points and fluxes for everything
    if 'zp' not in data.colnames:
        data['zp']=[zpMag.band_flux_to_mag(1/sncosmo.constants.HC_ERG_AA,data['band'][i]) for i in range(len(data))]
    if 'flux' not in data.colnames:
        data=snsed.mag_to_flux(data,bandDict,'vega')
    if 'zpsys' not in data.colnames:
        data['zpsys']=['vega' for i in range(len(data))]
    return(data)

def bbColorcurve(model,c1,c2,typ,zpsys='AB'):
    col=model.color(c1,_filters[c2],zpsys,np.arange(-10,40,1))
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(np.arange(-10,40,1),col)
    plt.savefig(os.path.join('figs','bbCurves',typ,c2,model._source.name+'.pdf'),format='pdf',overwrite=True)
    #sys.exit()

def main():
    snType=['Ib','Ic','II']
    bands=['J','H','K']
    
    for typ in snType:
        seds=glob.glob('../SED_Repository/Type_'+typ+'/*.SED')
        for sed in seds:
            phase,wave,flux=sncosmo.read_griddata_ascii(sed)
            bb=getBB(phase,wave,flux,os.path.basename(sed)[:-4])
            ePhase,eWave,eFlux=bbExtrap(bb,phase,wave,flux,os.path.basename(sed)[:-4],typ)
            tempSource=sncosmo.TimeSeriesSource(ePhase,eWave,eFlux,name=os.path.basename(sed)[:-4])
            for band in bands:
                bbColorcurve(sncosmo.Model(tempSource),'sdss::r',band,typ,zpsys='vega')
    sys.exit()
    
    bands=['J','H','K']
    curveDict={'Ib':['lc_2005hg.dat','lc_2006fo.dat','lc_2007C.dat'],'Ic':['lc_2006aj.dat'],'II':['lc_2010bq.dat']}
    #typeDict={'Ib':['SDSS-019323.SED'],'Ic':'CSP-2004fe.SED','II':'SDSS-018457.SED'}
    mod,sed=np.loadtxt('sedSN.ref',dtype='str',unpack=True)
    sedDict={mod[i]:sed[i]+'.SED' for i in range(len(mod))}
    
    dirDict={'Ib':'bianco','Ic':'bianco','II':'hicken'}
    for band in bands:
        for typ in snType:
            directory=os.path.join('..','snsedextend',dirDict[typ])
            for lc in curveDict[typ]:
                sn,redshift=np.loadtxt(os.path.join(directory,'redshift.ref'),dtype='str',unpack=True)
                redshift={sn[i]:redshift[i] for i in range(len(sn))}
                sn,colors=np.loadtxt(os.path.join(directory,'color.ref'),dtype='str',unpack=True)
                colors={sn[i]:colors[i] for i in range(len(sn))}
                sn,dusts=np.loadtxt(os.path.join(directory,'dust.ref'),dtype='str',unpack=True)
                dust={sn[i]:dusts[i] for i in range(len(sn))}
                sn,types=np.loadtxt(os.path.join(directory,'type.ref'),dtype='str',unpack=True)
                alltyp={sn[i]:types[i] for i in range(len(sn))}
                sne,peaks=np.loadtxt(os.path.join(directory,'timeOfPeak.ref'),dtype='str',unpack=True)
                peaks={sne[i]:float(peaks[i]) for i in range(len(sne))}
                fit,res,t0,optical,nir=snsedextend.curveToColor(os.path.join(directory,lc),['r-'+band],snType=alltyp[lc[:-4]],zpsys='Vega',
                                          bounds={'hostebv':(-1,1),'t0':(peaks[lc[:-4]]-5,peaks[lc[:-4]]+5)},
                                          constants={'z':redshift[lc[:-4]],'hostr_v':3.1,'mwr_v':3.1,'mwebv':dust[lc[:-4]]},
                                          dust='CCM89Dust',effect_frames=['rest','obs'],effect_names=['host','mw'])
                #sncosmo.plot_lc(model=fit,bands=['bessellb','bessellv','sdss::r'])
                #plt.show()
                #sys.exit()
            
                phase,wave,flux=sncosmo.read_griddata_ascii(os.path.join('..','SED_Repository','Type_'+typ,sedDict[fit._source.name]))
                bb=getBB(phase,wave,flux,os.path.join('bbDat',sedDict[fit._source.name]))
                ePhase,eWave,eFlux=bbExtrap(bb,phase,wave,flux)
                tempSource=sncosmo.TimeSeriesSource(ePhase,eWave,eFlux)
                #bbColorcurve(sncosmo.Model(tempSource),'sdss::r','paritel::j')
                effect_frames=['rest','obs']
                effect_names=['host','mw']
                dust_dict={'SFD98Map':sncosmo.SFD98Map,'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
                dust=dust_dict['CCM89Dust']()
                effects=[dust for i in range(len(effect_names))] 
                eMod=sncosmo.Model(source=tempSource,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
          
                
                eMod.set(**{x:y for x,y in zip(fit.param_names,fit.parameters)})
                #eMod.set(t0=eMod.get('t0')+5)
                #finalres,final=sncosmo.fit_lc(nir,eMod,['t0'])
                sncosmo.plot_lc(data=nir,model=eMod)
                plt.savefig(os.path.join('figs','lcs',typ,band,lc[:-4]+'_bb.pdf'),format='pdf',overwrite=True)

                for f in ['red','blue','median']:
                    phase,wave,flux=sncosmo.read_griddata_ascii(os.path.join('bbDat',sedDict[fit._source.name]+'_'+f))
                    tempSource=sncosmo.TimeSeriesSource(phase,wave,flux)
                    effect_frames=['rest','obs']
                    effect_names=['host','mw']
                    dust_dict={'SFD98Map':sncosmo.SFD98Map,'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
                    dust=dust_dict['CCM89Dust']()
                    effects=[dust for i in range(len(effect_names))] 
                    eMod=sncosmo.Model(source=tempSource,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
              
                    
                    eMod.set(**{x:y for x,y in zip(fit.param_names,fit.parameters)})
                    #eMod.set(t0=eMod.get('t0')+5)
                    #finalres,final=sncosmo.fit_lc(nir,eMod,['t0'])
                    sncosmo.plot_lc(data=nir,model=eMod)
                    plt.savefig(os.path.join('figs','lcs',typ,band,lc[:-4]+'_'+f+'.pdf'),format='pdf',overwrite=True)


if __name__=='__main__':
    main()