#! /usr/bin/env python
#S.rodney & J.R. Pierel
# 2017.03.31


import os,sys,sncosmo,snsedextend,warnings,importlib
from numpy import *
from scipy import interpolate as scint
from astropy.table import Table,Column,MaskedColumn,vstack
import pyParz

from .utils import *
from .helpers import *
from sncosmo.photdata import photometric_data

__all__=['mag_to_flux','flux_to_mag','curveToColor','colorTableCombine']






def _getErrorFromModel(args,zpsys,band):
    """
    (Private)
    Helper function to obtain a magnitude error estimate from an sncosmo model
    """
    model,curve,method,params,bounds,ignore,constants,dust,effect_names,effect_frames=args
    curve[_get_default_prop_name('flux')]=curve[_get_default_prop_name('flux')]+curve[_get_default_prop_name('fluxerr')]
    res,fit=_snFit((model,curve,method,params,bounds,ignore,constants,dust,effect_names,effect_frames))
    ugrid,umag=_snmodel_to_mag(fit,curve,zpsys,band)
    curve[_get_default_prop_name('flux')]=curve[_get_default_prop_name('flux')]-2*curve[_get_default_prop_name('fluxerr')]
    res,fit=_snFit((model,curve,method,params,bounds,ignore,constants,dust,effect_names,effect_frames))
    lgrid,lmag=_snmodel_to_mag(fit,curve,zpsys,band)
    return (ugrid,umag,lgrid,lmag)


def _getBandMaxTime(model,table,bands,band,zp,zpsys):
    """
    (Private)
     Helper function to obtain the time of maximum flux from an sncosmo model
    """

    tgrid,mflux=_snmodel_to_flux(model,table,zp,zpsys,bands[band])
    return(tgrid[where(mflux==max(mflux))])

def _snFit(args):
    """
    (Private)
    Function that does the fitting for the curveToColor method. It only uses the simple sncosmo fit_lc method, so it
    is assuming that the parameters you're feeding it are reasonably close to correct already and that the data is
    clipped to somethign reasonable. This function is used by the parallelization function as well.
    """
    warnings.simplefilter('ignore')
    sn_func={'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
    dust_dict={'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
    model,curve,method,params,bounds,ignore,constants,dust,effect_names,effect_frames=args
    if dust:
        dust=dust_dict[dust]()
    else:
        dust=[]
    bounds=bounds if bounds else {}
    ignore=ignore if ignore else []
    effects=[dust for i in range(len(effect_names))] if effect_names else []
    effect_names=effect_names if effect_names else []
    effect_frames=effect_frames if effect_frames else []
    if not isinstance(effect_names,(list,tuple)):
        effects=[effect_names]
    if not isinstance(effect_frames,(list,tuple)):
        effects=[effect_frames]
    model=sncosmo.Model(source=model,effects=effects,effect_names=effect_names,effect_frames=effect_frames)
    params=params if params else [x for x in model.param_names]
    no_bound = {x for x in params if x in _needs_bounds and x not in bounds.keys() and x not in constants.keys()}
    if no_bound:
        params=list(set(params)-no_bound)
    params= [x for x in params if x not in ignore and x not in constants.keys()]
    if constants:
        constants = {x:constants[x] for x in constants.keys() if x in model.param_names}
        model.set(**constants)

    if method=='nest' and 'amplitude' in params and 'amplitude' not in bounds.keys():
        res,fit=sn_func[method](curve, model, params,guess_amplitude_bound=True, bounds=bounds,verbose=False)
    else:
        res,fit=sn_func[method](curve, model, params, bounds=bounds,verbose=False)

    return(pyParz.parReturn((res,fit)))#parallel.parReturn makes sure that whatever is being returned is pickleable

def _snmodel_to_mag(model,table,zpsys,band):
    """
    (Private)
     Helper function that takes an sncosmo model in flux-space and a lightcurve astropy.Table object, returns a time
     grid and magnitude vector.
    """
    warnings.simplefilter("ignore")
    tmin = []
    tmax = []
    tmin.append(min(table[_get_default_prop_name('time')]) - 10)
    tmax.append(max(table[_get_default_prop_name('time')]) + 10)
    tmin.append(model.mintime())
    tmax.append(model.maxtime())
    tmin = min(tmin)
    tmax = max(tmax)
    tgrid = linspace(tmin, tmax, int(tmax - tmin) + 1)
    mMag = model.bandmag(band, zpsys,tgrid)
    return(tgrid,mMag)

def _snmodel_to_flux(model,table,zp,zpsys,band):
    """
    (Private)
     Helper function that takes an sncosmo model in flux-space and a lightcurve astropy.Table object, returns a time
     grid and flux vector.
    """
    warnings.simplefilter("ignore")
    tmin = []
    tmax = []
    tmin.append(min(table[_get_default_prop_name('time')]) - 10)
    tmax.append(max(table[_get_default_prop_name('time')]) + 10)
    tmin.append(model.mintime())
    tmax.append(model.maxtime())
    tmin = min(tmin)
    tmax = max(tmax)
    tgrid = linspace(tmin, tmax, int(tmax - tmin) + 1)
    mflux = model.bandflux(band, tgrid, zp=zp, zpsys=zpsys)
    return(tgrid,mflux)


def _unredden(color, bands, ebv, r_v=3.1):
    """
    (Private)
    Helper function that will get a magnitude correction from reddening coefficients using the CCM 1989 (and O'Donnell 1994)parameterization.
    """
    blue=bands[color[0]]
    red=bands[color[-1]]
    if not isinstance(blue,sncosmo.Bandpass):
        blue=sncosmo.get_bandpass(blue)
    if not isinstance(red,sncosmo.Bandpass):
        red=sncosmo.get_bandpass(red)


    blue_wave = array(blue.wave_eff, float)
    blue_a_lambda = _ccm_extinction(blue_wave, ebv, r_v)
    red_wave=array(red.wave_eff,float)
    red_a_lambda=_ccm_extinction(red_wave,ebv,r_v)

    return(blue_a_lambda-red_a_lambda)


def _ccm_extinction(wave, ebv, r_v=3.1):
    """
    (Private)
    Heler function for dereddening.

    """
    scalar = not iterable(wave)
    if scalar:
        wave = array([wave], float)
    else:
        wave = array(wave, float)

    x = 10000.0/wave
    npts = wave.size
    a = zeros(npts, float)
    b = zeros(npts, float)

    #Infrared
    good = where( (x > 0.3) & (x < 1.1) )
    a[good] = 0.574 * x[good]**(1.61)
    b[good] = -0.527 * x[good]**(1.61)

    # Optical & Near IR
    good = where( (x  >= 1.1) & (x < 3.3) )
    y = x[good] - 1.82

    c1 = array([ 1.0 , 0.104,   -0.609,    0.701,  1.137,
                    -1.718,   -0.827,    1.647, -0.505 ])
    c2 = array([ 0.0,  1.952,    2.908,   -3.989, -7.985,
                    11.102,    5.491,  -10.805,  3.347 ] )

    a[good] = polyval(c1[::-1], y)
    b[good] = polyval(c2[::-1], y)

    # Mid-UV
    good = where( (x >= 3.3) & (x < 8) )
    y = x[good]
    F_a = zeros(size(good),float)
    F_b = zeros(size(good),float)
    good1 = where( y > 5.9 )

    if size(good1) > 0:
        y1 = y[good1] - 5.9
        F_a[ good1] = -0.04473 * y1**2 - 0.009779 * y1**3
        F_b[ good1] =   0.2130 * y1**2  +  0.1207 * y1**3

    a[good] =  1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a
    b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b

    # Far-UV
    good = where( (x >= 8) & (x <= 11) )
    y = x[good] - 8.0
    c1 = [ -1.073, -0.628,  0.137, -0.070 ]
    c2 = [ 13.670,  4.257, -0.420,  0.374 ]
    a[good] = polyval(c1[::-1], y)
    b[good] = polyval(c2[::-1], y)

    # Defining the Extinction at each wavelength
    a_v = r_v * ebv
    a_lambda = a_v * (a + b/r_v)
    if scalar:
        a_lambda = a_lambda[0]
    return a_lambda


def curveToColor(lc,colors,snName,bandFit=None,snType='II',bandDict=_filters,color_bands=_opticalBands,zpsys='AB',model=None,singleBand=False,verbose=True, **kwargs):
    """
    Function takes a lightcurve file and creates a color table for it.

    :param lc: Name of lightcurve file you want to read, or astropy Table containing data
    :type lc: str or astropy.Table
    :param colors: Colors you want to calculate for the given SN (i.e U-B, r'-J)
    :type colors: str or list of strings
    :param snName: Name of the SN being used
    :type snName: str
    :param bandFit: If there is a specific band you would like to fit instead of default
    :type bandFit: str,optional
    :param snType: Classification of SN
    :type snType: str,optional
    :param bandDict: sncosmo bandpass for each band used in the fitting/table generation
    :type bandDict: dict,optional
    :param color_bands: bands making up the known component of chosen colors
    :type color_bands: list,optional
    :param zpsys: magnitude system (i.e. AB or Vega)
    :type zpsys: str,optional
    :param model: If there is a specific sncosmo model you would like to fit with, otherwise all models mathcing the SN classification will be tried
    :type model: str,optional
    :param singleBand: If you would like to only fit the bands in the color
    :type singleBand: Boolean,optional
    :param verbose: If you would like printing information to be turned on/off
    :type verbose: Boolean,optional
    :param kwargs: Catches all SNCOSMO fitting parameters here
    :returns: colorTable: Astropy Table object containing color information
    """
    if not isinstance(colors,(tuple,list)):
        colors=[colors]
    bands=append([col[:col.find('-')] for col in colors],[col[col.find('-')+1:] for col in colors])
    if bandDict is not None:
        for band in _filters:
            if band not in bandDict.keys() or band in bands or _filters[band] in bands:
                bandDict[band]=sncosmo.get_bandpass(_filters[band])
    else:
        bandDict={k:sncosmo.get_bandpass(k) for k in np.unique(lc[_get_default_prop_name('band')])}
        for band in _filters:
            if band not in bandDict.keys() or band in bands or _filters[band] in bands:
                bandDict[band]=sncosmo.get_bandpass(_filters[band])
    
    zpMag=sncosmo.get_magsystem(zpsys)

    if isinstance(lc,str):
        curve=_standardize(sncosmo.read_lc(lc))
    else:
        try:
            curve=_standardize(lc)
        except:
            raise RuntimeError("Can't understand your lightcurve.")
    if _get_default_prop_name('zpsys') not in curve.colnames:
        curve[_get_default_prop_name('zpsys')]=zpsys
    colorTable=Table(masked=True)
    colorTable.add_column(Column(data=[],name=_get_default_prop_name('time')))

    for band in bandDict:
        if not isinstance(bandDict[band],sncosmo.Bandpass):
            bandDict=_bandCheck(bandDict,band)
    t0=None
    if verbose:
        print('Getting best fit for: '+','.join(colors))

    args=[]
    for p in _fittingParams:
        args.append(kwargs.get(p,_fittingParams[p]))

        if p == 'method':
            try:
                importlib.import_module(_fittingPackages[args[-1]])
            except RuntimeError:
                sys.exit()
    for color in colors: #start looping through desired colors
        if bandDict[color[0]].wave_eff<_UVrightBound: #then extrapolating into the UV from optical
            if not bandFit:
                bandFit=color[-1]
            if singleBand:
                color_bands=[color[-1]]
            blue=curve[curve[_get_default_prop_name('band')]==color[:color.find('-')]] #curve on the blue side of current color
            red=curve[[x in color_bands for x in curve[_get_default_prop_name('band')]]] #curve on the red side of current color
        else: #must be extrapolating into the IR
            if not bandFit:
                bandFit=color[0]
            if singleBand:
                color_bands=[color[0]]
            blue=curve[[x in color_bands for x in curve[_get_default_prop_name('band')]]]
            red=curve[curve[_get_default_prop_name('band')]==color[color.find('-')+1:]]
        if len(blue)==0 or len(red)==0:
            if verbose:
                print('Asked for color %s-%s but missing necessary band(s)'%(color[:color.find('-')],color[color.find('-')+1:]))
            bandFit=None
            continue

        btemp=[bandDict[blue[_get_default_prop_name('band')][i]].name for i in range(len(blue))]
        rtemp=[bandDict[red[_get_default_prop_name('band')][i]].name for i in range(len(red))]
        blue.remove_column(_get_default_prop_name('band'))
        blue[_get_default_prop_name('band')]=btemp
        red.remove_column(_get_default_prop_name('band'))
        red[_get_default_prop_name('band')]=rtemp
        #now make sure we have zero-points and fluxes for everything
        if _get_default_prop_name('zp') not in blue.colnames:
            blue[_get_default_prop_name('zp')]=[zpMag.band_flux_to_mag(1/sncosmo.constants.HC_ERG_AA,blue[_get_default_prop_name('band')][i]) for i in range(len(blue))]
        if _get_default_prop_name('zp') not in red.colnames:
            red[_get_default_prop_name('zp')]=[zpMag.band_flux_to_mag(1/sncosmo.constants.HC_ERG_AA,red[_get_default_prop_name('band')][i]) for i in range(len(red))]
        if _get_default_prop_name('flux') not in blue.colnames:
            blue=mag_to_flux(blue,bandDict,zpsys)
        if _get_default_prop_name('flux') not in red.colnames:
            red=mag_to_flux(red,bandDict,zpsys)

        if not t0: #this just ensures we only run the fitting once
            if not model:
                if verbose:
                    print('No model provided, running series of models.')

                mod,types=loadtxt(os.path.join(__dir__,'data','sncosmo','models.ref'),dtype='str',unpack=True)
                modDict={mod[i]:types[i] for i in range(len(mod))}
                if snType!='Ia':
                    mods = [x for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and modDict[x[0]][:len(snType)]==snType]
                elif snType=='Ia':
                    mods = [x for x in sncosmo.models._SOURCES._loaders.keys() if 'salt2' in x[0]]

                mods = {x[0] if isinstance(x,(tuple,list)) else x for x in mods}
                if  bandFit==color[0] or len(blue)>len(red):
                    args[0]=blue
                    if len(blue)>60:
                        fits=[]
                        for mod in mods:
                            fits.append(_snFit([mod]+args))
                    else:
                        fits=pyParz.foreach(mods,_snFit,args)
                    fitted=blue
                    notFitted=red
                    fit=color[0]
                elif bandFit==color[-1] or len(blue)<len(red):
                    args[0]=red
                    '''
                    data,temp= sncosmo_fitting.cut_bands(photometric_data(red), sncosmo.Model(tempMod),
                                                             z_bounds=all_bounds.get('z', None),
                                                             warn=True)
                    print(data.fluxerr)
                    cov = diag(data.fluxerr**2) if data.fluxcov is None else data.fluxcov
                    invcov = linalg.pinv(cov)

                    args.append(invcov)
                    
                    sys.exit()
                    
                    fits=pyParz.foreach(mods,_snFit,args)
                    args.pop()
                    '''
                    if len(red)>60:
                        fits=[]
                        for mod in mods:
                            fits.append(_snFit([mod]+args))
                    else:
                        fits=pyParz.foreach(mods,_snFit,args)
                    fitted=red
                    notFitted=blue
                    fit=color[-1]
                else:
                    raise RuntimeError('Neither band "%s" nor band "%s" has more points, and you have not specified which to fit.'%(color[0],color[-1]))
                bestChisq=inf
                for f in fits:
                    if f:
                        res,mod=f
                        if res.chisq <bestChisq:
                            bestChisq=res.chisq
                            bestFit=mod
                            bestRes=res
                if verbose:
                    print('Best fit model is "%s", with a Chi-squared of %f'%(bestFit._source.name,bestChisq))
            elif bandFit==color[0] or len(blue)>len(red):
                args[0]=blue
                bestRes,bestFit=_snFit(append(model,args))
                fitted=blue
                notFitted=red
                fit=color[0]
                if verbose:
                    print('The model you chose (%s) completed with a Chi-squared of %f'%(model,bestRes.chisq))
            elif bandFit==color[-1] or len(blue)<len(red):
                args[0]=red
                bestRes,bestFit=_snFit(append(model,args))
                fitted=red
                notFitted=blue
                fit=color[-1]
                if verbose:
                    print('The model you chose (%s) completed with a Chi-squared of %f'%(model,bestRes.chisq))
            else:
                raise RuntimeError('Neither band "%s" nor band "%s" has more points, and you have not specified which to fit.'%(color[0],color[-1]))

            t0=_getBandMaxTime(bestFit,fitted,bandDict,'B',zpMag.band_flux_to_mag(1,bandDict['B']),zpsys)
            if len(t0)==1:
                t0=t0[0]
            else:
                raise RuntimeError('Multiple global maxima in best fit.')
        else:
            if len(blue)>len(red) or bandFit==color[0]:
                fitted=blue
                notFitted=red
                fit=color[0]
            elif len(blue)<len(red) or bandFit==color[-1]:
                fitted=red
                notFitted=blue
                fit=color[-1]
            else:
                raise RuntimeError('Neither band "%s" nor band "%s" has more points, and you have not specified which to fit.'%(color[0],color[-1]))
        #return (bestFit,bestRes,t0,fitted,notFitted)

        tGrid,bestMag=_snmodel_to_mag(bestFit,fitted,zpsys,bandDict[fit])

        ugrid,UMagErr,lgrid,LMagErr=_getErrorFromModel(append([bestFit._source.name,fitted],args[1:]),zpsys,bandDict[fit])
        tempTable=Table([tGrid-t0,bestMag,bestMag*.1],names=(_get_default_prop_name('time'),_get_default_prop_name('mag'),_get_default_prop_name('magerr')))#****RIGHT NOW THE ERROR IS JUST SET TO 10%*****
        interpFunc=scint.interp1d(array(tempTable[_get_default_prop_name('time')]),array(tempTable[_get_default_prop_name('mag')]))
        minterp=interpFunc(array(notFitted[_get_default_prop_name('time')]-t0))
        interpFunc=scint.interp1d(ugrid-t0,UMagErr)
        uinterp=interpFunc(array(notFitted[_get_default_prop_name('time')]-t0))
        interpFunc=scint.interp1d(lgrid-t0,LMagErr)
        linterp=interpFunc(array(notFitted[_get_default_prop_name('time')]-t0))
        magerr=mean([minterp-uinterp,linterp-minterp],axis=0)


        for i in range(len(minterp)):
            colorTable.add_row(append(notFitted[_get_default_prop_name('time')][i]-t0,[1 for j in range(len(colorTable.colnames)-1)]),mask=[True if j>0 else False for j in range(len(colorTable.colnames))])
        if fit==color[0]:
            colorTable[color]=MaskedColumn(append([1 for j in range(len(colorTable)-len(minterp))],minterp-array(notFitted[_get_default_prop_name('mag')])),mask=[True if j<(len(colorTable)-len(minterp)) else False for j in range(len(colorTable))])
        else:
            colorTable[color]=MaskedColumn(append([1 for j in range(len(colorTable)-len(minterp))],array(notFitted[_get_default_prop_name('mag')])-minterp),mask=[True if j<(len(colorTable)-len(minterp)) else False for j in range(len(colorTable))])
        colorTable[color[0]+color[-1]+'_err']=MaskedColumn(append([1 for j in range(len(colorTable)-len(magerr))],magerr+array(notFitted[_get_default_prop_name('magerr')])),mask=[True if j<(len(colorTable)-len(magerr)) else False for j in range(len(colorTable))])

        bandFit=None
    colorTable['SN']=snName
    colorTable.sort(_get_default_prop_name('time'))
    return(colorTable)




def mag_to_flux(table,bandDict,zpsys='AB'):
    """
    Accepts an astropy table of magnitude data and does the conversion to fluxes (flux in ergs/s/cm^2/AA)

    :param table: Table containing mag, mag error error, and band columns
    :type table: astropy.Table
    :param bandDict: translates band to sncosmo.Bandpass object (i.e. 'U'-->bessellux)
    :type bandDict: dict
    :param zpsys: magnitude system
    :type zpsys: str,optional
    :returns: astropy.Table object with flux and flux error added (flux in ergs/s/cm^2/AA)
    """
    ms=sncosmo.get_magsystem(zpsys)

    table[_get_default_prop_name('flux')]=asarray(list(map(lambda x,y: ms.band_mag_to_flux(x,y)*sncosmo.constants.HC_ERG_AA,asarray(table[_get_default_prop_name('mag')]),
                                                         table[_get_default_prop_name('band')])))
    table[_get_default_prop_name('fluxerr')] = asarray(
        list(map(lambda x, y: x * y / (2.5 * log10(e)), table[_get_default_prop_name('magerr')],
            table[_get_default_prop_name('flux')])))
    return(table)

def flux_to_mag(table,bandDict,zpsys='AB'):
    """
    Accepts an astropy table of flux data and does the conversion to mags (flux in ergs/s/cm^2/AA)

    :param table: Table containing flux, flux error error, and band columns
    :type table: astropy.Table
    :param bandDict: translates band to sncosmo.Bandpass object (i.e. 'U'-->bessellux)
    :type bandDict: dict
    :param zpsys: magnitude system
    :type zpsys: str,optional
    :returns: astropy.Table object with mag and mag error added
    """
    ms=sncosmo.get_magsystem(zpsys)
    table[_get_default_prop_name('mag')] = asarray(map(lambda x, y: ms.band_flux_to_mag(x/sncosmo.constants.HC_ERG_AA,y), table[_get_default_prop_name('flux')],
                                                          bandDict[table[_get_default_prop_name('band')]]))
    table[_get_default_prop_name('magerr')] = asarray(map(lambda x, y: 2.5 * log10(e) * y / x, table[_get_default_prop_name('flux')],
                                                             table[_get_default_prop_name('fluxerr')]))
    return(table)

def colorTableCombine(tableList):
    """
    Takes a list of color tables generated from curveToTable and returns a concatenated table.

    :param tableList: List of astropy.Table objects
    :type tableList: list
    :returns: Concatenated astropy.Table object
    """
    result=None
    for table in tableList:
        if result:
            result=vstack([result,table])
        else:
            result=table
    import math
    remove=[]
    i=0

    while i<len(result)-1:
        i+=1
        start=i-1
        while i<=len(result)-1 and math.fabs(result[_get_default_prop_name('time')][i]-result[_get_default_prop_name('time')][start])<.01:
            remove.append(i)
            for col in [x for x in result.colnames if x != _get_default_prop_name('time') and x!='SN']:
                    if result[i][col] and not result[start][col]:
                        result[start][col]=result[i][col]
                    elif result[i][col] and result[start][col]:
                        result[start][col]=average([array(result[start][col]),array(result[i][col])])

            i+=1

    result.remove_rows(remove)
    result.sort(_get_default_prop_name('time'))
    return(result)
