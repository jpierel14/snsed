import sncosmo
import glob,os,sys,multiprocessing,warnings
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

_filters={
	'U':'bessellux',
	'B':'bessellb',
	'V':'bessellv',
	'R':'bessellr',
	'I':'besselli',
	'J':'paritel::j',
	'H':'paritel::h',
	'K':'paritel::ks',
	'Ks':'paritel::ks',
	'u':'sdss::u',
	'g':'sdss::g',
	'r':'sdss::r',
	'i':'sdss::i',
	'Z':'sdss::z'
}

_needs_bounds={'z'}

method='minuit'
model=None
params=None
ignore=None
dust=None
effect_names=None
effect_frames=None
for f in ['j','h','ks']:
	wave,trans=np.loadtxt(os.path.join('/Users','jpierel','rodney','snsedextend','filters',f+'.dat'),unpack=True)
	wave*=10000
	sncosmo.registry.register(sncosmo.Bandpass(wave,trans,name='paritel::'+f),force=True)

sn,redshift=np.loadtxt('redshift.ref',dtype='str',unpack=True)
redshift={sn[i]:redshift[i] for i in range(len(sn))}

sn,types=np.loadtxt('type.ref',dtype='str',unpack=True)
typ={sn[i]:types[i] for i in range(len(sn))}

mod,types=np.loadtxt('models.ref',dtype='str',unpack=True)
modDict={mod[i]:types[i] for i in range(len(mod))}

sn,dusts=np.loadtxt('dust.ref',dtype='str',unpack=True)
dusts={sn[i]:dusts[i] for i in range(len(sn))}
			
def _mag_to_flux(table):

	table['flux'] = np.asarray(
		map(lambda x, y: 10 ** (-.4 * (x -y)), table['Mag'],
			table['zp']))
	table['fluxerr'] = np.asarray(
		map(lambda x, y: x * y / (2.5 * np.log10(np.e)), table['Magerr'],
			table['flux']))
	return table

def _snFit(args):
	#warnings.simplefilter('ignore')
	sn_func={'minuit': sncosmo.fit_lc, 'mcmc': sncosmo.mcmc_lc, 'nest': sncosmo.nest_lc}
	dust_dict={'SFD98Map':sncosmo.SFD98Map,'CCM89Dust':sncosmo.CCM89Dust,'OD94Dust':sncosmo.OD94Dust,'F99Dust':sncosmo.F99Dust}
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
	

	res,fit=sn_func[method](curve, model, params, bounds=bounds,verbose=False)
	return((res,fit))

def _snFitWrap(args):
    try:
        return(_snFit(args))
    except:
        return(None)

files=glob.glob('lc_*')
zpMag=sncosmo.get_magsystem('Vega')
sne=[]
times=[]
#plt.figure()
for f in files:
	if f not in ['lc_2006jc.dat']:
		continue
	print(f)
	t0=0
	lc=sncosmo.read_lc(f)
	if len(lc[lc['Band']=='U'])==0 and len(lc[lc['Band']=='u'])==0 and len(lc[lc['Band']=='J'])==0 and len(lc[lc['Band']=='H'])==0 and len(lc[lc['Band']=='Ks'])==0:
		os.remove(f)
		continue

	lc=lc[lc['Band']!='U']
	lc=lc[lc['Band']!='u']
	lc=lc[lc['Band']!='H']
	lc=lc[lc['Band']!='J']
	lc=lc[lc['Band']!='K']
	lc=lc[lc['Band']!='Ks']
	lc=lc[lc['Band']!='i']
	lc=lc[lc['Band']!='I']
	if not lc:
		os.remove(f)
		continue
	lc['zpsys']='Vega'
	lc['zp']=[zpMag.band_flux_to_mag(1,_filters[lc['Band'][i]]) for i in range(len(lc))]
	lc=_mag_to_flux(lc)
	lc['Band']=[_filters[lc['Band'][i]] for i in range(len(lc))]
	minchisq=np.inf
	fits=[]
	'''
	if f=='lc_2005kd.dat':
		lc=lc[lc['MJD']>np.min(lc['MJD'])-30]
		lc=lc[lc['MJD']<np.min(lc['MJD'])+30]
	'''
	print(np.min(lc['MJD']))
	constants={'z':redshift[f[:-4]],'hostr_v':3.1,'mwr_v':3.1,'mwebv':dusts[f[:-4]]}
	bounds={'t0':(np.min(lc['MJD'])-20,np.min(lc['MJD'])+20),'hostebv':(-1,1)}
	tempType=typ[f[:-4]]
	if '/' in tempType:
		temp1=tempType[:tempType.rfind('/')]
		temp2=tempType[tempType.rfind('/')+1:]
		mods= [x for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and (modDict[x[0]][:len(temp1)]==temp1 or modDict[x[0]][:len(temp2)]==temp2)]
	else:
		mods = [x for x in sncosmo.models._SOURCES._loaders.keys() if x[0] in modDict.keys() and modDict[x[0]][:len(typ[f[:-4]])]==typ[f[:-4]]]
	mods = {x[0] if isinstance(x,(tuple,list)) else x for x in mods}
	dust='CCM89Dust'
	effect_frames=['rest','obs']
	effect_names=['host','mw']
	'''
	if len(mods)>3 and f not in ['lc_2007aa.dat']:
		p = Pool(processes=multiprocessing.cpu_count())
		for x in p.imap_unordered(_snFitWrap,[(x,lc,method,params,bounds,ignore,constants,dust,effect_names,effect_frames) for x in mods]):
			fits.append(x)
		p.close()
	else:
	'''
	for mod in mods:
		#print(mod)
		fits.append(_snFit((mod,lc,method,params,bounds,ignore,constants,dust,effect_names,effect_frames)))
	'''
	for model in models[typ[f[:-4]]]:
		mod=sncosmo.Model(source=model)
		mod.set(z=redshift[f[:-4]])
		res,fit=sncosmo.fit_lc(lc,mod,[x for x in mod.param_names if x != 'z'],bounds={'t0':(np.min(lc['MJD'])-20,np.min(lc['MJD'])+20)})
	'''
	bestChisq=np.inf

	for fit in fits:
		if fit:
			res,mod=fit
			if res.chisq <bestChisq:
				bestChisq=res.chisq
				bestfit=mod
				bestres=res
	t0=bestres.parameters[bestres.param_names.index('t0')]
	#print(bestChisq,bestfit._source.name,t0,np.min(lc['MJD']))

	sncosmo.plot_lc(lc,model=bestfit,errors=bestres.errors)
	plt.show()
	sne.append(f)
	times.append(t0)
np.savetxt('sne.dat',sne,fmt='%s')
np.savetxt('peaks.dat',times)











