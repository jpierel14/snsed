import glob,sncosmo,sys,os
import numpy as np
from copy import deepcopy

files=glob.glob('lc_*')
remove=[]

sne,peaks=np.loadtxt('timeOfPeak.dat',dtype='str',unpack=True)
peaks={sne[i]:float(peaks[i]) for i in range(len(sne))}

for f in files:
	if f not in ['lc_2006jc.dat']:
		continue
	maxBand=None
	maxBand2=None
	lc=sncosmo.read_lc(f)
	for i in range(len(lc)):
		if lc['Band'][i]=='Ks':
			lc['Band'][i]='K'
	#print(lc)
	
	if f[:-4] in ['lc_2008aq']:
		lc=lc[lc['MJD']<peaks[f[:-4]]+40]
		lc=lc[lc['MJD']>peaks[f[:-4]]-40]
	elif f[:-4] in ['lc_2006jc']:
		lc=lc[lc['MJD']<peaks[f[:-4]]+30]
		lc=lc[lc['MJD']>peaks[f[:-4]]-30]
	elif f[:-4] in ['lc_2009iz']:
		lc=lc[lc['MJD']<peaks[f[:-4]]+20]
		lc=lc[lc['MJD']>peaks[f[:-4]]-20]
	else:
		lc=lc[lc['MJD']<peaks[f[:-4]]+50]
		lc=lc[lc['MJD']>peaks[f[:-4]]-50]
	'''
	
	elif f[:-4] in ['lc_2004aw','lc_2007uy','lc_2007gr']:
		lc=lc[lc['MJD']<peaks[f[:-4]]+25]
		lc=lc[lc['MJD']>peaks[f[:-4]]-25]
	elif f[:-4] in ['lc_2005az','lc_2005kl','lc_2006fo','lc_2006jc','lc_2009iz']:
		lc=lc[lc['MJD']<peaks[f[:-4]]+20]
		lc=lc[lc['MJD']>peaks[f[:-4]]-20]
	elif f[:-4] in ['lc_2005bf']:
		lc=lc[lc['MJD']<peaks[f[:-4]]+15]
		lc=lc[lc['MJD']>peaks[f[:-4]]-15]
	'''
	
	
	if 'U' not in lc['Band'] and 'u' not in lc['Band'] and 'J' not in lc['Band'] and 'H' not in lc['Band'] and 'Ks' not in lc['Band']:
		#os.remove(f)
		continue
	elif 'U' in lc['Band'] or 'u' in lc['Band']:
		if 'U' in lc['Band']:
			blue='U'
		else:
			blue='u'
		num=0
		maxBand=None
		maxBand2=None
		for Band in [x for x in np.unique(lc['Band']) if x in ['B','V','R','g','r']]:
			temp=len(lc[lc['Band']==Band])
			if temp>num:
				maxBand=Band
				num=temp
		red=maxBand
		if 'K' in lc['Band'] or 'J' in lc['Band'] or 'H' in lc['Band']:
			num=0
			maxBand=None
			maxBand2=None	
			for Band in [x for x in np.unique(lc['Band']) if x in ['B','V','R','g','r']]:
				temp=len(lc[lc['Band']==Band])
				if temp>num:
					maxBand=Band
					num=temp
			blue2=maxBand
			num=0
			maxBand=None
			maxBand2=None
			for Band in [x for x in np.unique(lc['Band']) if x in ['J','H','Ks']]:
				temp=len(lc[lc['Band']==Band])
				if temp>num:
					maxBand=Band
					num=temp
			red2=maxBand
			print(os.path.basename(f),blue+'-'+red,blue2+'-'+red2)
		else:
			print(os.path.basename(f),blue+'-'+red)

	else:
		num=0
		maxBand=None
		maxBand2=None
		for Band in [x for x in np.unique(lc['Band']) if x in ['B','V','R','g','r']]:
			temp=len(lc[lc['Band']==Band])
			if temp>num:
				maxBand=Band
				num=temp
		blue=maxBand
		num=0
		maxBand=None
		maxBand2=None
		for Band in [x for x in np.unique(lc['Band']) if x in ['J','H','Ks']]:
			temp=len(lc[lc['Band']==Band])
			if temp>num:
				maxBand2=Band
				num=temp
		red=maxBand2
		if red and blue:
			print(os.path.basename(f),blue+'-'+red)

	if not blue or not red:
		remove.append(f)
	sncosmo.write_lc(lc,f[:-4]+'_clipped.dat')
#print(remove)
