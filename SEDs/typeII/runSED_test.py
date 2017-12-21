from __future__ import print_function
#from allColors import *
from tempCustom import createSncosmoSED,plotSncosmoSED
import os,sncosmo,sys
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
from scipy import interpolate as scint
def _getZP(band,zpsys):
	"""
	Helper function that calculates the zero-point of a given band
	:param band: filename of band transmission function to get zero-point of 
	"""
	mag=sncosmo.get_magsystem(zpsys)
	return(mag.band_flux_to_mag(1,band))

def mag_to_flux(mag,band,zpsys):
	return(10**(-.4*(mag-_getZP(band,zpsys))))

def flux_to_mag(flux,band,zpsys):
	return(-2.5*np.log10(flux)+_getZP(band,zpsys))

import os
for f in ['j','h','ks']:
    wave,trans=np.loadtxt(os.path.join('/Users','jpierel','rodney','snsedextend','filters',f+'.dat'),unpack=True)
    wave*=10000
    sncosmo.registry.register(sncosmo.Bandpass(wave,trans,name='paritel::'+f),force=True)
def testSED(sedFile,amplitude=1.0,phase=0,minwave=2000.0,maxwave=20000.0,wavestep=500,showfig=True,outFilename="mySED.pdf",savefig=False,fontsize=18):
	sed=createSncosmoSED(sedFile) #should come back as photons not ergs
	fig=plotSncosmoSED(sed,phase=phase,minwave=minwave,maxwave=maxwave,showfig=showfig,outFilename=outFilename)
	ax=fig.gca()

	model=sncosmo.Model(sed)
	#print(model.flux(phase,sncosmo.get_bandpass('sdss::r').wave[0]),model.bandflux('sdss::r',phase,zp=_getZP('sdss::r','vega'),zpsys='vega'))
	#vcolors=np.array(ascii.read('vColors.dat')['V-r'])
	#phases=np.array(ascii.read('vColors.dat')['time'])
	#func=scint.interp1d(phases,vcolors)
	#vcolor=func(phase)
	#otherColors=ascii.read('allIIColors.dat')
	#jcolors=otherColors['r-H']
	#phases2=otherColors['time'][~jcolors.mask]
	
	#jcolors=jcolors[~jcolors.mask]
	#print(np.array(jcolors))
	#print()
	#print(phases2)
	#func2=scint.interp1d(phases2,jcolors)
	#jcolor=func2(phase)
	#color1=model.bandmag('bessellv','vega',phase)-model.bandmag('sdss::r','vega',phase)
	#vBand=flux_to_mag(model.flux(phase,(sncosmo.get_bandpass('bessellv').wave[0]+sncosmo.get_bandpass('bessellv').wave[-1])/2),'bessellv','vega')#model.bandmag('bessellv','vega',phase)
	#print(mag_to_flux(vBand,'bessellv','vega'))

	#vBand=model.bandmag('bessellv','vega',phase)
	#rBand=model.bandmag('sdss::r','vega',phase)
	#print(model.bandflux('bessellv',phase,_getZP('bessellv','vega'),'vega'))
	#print(model.bandmag('paritel::ks','vega',phase))
	print(model.color('bessellux','bessellb','vega',phase))
	sys.exit()
	#print(sncosmo.constants.HC_ERG_AA)
	#print(model.color('bessellv','sdss::r','vega',0))
	rBand=vBand/vcolor
	c_lambda=(vBand/(sncosmo.get_bandpass('bessellv').wave[-1]-sncosmo.get_bandpass('bessellv').wave[0]))*((sncosmo.get_bandpass('sdss::r').wave[-1]-sncosmo.get_bandpass('sdss::r').wave[0])/rBand)
	#r_lambda=
	print(rBand)
	#rBand=10**(-.4*(rBand.astype(float)-_getZP('sdss::r','vega')))
	ax.scatter((sncosmo.get_bandpass('sdss::r').wave[0]+sncosmo.get_bandpass('sdss::r').wave[-1])/2,rBand)
	plt.show()
	sys.exit()
	
	ax=fig.gca()
	#ax.text(sncosmo.get_bandpass('bessellv').wave[0],.2*ax.get_ylim()[1],'V-r={:3.3f}'.format(color))
	
	print(sncosmo.get_magsystem('vega').band_mag_to_flux(model.bandmag('bessellv','vega',phase),'bessellv'))
	ax.plot([sncosmo.get_bandpass('bessellv').wave[0],sncosmo.get_bandpass('sdss::r').wave[0]],[sncosmo.get_magsystem('vega').band_flux_to_mag(model.bandmag('bessellv','vega',phase),'bessellv'),sncosmo.get_magsystem('vega').band_flux_to_mag(model.bandmag('sdss::r','vega',phase),'sdss::r')])
	plt.savefig(outFilename,format='pdf')


#models.py->_bandflux_single....constants.py

def main():
	#filename='SDSS-018700.SED'
	filename='SDSS-000018.SED'
	testSED(filename,phase=0.96319999999999995,minwave=2500,maxwave=19600.0,showfig=True,outFilename=os.path.join('test_output',filename[:-4]+'.pdf'))


if __name__ == "__main__":
	main()