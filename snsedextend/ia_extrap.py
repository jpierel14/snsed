"""
2016   S. Rodney
Module for generating the SALT2 model data files with extrapolation to
IR wavelengths, guided by a sample of low-z Type Ia SN templates provided by
Arturo Avelino and Andy Friedman
"""

__all__=['WavelengthGrid','TimeSeriesGrid','extend_dispersion_files',
         'extend_variance_covariance', 'extend_info_file','extendIa']

from matplotlib import pyplot as pl, ticker

from scipy import integrate as scint, interpolate as scinterp, optimize as scopt
import sncosmo,os,glob
import numpy as np
from astropy.table import Table
from astropy.io import ascii
from scipy.interpolate import interp1d

from .utils import __dir__

_B_WAVELENGTH = 4302.57
_V_WAVELENGTH = 5428.55

__MODELDIR__ = os.path.join(__dir__,'data')

__SALT2_ERROR_FILES__=[
    'salt2_color_correction.dat',
    'salt2_color_dispersion.dat',
    'salt2_lc_dispersion_scaling.dat',
    'salt2_lc_relative_covariance_01.dat',
    'salt2_lc_relative_variance_0.dat',
    'salt2_lc_relative_variance_1.dat',
    'salt2_spec_covariance_01.dat',
    'salt2_spec_variance_0.dat',
    'salt2_spec_variance_1.dat'
]
class WavelengthGrid(object):
    """A class for handling a SN SED component (i.e., a single 2-D grid of
    wavelength and value (no variaition with phase). The input SEDfile should
     have two columns giving wavelength and value. """
    def __init__(self, datafile):
        self.datafile = datafile
        fin = open(self.datafile, 'r')
        all_lines = fin.readlines()
        self.headerlines = np.array([line for line in all_lines
                                     if len(line.strip()) > 0 and
                                     line.lstrip().startswith('#')])
        self.wave, self.value = self.read_wavelengthgrid_data()

    def read_wavelengthgrid_data(self) :
        """ read in the wavelength and value grids from the given data file
        (e.g. from a salt2_color_dispersion.dat file)
        :rtype: np.ndarray, np.ndarray
        :param file:
        :return phase, wave, flux: phase is a 1D array with each entry giving the
           rest-frame phase relative to B band max. wave and flux are 2D arrays,
           each with a 1D array for every day in the phase array.
        """
        w,f = np.loadtxt( self.datafile, unpack=True )
        return w, f

    def plot(self, fluxoffset=0, **kwargs):
        """ plot the SED template data from the given set of data arrays
        :param fluxoffset: additive offset to the flux (for separating multiple
        curves on the figure)
        :param kwargs: passed on to matplotlib.pyplot.plot()
        :return:
        """
        pl.plot(self.wave, self.value + fluxoffset, **kwargs)
        return

    def extrapolate_flatline(self, outfilename, refwave=8500,
                             maxwave=25000):
        """ use a super-crude flatline extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        w = self.wave
        v = self.value
        wavestep = w[1] - w[0]

        # redward flatline extrapolation from last point
        wavenew = np.arange(w[0], maxwave, wavestep)
        iref = np.argmin(np.abs(w-refwave))
        valnew = np.append(v[:iref],
                           np.zeros(len(wavenew)-len(w[:iref])) + v[iref])

        # append to the list of output data lines
        for j in range(len(wavenew)):
            outlines.append(" %12i  %12.7e\n" % (
                wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

    def extrapolate_flatline_uv(self, outfilename, refwave=3500,
                                minwave=1700,maxwave=25000):
        """ use a super-crude flatline extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        w = self.wave

        v = self.value
        wavestep = w[1] - w[0]

        # redward flatline extrapolation from last point
        wavenew = np.arange(minwave, maxwave, wavestep)
        iref = np.argmin(np.abs(w-refwave))
        valnew = np.append(np.zeros(len(np.arange(minwave,refwave,wavestep))) + v[iref],v[iref:])

        # append to the list of output data lines
        for j in range(len(wavenew)):
            outlines.append(" %12i  %12.7e\n" % (
                wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

    def extrapolate_linearfit(self, outfilename, refwave=8500, maxwave=25000):
        """ use a super-crude linear-fit extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        w = self.wave
        v = self.value
        wavestep = w[1] - w[0]

        # red end of the array, to be fit
        iref = np.argmin(np.abs(w-refwave))

        linfitres = scopt.curve_fit(linear, w[iref:], v[iref:])
        linfitA = linfitres[0][0]
        linfitB = linfitres[0][1]

        # redward flatline extrapolation from last point
        waveir = np.arange(w[-1]+wavestep, maxwave, wavestep)
        wavenew = np.append(w, waveir)
        valnew = np.append( v, linear(waveir, linfitA, linfitB))

        # append to the list of output data lines
        for j in range(len(wavenew)):
            outlines.append("%12i  %12.7e\n" % (
                wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

    def extrapolate_linearto0(self, outfilename, refwave=8500,
                              maxwave=25000, maxwaveval=0.0):
        """ use a super-crude linear-fit extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        w = self.wave
        v = self.value
        wavestep = w[1] - w[0]

        # red end of the array, to be fit
        iref = np.argmin(np.abs(w-refwave))

        linfitB = (maxwaveval-v[iref]) / (maxwave-refwave)
        linfitA = maxwaveval - linfitB * maxwave

        # redward flatline extrapolation from last point
        waveir = np.arange(refwave, maxwave, wavestep)
        wavenew = np.append(w[:iref], waveir)
        valnew = np.append( v[:iref],
                            linear(waveir, linfitA, linfitB))

        # append to the list of output data lines
        for j in range(len(wavenew)):
            outlines.append("%12i  %12.7e\n" % (
                wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return


class TimeSeriesGrid(object):
    """A class for handling a SN SED time series (i.e., a sequence of spectro-
     photometric template arrays at many phases. The input SEDfile should 
     have three columns giving phase, wavelength and flux. """
    def __init__(self, sedfile):
        self.sedfile = sedfile
        fin = open(self.sedfile,'r')
        all_lines = fin.readlines()
        self.headerlines = np.array([line for line in all_lines
                                     if len(line.strip()) > 0 and
                                     line.lstrip().startswith('#')])
        self.parse_lowz_template_header()
        self.phase, self.wave, self.value = self.read_timeseriesgrid_data()

        # check that all phases have the same number of wavelength steps
        # If any phases do not, then it is probably a duplicate
        # data entry, so we just trim it to the correct size.
        nwavestepslist = [len(w) for w in self.wave]
        if len(np.unique(nwavestepslist))>1:
            truensteps = np.min(nwavestepslist)
            ibad = np.where(nwavestepslist>truensteps)[0]
            for i in ibad:
                self.wave[i] = self.wave[i][:truensteps]
                self.value[i] = self.value[i][:truensteps]
            # igood = np.where(nwavestepslist==truensteps)[0]
            # self.phase = self.phase[igood]
            # self.wave = self.wave[igood]
            # self.value = self.value[igood]

        nwavestepslist = [len(w) for w in self.wave]
        assert (len(np.unique(nwavestepslist))==1)



    def parse_lowz_template_header(self):
        """ read in metadata from the sed data file header
        :return:
        """
        for hdrline in self.headerlines:
            hdrline = hdrline.strip()
            if '=' not in hdrline:
                continue
            hdrline = hdrline.strip('# ')
            key = hdrline.split('=')[0].strip()
            value = hdrline.split('=')[1].strip()
            if key not in ['name','survey']:
                value = float(value)
            self.__dict__[key] = value

    def read_timeseriesgrid_data(self) :
        """ read in the phase, wavelength and flux grids from the given SN SED
        template data file (e.g. from a SALT2 template0.dat file)
        :rtype: np.ndarray, np.ndarray, np.ndarray
        :param sedfile:
        :return phase, wave, flux: phase is a 1D array with each entry giving the
           rest-frame phase relative to B band max. wave and flux are 2D arrays,
           each with a 1D array for every day in the phase array.
        """
        p,w,f = np.loadtxt( self.sedfile, unpack=True )
        phaselist = np.unique( p )
        phasearray = np.array(phaselist)
        wavearray = np.array([ w[ np.where( p == day ) ] for day in phaselist ])
        fluxarray = np.array([ f[ np.where( p == day ) ] for day in phaselist ])
        return phasearray, wavearray, fluxarray

    def plot_at_phase(self, phase=0, fluxoffset=0, **kwargs):
        """ plot the SED template data from the given set of data arrays, only at
        the specified phase.
        :param phase: phase of the SED model to plot (rel. to B band max)
        :param fluxoffset: additive offset to the flux (for separating multiple
        curves on the figure)
        :param kwargs: passed on to matplotlib.pyplot.plot()
        :return:
        """
        ithisphase = np.where(np.abs(self.phase-phase)<1)[0]
        wave = self.wave[ithisphase, :][0]
        flux = self.value[ithisphase, :][0]
        pl.plot(wave, flux + fluxoffset, **kwargs)
        return wave, flux


    def extrapolate_flatline(self, outfilename, refwave=8500,
                             maxwave=25000):
        """ use a super-crude flatline extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        for iphase in range(len(self.phase)):
            thisphase = self.phase[iphase]
            w = self.wave[iphase]
            v = self.value[iphase]
            wavestep = w[1] - w[0]

            # redward flatline extrapolation from last point
            wavenew = np.arange(w[0], maxwave, wavestep)
            iref = np.argmin(np.abs(w-refwave))
            valnew = np.append(v[:iref],
                               np.zeros(len(wavenew)-len(w[:iref])) + v[iref])

            # append to the list of output data lines
            for j in range(len(wavenew)):
                outlines.append("%6.2f    %12i  %12.7e\n" % (
                    thisphase, wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

    def extrapolate_flatline_uv(self, outfilename, refwave=3500,
                                minwave=1700,maxwave=25000):
        """ use a super-crude flatline extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        for iphase in range(len(self.phase)):
            thisphase = self.phase[iphase]
            w = self.wave[iphase]

            v = self.value[iphase]
            wavestep = w[1] - w[0]

            # redward flatline extrapolation from last point
            wavenew = np.arange(minwave, maxwave, wavestep)
            iref = np.argmin(np.abs(w-refwave))
            valnew = np.append(np.zeros(len(np.arange(minwave,refwave,wavestep))) + v[iref],v[iref:])

            # append to the list of output data lines
            for j in range(len(wavenew)):
                outlines.append("%6.2f    %12i  %12.7e\n" % (
                    thisphase, wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

    def extrapolate_linearfit(self, outfilename, refwave=8500,
                              maxwave=25000, showplots=False):
        """ use a super-crude linear-fit extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        for iphase in range(len(self.phase)):
            thisphase = self.phase[iphase]
            w = self.wave[iphase]
            v = self.value[iphase]
            wavestep = w[1] - w[0]

            # red end of the array, to be fit
            iref = np.argmin(np.abs(w-refwave))

            linfitres = scopt.curve_fit(linear, w[iref:], v[iref:])
            linfitA = linfitres[0][0]
            linfitB = linfitres[0][1]

            # redward flatline extrapolation from last point
            waveir = np.arange(w[-1]+wavestep, maxwave, wavestep)
            wavenew = np.append(w, waveir)
            valnew = np.append( v, linear(waveir, linfitA, linfitB))

            # append to the list of output data lines
            for j in range(len(wavenew)):
                outlines.append("%6.2f    %12i  %12.7e\n" % (
                    thisphase, wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

    def extrapolate_linearto0(self, outfilename, refwave=8500,
                              maxwave=25000, maxwaveval=0):
        """ use a super-crude linear-fit extension of the value array
        to extrapolate this model component to the IR at all phases """
        outlines = []
        for iphase in range(len(self.phase)):
            thisphase = self.phase[iphase]
            w = self.wave[iphase]
            v = self.value[iphase]
            wavestep = w[1] - w[0]

            # red end of the array, to be fit
            iref = np.argmin(np.abs(w-refwave))

            linfitB = (maxwaveval-v[iref]) / (maxwave-refwave)
            linfitA = maxwaveval - linfitB * maxwave

            # redward flatline extrapolation from last point
            waveir = np.arange(refwave, maxwave, wavestep)
            wavenew = np.append(w[:iref], waveir)
            valnew = np.append( v[:iref],
                                linear(waveir, linfitA, linfitB))

            # append to the list of output data lines
            for j in range(len(wavenew)):
                outlines.append("%6.2f    %12i  %12.7e\n" % (
                    thisphase, wavenew[j], valnew[j]))

        # write it out to the new .dat file
        fout = open(outfilename, 'w')
        fout.writelines(outlines)
        fout.close()
        return

###UV Extrapolation functions

def _getFoleyData():
    foleyWave,fLambda,sLambda=np.loadtxt(os.path.join(__MODELDIR__,'salt2-uv','uvmodel.data'),unpack=True)
    return(foleyWave,fLambda,sLambda)

def _getSalt2Data():
    m0phase,m0wave,m0val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-4','salt2_template_0.dat'))
    m1phase,m1wave,m1val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-4','salt2_template_1.dat'))
    return(m0phase,m0wave,m0val,m1val)

def _getRodneyExtrap():
    m0phase,m0wave,m0val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-ir','salt2_template_0.dat'))
    m1phase,m1wave,m1val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-ir','salt2_template_1.dat'))
    return(m0phase,m0wave,m0val,m1val)

def _getUVExtrap():
    m0phase,m0wave,m0val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-uv','salt2_template_0.dat'))
    m1phase,m1wave,m1val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-uv','salt2_template_1.dat'))
    return(m0phase,m0wave,m0val,m1val)

def _getSalt2extended():
    m0phase,m0wave,m0val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-extended','salt2_template_0.dat'))
    m1phase,m1wave,m1val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-extended','salt2_template_1.dat'))
    return(m0phase,m0wave,m0val,m1val)

def createSalt2(x1Param):
    if not isinstance(x1Param,(list,tuple,np.ndarray)):
        x1Param=[x1Param]
    salt2=[]
    for x1 in x1Param:
        model=sncosmo.Model('salt2')
        model.set(c=0,z=0,x0=1,x1=x1)
        salt2.append(model)
    return(salt2)

def getDeltaM15(models):
    if not isinstance(models,(list,tuple,np.ndarray)):
        models=[models]
    deltaM15=[]
    bBand=sncosmo.get_bandpass('bessellb')
    for model in models:
        deltaM15.append(model.bandmag(bBand,'Vega',15)-model.bandmag(bBand,'Vega',0))
    return(deltaM15)


def deltaM15toX1(deltaM15):
    x1Params=np.linspace(-4,4,40)
    salt2Models=createSalt2(x1Params)
    deltaMs=getDeltaM15(salt2Models)
    result=Table(([model.get('x1') for model in salt2Models],deltaMs),names=('x1',"deltaM15"))
    interpFunc=interp1d(result['deltaM15'],result['x1'])
    return(interpFunc(deltaM15))




def foley2Salt(deltaM15,x0=1):
    x1=deltaM15toX1(deltaM15)
    m0phase,m0wave,m0val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-4','salt2_template_0.dat'))
    m1phase,m1wave,m1val=sncosmo.read_griddata_ascii(os.path.join(__MODELDIR__,'salt2-4','salt2_template_1.dat'))
    foleyWave,fLambda,sLambda=_getFoleyData()
    findex=list(foleyWave).index(3500)
    findex2=list(foleyWave).index(m0wave[0])
    totalFunc1=interp1d([-6,0],[0,1])
    totalFunc2=interp1d([0,11],[1,0])
    fweighting=np.append(np.append([1 for i in range(len(foleyWave[foleyWave<m0wave[0]]))],np.linspace(1,0,len(foleyWave[foleyWave>=m0wave[0]][:findex-findex2])-1)),[0])

    interpFunc=interp1d(foleyWave,fLambda)
    interpFunc2=interp1d(foleyWave[:findex],fweighting)
    interpFunc3=interp1d(foleyWave,sLambda)
    newFoleyWave=np.arange(foleyWave[0],foleyWave[findex],(m0wave[1]-m0wave[0]))
    newfLambdaVals=interpFunc(newFoleyWave)
    newSLambdaVals=interpFunc3(newFoleyWave)
    finalM0=np.ndarray((len(m0phase),len(newFoleyWave)))
    salt2x0=np.zeros(len(m0phase))
    finalM1=np.ndarray((len(m0phase),len(newFoleyWave)))
    finalFlux=np.ndarray((len(m0phase),len(newFoleyWave)))
    ind=list(newFoleyWave).index(m0wave[0])
    waveWeight=interpFunc2(newFoleyWave[ind:])
    foleyFlux=newfLambdaVals+newSLambdaVals*(deltaM15-1.1)

    for i in range(len(m0phase)):
        if m0phase[i]>=-5:
            if m0phase[i]<=0:
                phaseWeight=float(totalFunc1(m0phase[i]))
            elif m0phase[i]<=10:
                phaseWeight=float(totalFunc2(m0phase[i]))
            else:
                phaseWeight=0
        else:
            phaseWeight=0
        m0interpFunc=interp1d(m0wave,m0val[i])
        m1interpFunc=interp1d(m1wave,m1val[i])

        salt2x0[i]=(m0interpFunc(newFoleyWave[-1])+x1*m1interpFunc(newFoleyWave[-1]))/foleyFlux[-1] if foleyFlux[-1]>0 else 1#this is the only line i'm not sure about

        finalM0[i][:ind]=phaseWeight*salt2x0[i]*(newfLambdaVals[:ind])
        finalM1[i][:ind]=phaseWeight*salt2x0[i]*(newSLambdaVals[:ind]*(deltaM15-1.1))/x1 if x1 !=0 else 0
        finalM0[i][ind:]=(1-phaseWeight)*m0interpFunc(newFoleyWave[ind:])+phaseWeight*((1-waveWeight)*m0interpFunc(newFoleyWave[ind:])+
                                                                                       waveWeight*salt2x0[i]*newfLambdaVals[ind:])
        finalM1[i][ind:]=(1-phaseWeight)*m1interpFunc(newFoleyWave[ind:])+phaseWeight*((1-waveWeight)*m1interpFunc(newFoleyWave[ind:])+
                                                                                       waveWeight*salt2x0[i]*newSLambdaVals[ind:]*(deltaM15-1.1)/x1)
       
        actualSalt2=_getSalt2Flux(x0,x1,m0interpFunc(newFoleyWave[ind:]),m1interpFunc(newFoleyWave[ind:]))

        if phaseWeight>0:
            finalFlux[i]=x0*(finalM0[i]+x1*finalM1[i])
        else:
            finalFlux[i][ind:]=actualSalt2

    return(m0phase,newFoleyWave,finalM0,finalM1,x1,finalFlux,salt2x0)

def _getSalt2Flux(x0,x1,m0,m1):
    return(x0*(m0+x1*m1))

def outputM0M1(wave,x0,m0,m1,m0Filename="salt2_template_0.dat",m1Filename="salt2_template_1.dat"):
    phase,salt2wave,salt2m0,salt2m1=_getSalt2Data()
    #newwave=np.unique(np.append(wave,salt2wave))
    #ind=list(wave).index(salt2wave[0])
    fout=open(os.path.join(__MODELDIR__,'salt2-uv',os.path.basename(m0Filename)),'w')
    fout2=open(os.path.join(__MODELDIR__,'salt2-uv',os.path.basename(m1Filename)),'w')
    for i in range(len(phase)):
        for j in range(len(wave)):
            fout.write("%5.6f  %10.6f  %12.7e \n"%(phase[i],wave[j],m0[i][j]))
            fout2.write("%5.6f  %10.6f  %12.7e \n"%(phase[i],wave[j],m1[i][j]))
    fout.close()
    fout2.close()
    return

def combineUV_IR():
    irphase,irwave,irm0,irm1=_getRodneyExtrap()
    uvphase,uvwave,uvm0,uvm1=_getUVExtrap()
    newwave=np.unique(np.append(uvwave,irwave))
    ind=list(newwave).index(irwave[0])
    fout=open(os.path.join(__MODELDIR__,'salt2-extended',os.path.basename('salt2_template_0.dat')),'w')
    fout2=open(os.path.join(__MODELDIR__,'salt2-extended',os.path.basename('salt2_template_1.dat')),'w')
    for i in range(len(irphase)):
        for j in range(len(newwave)):
            if newwave[j]>uvwave[-1]:
                fout.write("%5.6f  %10.6f  %12.7e \n"%(irphase[i],newwave[j],irm0[i][j-ind]))
                fout2.write("%5.6f  %10.6f  %12.7e \n"%(irphase[i],newwave[j],irm1[i][j-ind]))
            else:
                fout.write("%5.6f  %10.6f  %12.7e \n"%(irphase[i],newwave[j],uvm0[i][j]))
                fout2.write("%5.6f  %10.6f  %12.7e \n"%(irphase[i],newwave[j],uvm1[i][j]))
    fout.close()
    fout2.close()


def extendIa(deltaM15=1.1,showplots=False):
    phase,wave,m0,m1,x1,salt2Flux,x0=foley2Salt(deltaM15)
    outputM0M1(wave,x0,m0,m1)
    extend_dispersion_files(showplots=showplots)
    extend_variance_covariance(showplots=showplots)
    extend_info_file(showplots=showplots)
    fit_salt2colorlaw_to_ccm(c=0.1, wmin=2800, wmax=24990, wjoin=6750, salt2_colorlaw_range=[2800,10000], showfit=showplots)
    seddir = os.path.join(__MODELDIR__,'salt2-ir','salt2irdata','SNeIa_seds_v2')
    subdirlist = os.listdir(seddir)

    sedfilelist = []
    for subdir in subdirlist:
        fullpath = os.path.join(seddir,subdir)
        filenames = os.listdir(fullpath)
        sedfilelist += [os.path.join(fullpath, f) for f in filenames]

    deredden_template_list(sedfilelist, showfig=showplots, savefig=False)
    modeldict = load_sncosmo_models(x1min=-1, x1max=1, cmin=-0.1, cmax=0.3)

    extend_template0_ir(modeldict, verbose=False)
    if showplots:
        plot_extended_template0()
        pl.show()

    extend_template1_ir(modeldict=modeldict, showplots=showplots)
    combineUV_IR()


###IR Extrapolation functions


def createSncosmoSource(filename):
    """Create an sncosmo Source from a suitable SED file containing a 
    time series of spectrophotometric templates. """
    phase1, wave1, flux1 = sncosmo.read_griddata_ascii(filename)
    snname = os.path.basename(filename).split('.')[0]
    source = sncosmo.TimeSeriesSource(phase1, wave1, flux1, name=snname,
                                      zero_before=True)
    return source


def load_sncosmo_models(modeldir=__MODELDIR__,
                        salt2indir='salt2-4', salt2outdir='salt2-ir',
                        lowzsubdir='lowzIa',
                        x1min=-1, x1max=1, cmin=-0.1, cmax=0.3):
    """
    Load all the lowz template SEDs from sncosmo, along with the original
    salt2 model and the Hsiao Type Ia model.
    :param datadir:
    :return:
    """
    # read in all the low-z mangled SEDs
    lowzmodeldir = os.path.join(modeldir, 'salt2-ir','salt2irdata',lowzsubdir)
    lowzsedfilelist = glob.glob("%s/*.dat" % lowzmodeldir)
    modeldict = {}
    for lowzsedfile in lowzsedfilelist:
        source = createSncosmoSource(lowzsedfile)
        model = sncosmo.Model(source)

        # read the header info of the sed template file to
        # determine the x1,c, Delta, Av, z  and
        # store these as properties of the sncosmo Model object
        fin = open(lowzsedfile,'r')
        all_lines = fin.readlines()
        for hdrline in all_lines:
            hdrline = hdrline.strip()
            if len(hdrline)>0 and not hdrline.startswith('#'):
                break
            if '=' not in hdrline:
                continue
            hdrline = hdrline.strip('# ')
            key = hdrline.split('=')[0].strip()
            value = hdrline.split('=')[1].strip()
            if key not in ['name','survey']:
                value = float(value)
            model.__dict__[key] = value
        if model.salt2x1 < x1min:
            continue
        if model.salt2x1 > x1max:
            continue
        if model.salt2c < cmin:
            continue
        if model.salt2c > cmax:
            continue
        modeldict[source.name] = model

    # read in the original and the revised salt2 models:
    modeldict['hsiao'] = sncosmo.Model('hsiao')

    salt2modeldir = os.path.join(modeldir, salt2indir)
    salt2 = sncosmo.models.SALT2Source(modeldir=salt2modeldir, name='salt2')
    modeldict['salt2'] = salt2

    salt2irmodeldir = os.path.join(modeldir, salt2outdir)
    try:
        salt2ir = sncosmo.models.SALT2Source(modeldir=salt2irmodeldir,
                                            name='salt2-ir')
        modeldict['salt2-ir'] = salt2ir
    except:
        pass

    return modeldict


def plot_template0_data(modeldict=None, phase=0, x1=0, c=0):
    """
    Load all the lowz template SEDs from sncosmo

    Plot the template0 data at the given phase with

    :param datadir:
    :return:
    """
    if modeldict == None:
        modeldict = load_sncosmo_models()

    salt2mod = modeldict['salt2']
    salt2mod.set(x1=x1, c=c)

    salt2irmod = modeldict['salt2-ir']
    salt2irmod.set(x1=x1, c=c)

    wave0 = np.arange(salt2mod.minwave(),salt2mod.maxwave(),10)
    wavelowz = np.arange(salt2mod.minwave(),18000,10)
    waveir = np.arange(2000, salt2irmod.maxwave())

    # plot all the low-z mangled SEDs at the given phase
    normwavemin = 3500
    normwavemax = 8500

    # plot all the low-z mangled SEDs at the given phase
    pl.clf()
    fluxlowzarray = []
    for name in modeldict.keys():
        if not name.startswith('sn'):
            continue
        lowzsn = modeldict[name]
        fluxlowz = lowzsn.flux(phase, wavelowz)
        # iwavemax = np.where(wavelowz>=np.max(wave0))[0][0]
        # fluxlowz_norm = scint.trapz(fluxlowz[:iwavemax], wavelowz[:iwavemax])
        iwavemin = np.argmin(np.abs(wavelowz-normwavemin))
        iwavemax = np.argmin(np.abs(wavelowz-normwavemax))
        fluxlowz_norm = scint.trapz(fluxlowz[iwavemin:iwavemax],
                                    wavelowz[iwavemin:iwavemax])
        pl.plot(wavelowz, fluxlowz/fluxlowz_norm,
                color='b',ls='-', alpha=0.1, lw=2 )
        # pl.plot(waveir, snflux/snflux.sum(), 'b-', alpha=0.3, lw=1 )
        fluxlowzarray.append(fluxlowz/fluxlowz_norm)

    # plot the normalized template0 data from the original SALT2 model
    # at the given phase:
    flux0 = salt2mod.flux(phase, wave0)
    iwavemin0 = np.argmin(np.abs(wave0 - normwavemin))
    iwavemax0 = np.argmin(np.abs(wave0 - normwavemax))
    flux0_norm = scint.trapz(flux0[iwavemin0:iwavemax0],
                             wave0[iwavemin0:iwavemax0])
    pl.plot(wave0, flux0/flux0_norm, color='0.5', ls='-', lw=3.5,
            label='SALT2' )

    # plot the normalized template0 flux from the modified SALT2-IR model
    # at the given phase:
    fluxir = salt2irmod.flux(phase, waveir)
    iwaveminir = np.argmin(np.abs(waveir - normwavemin))
    iwavemaxir = np.argmin(np.abs(waveir - normwavemax))
    fluxir_norm = scint.trapz(fluxir[iwaveminir:iwavemaxir],
                              waveir[iwaveminir:iwavemaxir])
    pl.plot(waveir, fluxir/fluxir_norm, marker=' ', color='darkorange',
            ls='-', lw=1.2, label='SALT2-IR', zorder=90)

    ax = pl.gca()
    ax.set_xlabel(r'Wavelength ($\rm{\AA}$)',fontsize=14)
    ax.set_ylabel('SALT2 M0 value (arbitrary units)',fontsize=14)
    ax.set_title('Template SEDs Used in SALT2-IR Extrapolation',fontsize=15)
    pl.setp(ax.get_yticklabels(), visible=False)
    # return temp0, modeldict
    pl.legend(loc='upper right')
    ax.set_xlim(5500,15000)
    ax.set_ylim(0,0.00021)


def check_for_salt2_fits(sedfiledir='lowzIa',
                         fitresfilename='lowz_salt2.fitres'):
    """ Check how many of our Low-z Ia sample have SALT2 fit parameters in
    the .fitres file from D.Scolnic.
    """
    # get a list of SN names from te SED file directory.
    sedfilelist = glob.glob("%s/*.dat"%sedfiledir)
    snnamelist = [os.path.basename(s).split('_')[0].lower() for s in sedfilelist]

    # read in the low-z SN salt2 fit parameters from the file provided by Dan
    salt2fitfile = os.path.join(__MODELDIR__,'salt2-ir','salt2irdata', fitresfilename)
    salt2fitdata = ascii.read(salt2fitfile,
                              format='commented_header', data_start=0,
                              header_start=-1)

    salt2fitnames = ['sn' + cid.lower() for cid in salt2fitdata['CID']]
    gotsalt2fit = [sn for sn in snnamelist if sn in salt2fitnames]
    nosalt2fit = [sn for sn in snnamelist if sn not in salt2fitnames]
    return gotsalt2fit, nosalt2fit


def get_mlcs_to_salt2_parameter_conversion_functions(
        fitresfilename='lowz_salt2.fitres', showfits=False, verbose=False):
    """ NOTE: this is a really crude kludge of a solution.

    Get the SALT2 x1,c and MLCS delta, Av values for all SNe for which we have
    both.  Fit a simple quadratic to each pair of corresponding parameters.
    :returns delta2x1, av2c: functions that convert from the MLCS parameter
        delta or Av to the SALT2 parameter x1 or c, respectively.
    """
    # read in the low-z SN metadata from the file provided by Arturo
    metadata = load_metadata()
    # read in the low-z SN salt2 fit parameters from the file provided by Dan
    salt2fitfile = os.path.join(__MODELDIR__,'salt2-ir','salt2irdata', fitresfilename)
    salt2fitdata = ascii.read(salt2fitfile,
                              format='commented_header', data_start=0,
                              header_start=-1)
    x1list, x1errlist, clist, cerrlist = [],[],[], []
    deltalist, deltaerrlist, avlist, averrlist, namelist = [],[],[],[], []

    for snname in metadata['snname']:
        imeta = np.where(metadata['snname']==snname)[0]
        snname_stripped = snname.lstrip('sn')
        delta = float(metadata['Delta'][imeta])
        deltaerr = float(metadata['dDelta'][imeta])
        av = float(metadata['Av_mlcs'][imeta])
        averr = float(metadata['dAv'][imeta])
        if delta <= -0.5:
            continue
        if delta > 1:
            continue

        if av > 1.8:
            continue

        if snname_stripped not in salt2fitdata['CID']:
            if verbose:
                print("missing %s in salt2 fit data. Skipping" %
                      snname_stripped)
            continue

        isalt2 = np.where(salt2fitdata['CID']==snname_stripped)[0]
        x1 = np.median(salt2fitdata['x1'][isalt2])
        x1err = np.median(salt2fitdata['x1ERR'][isalt2])
        c = np.median(salt2fitdata['c'][isalt2])
        cerr = np.median(salt2fitdata['cERR'][isalt2])
        if x1<-3:
            continue

        x1list.append(x1)
        x1errlist.append(c)
        clist.append(c)
        cerrlist.append(cerr)
        deltalist.append(delta)
        deltaerrlist.append(deltaerr)
        avlist.append(av)
        averrlist.append(averr)
        namelist.append(snname)

    x1 = np.array(x1list)
    x1err = np.array(x1errlist)
    c = np.array(clist)
    cerr = np.array(cerrlist)
    delta = np.array(deltalist)
    deltaerr = np.array(deltaerrlist)
    av = np.array(avlist)
    averr = np.array(averrlist)

    # TODO : switch to using scipy.odr for fitting with errors in both
    #  dimensions.
    x1fit = scopt.curve_fit(quadratic, delta, x1,
                            p0=None, sigma=x1err,
                            absolute_sigma=True,
                            check_finite=True, )
    x1fitparam = x1fit[0]
    x1fitcov = np.sqrt(np.diag(x1fit[1]))

    cfit = scopt.curve_fit(linear, av, c,
                           p0=None, sigma=cerr,
                           absolute_sigma=True,
                           check_finite=True, )
    cfitparam = cfit[0]
    cfitcov = np.sqrt(np.diag(cfit[1]))

    if showfits:
        fig = pl.gcf()
        fig.clf()
        ax1 = fig.add_subplot(2,1,1)
        pl.errorbar(deltalist, x1list, x1errlist, deltaerrlist, marker='o',
                    color='k', ls=' ')
        ax = pl.gca()
        ax.set_xlabel('MLCS $\Delta$')
        ax.set_ylabel('SALT2 x$_1$')

        ax2 = fig.add_subplot(2,1,2)
        pl.errorbar(avlist, clist, cerrlist, averrlist, marker='d', color='g',
                    ls=' ')
        ax = pl.gca()
        ax.set_xlabel('MLCS $A_V$')
        ax.set_ylabel('SALT2 c')

        deltarange = np.arange(-0.4, 1.0, 0.01)
        ax1.plot( deltarange,
                  quadratic(deltarange, x1fitparam[0],
                            x1fitparam[1], x1fitparam[2]),
                  ls='-', color='r', marker=' ')


        avrange = np.arange(-0.1, 1.9, 0.01)
        ax2.plot( avrange,
                  linear(avrange, cfitparam[0],  cfitparam[1]),
                  ls='-', color='r', marker=' ')

        ax2.set_xlim(-0.1, 1.9)
        pl.draw()

    def mlcsdelta_to_salt2x1(delta):
        return quadratic(delta, x1fitparam[0], x1fitparam[1], x1fitparam[2])

    def mlcsav_to_salt2c(c):
        return linear(c, cfitparam[0], cfitparam[1])

    return mlcsdelta_to_salt2x1, mlcsav_to_salt2c


def linear(x, A, B):
    return A + B * x

def quadratic(x, A, B, C):
    return A + B * x + C * x * x

def cubic(x, A, B, C, D):
    return A + B * x + C * x * x + D * x * x * x


def ccm_unred(wave, flux, ebv, r_v=3.1):
    """ccm_unred(wave, flux, ebv, r_v="")
    Deredden a flux vector using the CCM 1989 (and O'Donnell 1994)
    parameterization. Returns an array of the unreddened flux
    """
    wave = np.array(wave, float)
    flux = np.array(flux, float)
    if wave.size != flux.size:
        raise TypeError('ERROR - wave and flux vectors ' +
                        'must be the same size')

    a_lambda = ccm_extinction(wave, ebv, r_v)

    funred = flux * 10.0**(0.4*a_lambda)

    return funred


def ccm_extinction(wave, ebv, r_v=3.1):
    """
    The extinction (A_lambda) for given wavelength (or vector of wavelengths)
    from the CCM 1989 (and O'Donnell 1994) parameterization. Returns an
    array of extinction values for each wavelength in 'wave'

    INPUTS:
    wave - array of wavelengths (in Angstroms)
    ebv - colour excess E(B-V) float. If a negative ebv is supplied
          fluxes will be reddened rather than dereddened

    OPTIONAL INPUT:
    r_v - float specifying the ratio of total selective
          extinction R(V) = A(V)/E(B-V). If not specified,
          then r_v = 3.1

    OUTPUTS:
    funred - unreddened calibrated flux array, same number of
             elements as wave

    NOTES:
    1. This function was converted from the IDL Astrolib procedure
       last updated in April 1998. All notes from that function
       (provided below) are relevant to this function

    2. (From IDL:) The CCM curve shows good agreement with the Savage & Mathis (1979)
       ultraviolet curve shortward of 1400 A, but is probably
       preferable between 1200 and 1400 A.

    3. (From IDL:) Many sightlines with peculiar ultraviolet interstellar extinction
       can be represented with a CCM curve, if the proper value of
       R(V) is supplied.

    4. (From IDL:) Curve is extrapolated between 912 and 1000 A as suggested by
       Longo et al. (1989, ApJ, 339,474)

    5. (From IDL:) Use the 4 parameter calling sequence if you wish to save the
       original flux vector.

    6. (From IDL:) Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
       curve (3.3 -- 8.0 um-1).    But since their revised curve does
       not connect smoothly with longer and shorter wavelengths, it is
       not included here.

    7. For the optical/NIR transformation, the coefficients from
       O'Donnell (1994) are used

    >>> ccm_unred([1000, 2000, 3000], [1, 1, 1], 2 )
    array([9.7976e+012, 1.12064e+07, 32287.1])
    """
    import numpy as np
    scalar = not np.iterable(wave)
    if scalar:
        wave = np.array([wave], float)
    else:
        wave = np.array(wave, float)

    x = 10000.0/wave
    npts = wave.size
    a = np.zeros(npts, float)
    b = np.zeros(npts, float)

    #Infrared
    good = np.where( (x > 0.3) & (x < 1.1) )
    a[good] = 0.574 * x[good]**(1.61)
    b[good] = -0.527 * x[good]**(1.61)

    # Optical & Near IR
    good = np.where( (x  >= 1.1) & (x < 3.3) )
    y = x[good] - 1.82

    c1 = np.array([ 1.0 , 0.104,   -0.609,    0.701,  1.137,
                  -1.718,   -0.827,    1.647, -0.505 ])
    c2 = np.array([ 0.0,  1.952,    2.908,   -3.989, -7.985,
                  11.102,    5.491,  -10.805,  3.347 ] )

    a[good] = np.polyval(c1[::-1], y)
    b[good] = np.polyval(c2[::-1], y)

    # Mid-UV
    good = np.where( (x >= 3.3) & (x < 8) )
    y = x[good]
    F_a = np.zeros(np.size(good),float)
    F_b = np.zeros(np.size(good),float)
    good1 = np.where( y > 5.9 )

    if np.size(good1) > 0:
        y1 = y[good1] - 5.9
        F_a[ good1] = -0.04473 * y1**2 - 0.009779 * y1**3
        F_b[ good1] =   0.2130 * y1**2  +  0.1207 * y1**3

    a[good] =  1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a
    b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b

    # Far-UV
    good = np.where( (x >= 8) & (x <= 11) )
    y = x[good] - 8.0
    c1 = [ -1.073, -0.628,  0.137, -0.070 ]
    c2 = [ 13.670,  4.257, -0.420,  0.374 ]
    a[good] = np.polyval(c1[::-1], y)
    b[good] = np.polyval(c2[::-1], y)

    # Defining the Extinction at each wavelength
    a_v = r_v * ebv
    a_lambda = a_v * (a + b/r_v)
    if scalar:
        a_lambda = a_lambda[0]
    return a_lambda



def load_metadata(metadatafilename='lowz_metadata.txt'):
    """read in the low-z SN metadata from the file provided by Arturo"""
    metadatafile = os.path.join(__MODELDIR__,'salt2-ir','salt2irdata', metadatafilename)
    metadata = ascii.read(metadatafile)
    return metadata


def deredden_template_sed(sedfile, sedfileout=None, snname=None,
                          metadatafilename='lowz_metadata.txt',
                          fitresfilename='lowz_salt2.fitres'):
    """
    For the given low-z SN, modify the SED template by correcting for
    the redshift and dust extinction (both in the SN host galaxy and in the
    Milky Way along the line of sight).   The result is an un-reddened
    SED template in the rest-frame.

    read in the metadata file
    for each low-z SN data file
    get the metadata
    correct for host galaxy E(B-V)   :  use the MLCS Av and Rv
    correct for redshift  :  use z_helio  (though zcmb is ~ equivalent)
    correct for MW A_V  : these are the E(B-V) and AV_gal columns

    :param sedfile: The SN SED file in the template library.
        If snname is None, then we assume that the root of the file name
        corresponds to the SN ID in the metadata file.
        e.g.  'sn2010ju_cfa'
    :param snname: the SN ID in the metadata file.
    :param sedfileout: filename in which to dump the output arrays
    :return phase, wave, flux: These arrays provide the data for the
        de-reddened template SED.  The array phase is a 1D array with each
        entry giving the rest-frame phase relative to B band max. wave and
        flux are 2D arrays, each with a 1D array for every day in <phase>.
    """
    if snname is None:
        snname = os.path.basename(sedfile).split('_')[0]

    if not os.path.isfile(sedfile):
        raise RuntimeError("No such file %s"% sedfile)

    metadata = load_metadata(metadatafilename)

    if snname not in metadata['snname']:
        #raise RuntimeError(
        print('No SN named %s in the metadata file' % snname)
        return -1

    imeta = np.where(metadata['snname']==snname)[0]
    zhelio = metadata['z_helio'][imeta]
    zcmb = metadata['z_cmb'][imeta]
    EBVmw = metadata['E(B-V)'][imeta]
    Rvmw = 3.1
    Rvhost = metadata['Rv_mlcs'][imeta]
    Avhost = metadata['Av_mlcs'][imeta]
    Delta_mlcs = metadata['Delta'][imeta]
    TBmax = metadata['TBmax'][imeta]
    if Delta_mlcs==-999:
        Delta_mlcs=0
    if Avhost==-999:
        Avhost= 0
    if Rvhost==-999:
        Rvhost= 3.1
    EBVhost = Avhost / Rvhost

    #TODO: convert from Delta_mlcs to SALT2 x1

    # read in the low-z SN salt2 fit parameters from the file provided by Dan
    salt2fitfile = os.path.join(__MODELDIR__,'salt2-ir','salt2irdata', fitresfilename)
    salt2fitdata = ascii.read(salt2fitfile,
                              format='commented_header', data_start=0,
                              header_start=-1)
    snname_stripped = snname.lstrip('sn')
    if snname_stripped in salt2fitdata['CID']:
        isalt2 = np.where(salt2fitdata['CID']==snname_stripped)[0]
        x0 = np.median(salt2fitdata['x0'][isalt2])
        x1 = np.median(salt2fitdata['x1'][isalt2])
        c = np.median(salt2fitdata['c'][isalt2])
        zHD = np.median(salt2fitdata['zHD'][isalt2])
        zsalt2 = np.median(salt2fitdata['z'][isalt2])
        if np.abs(zhelio - zcmb) > 0.01:
            print("Note significant difference in redshift for %s" % snname +
                  " \\n zhelio = %.5f    zsalt2= %.5f" % (zhelio, zsalt2))
    else:
        print("No salt2 fit results available for %s"%snname)
        return -1
        #delta2x1, av2c = get_mlcs_to_salt2_parameter_conversion_functions(
        #    fitresfilename=fitresfilename)
        #x1 = delta2x1(Delta_mlcs)
        #c = av2c(Avhost)
        #zHD = zcmb
        #z = zhelio

    # read in the SED template file directly
    try:
        lowzsn = TimeSeriesGrid(sedfile)
    except AssertionError:
        print("Bad input SED data for %s"%sedfile)
        # return -1

    snphase, snwave, snflux = lowzsn.phase, lowzsn.wave, lowzsn.value
    #snphase=sorted(snphase.tolist())
    # Define new arrays to hold the de-reddened template data
    snphaseout, snwaveout, snfluxout = [], [], []

    if sedfileout is not None:
        fout = open(sedfileout, 'w')
        # print a header that carries all the relevant metadata
        print >> fout, '# name = %s' % snname_stripped
        print >> fout, '# survey = %s' % snname.split('_')[-1]
        print >> fout, '# z = %.5f' % zhelio
        print >> fout, '# salt2x0 = %.5f' % x0
        print >> fout, '# salt2x1 = %.5f' % x1
        print >> fout, '# salt2c = %.5f' % c
        print >> fout, '# mlcsDelta = %.5f' % Delta_mlcs
        print >> fout, '# E(B-V)mw = %.5f' % EBVmw
        print >> fout, '# RVmw = %.5f' % Rvmw
        print >> fout, '# AVhost = %.5f' % Avhost
        print >> fout, '# RVhost = %.5f' % Rvhost
        print >> fout, '# E(B-V)host = %.5f' % EBVhost
        print >> fout, '# TBmax  = %.1f' % TBmax
        print >> fout, '# phase                   wavelength               flux'

    for phase in sorted(snphase.tolist()) :
        iphase = np.where(snphase == phase)[0][0]
        # correct the host extinction, redshift and MW extinction
        phase = snphase[iphase]
        snflux1 = ccm_unred(snwave[iphase], snflux[iphase], EBVhost, Rvhost)
        snwave1 = snwave[iphase] / (1+zhelio)
        snflux2 = ccm_unred(snwave1, snflux1, EBVmw, Rvmw)
        snphaseout.append(phase)
        snwaveout.append(snwave1)
        snfluxout.append(snflux2)
        if sedfileout is not None:
            for w,f in zip(snwave1, snflux2):
                print >> fout, '%25.18e %25.18e %25.18e' % (phase, w, f)
    if sedfileout is not None:
        fout.close()

    return 1


def plot_dereddened_template_comparison(sedfile0, sedfile, phase=0):
    lowzsn0 = TimeSeriesGrid(sedfile0)
    lowzsn1 = TimeSeriesGrid(sedfile)

    snphase0, snwave0, snflux0 = lowzsn0.phase, lowzsn0.wave, lowzsn0.value
    snphase1, snwave1, snflux1 = lowzsn1.phase, lowzsn1.wave, lowzsn1.value
    salt2mod = sncosmo.Model('salt2')

    iphase = np.argmin(np.abs(snphase0-phase))
    bestphase0 = snphase0[iphase]

    minwave = 3500
    salt2wave = np.arange(salt2mod.minwave(),salt2mod.maxwave(),10)
    if minwave <= np.min(salt2wave):
        iwavemin = 0
    else:
        iwavemin = np.where(salt2wave <= minwave)[0][-1]

    salt2flux = salt2mod.flux(bestphase0, salt2wave)
    salt2flux_norm = scint.trapz(salt2flux[iwavemin:], salt2wave[iwavemin:])

    pl.plot(salt2wave, salt2flux/salt2flux_norm, 'k--', lw=1,
            label='SALT2-4 template0' )

    for snphase, snwave, snflux, color, label in zip(
            [snphase0,snphase1],[snwave0,snwave1],[snflux0,snflux1],
            ['r','c'], ['original', 'dereddened']):
        iwavemax = np.where(snwave[iphase]>=np.max(salt2wave))[0][0]
        if minwave <= np.min(snwave[iphase]):
            iwavemin = 0
        else:
            iwavemin = np.where(snwave[iphase] <= minwave)[0][-1]

        snflux_norm = scint.trapz(snflux[iphase][iwavemin:iwavemax],
                                  snwave[iphase][iwavemin:iwavemax])
        pl.plot(snwave[iphase], snflux[iphase]/snflux_norm,
                ls='-', color=color, alpha=1, lw=2, label=label )
    pl.legend(loc='upper right')
    ax = pl.gca()
    ax.set_xlabel('Wavelength ($\AA$)')
    ax.set_xlim(2000, 16000)
    ax.set_ylim(-0.00005, 0.00058)

    snid = os.path.basename(sedfile0).split('_')[0].lstrip('sn')
    ax.text(8500, 0.00015, 'SN ' + snid, fontsize='large')
    pl.draw()

def deredden_template_list(sedfilelist=None,
        showfig=True, savefig=True, clobber=False):
    import glob
    if sedfilelist is None:
        sedfilelist = glob.glob('lowzIaObsFrame/sn*_*.dat')
    for sedfile0 in sedfilelist:
        sedfile0basename = os.path.basename(sedfile0)
        sedfile1basename = sedfile0basename.split('_')[0] + '_' + \
                           sedfile0basename.split('__')[-1].split('_')[0] + '.dat'
        snname = sedfile0basename.split('_')[0]
        sedfile1 = os.path.join(__MODELDIR__,'salt2-ir','salt2irdata','lowzIa', sedfile1basename)
        if os.path.exists(sedfile1) and not clobber:
            print("%s exists. Not clobbering" % sedfile1)
            continue
        result = deredden_template_sed(sedfile0, sedfile1, snname=snname)
        if result<0:
            print("!!!  Failed to Deredden %s  !!!" % snname)
            continue
        else:
            print("Successfully Dereddened %s" % snname)

        if showfig or savefig:
            pl.clf()
            plot_dereddened_template_comparison(sedfile0, sedfile1, phase=0)
            if savefig:
                figfilename = 'lowzIa/%s.png' % sedfile1basename.split('.')[0]
                pl.savefig(figfilename)

def load_all_templates():
    templatedict = {}
    sedfilelist = os.listdir('lowzIa')
    for sedfile in sedfilelist:
        snname = sedfile.split('.')[0]
        sedfile = os.path.join('lowzIa', snname) + '.dat'
        templatedict[snname] = TimeSeriesGrid(sedfile=sedfile)
    return templatedict


def extend_template0_ir(modeldict = None, x1min=-1, x1max=1,
                        modeldir=__MODELDIR__,
                        salt2indir = 'salt2-4',
                        salt2outdir = 'salt2-ir',
                        wavejoin = 8500, wavemax = 24990, verbose=False):
    """ extend the salt2 Template_0 model component
    by adopting the IR tails from a collection of SN Ia template SEDs.
    Here we use the collection of CfA, CSP, and other low-z SNe provided by
    Arturo Avelino (2016, priv. comm.)
    The median of the sample is scaled and joined at the
    wavejoin wavelength, and extrapolated out to wavemax.
    """
    if modeldict == None:
        modeldict = load_sncosmo_models()
    salt2indir = os.path.join(modeldir, salt2indir)
    salt2outdir = os.path.join(modeldir, salt2outdir)

    temp0fileIN = os.path.join( salt2indir, 'salt2_template_0.dat' )
    temp0fileOUT = os.path.join( salt2outdir, 'salt2_template_0.dat' )
    temp0 = TimeSeriesGrid(temp0fileIN)
    wavestep = np.median(np.diff(temp0.wave[0]))
    waveir = np.arange(wavejoin, wavemax+wavestep, wavestep)

    outlines = []
    for iphase in range(len(temp0.phase)): # iphaselist:
        # get the SALT2 template0 SED for this day
        # phase0 = temp0.phase[iphase]
        wave0 = temp0.wave[iphase]
        flux0 = temp0.value[iphase]
        thisphase = temp0.phase[iphase]

        ijoin = np.argmin(np.abs(wave0-wavejoin))
        fluxjoin = flux0[ijoin]
        if verbose:
            print( 'splicing tail onto template for day : %i'%thisphase )

        # get the median of all the low-z mangled SEDs that have data at
        # this phase and satisfy the x1 range requirements
        fluxlowzarray = []
        for name in modeldict.keys():
            if not name.startswith('sn'):
                continue
            lowzsn = modeldict[name]
            if lowzsn.maxwave() <wavemax:
                continue

            if lowzsn.mintime()>thisphase: continue
            if lowzsn.maxtime()<thisphase: continue
            if lowzsn.salt2x1<x1min: continue
            if lowzsn.salt2x1>x1max: continue

            fluxlowz = lowzsn.flux(thisphase, waveir)
            if np.sum(fluxlowz)==0: continue

            # normalize the flux of this lowz mangled SED so that it matches
            # the flux of the salt2 template0 model at the join wavelength
            normalization_factor = (fluxjoin / fluxlowz[0])
            fluxlowzarray.append(fluxlowz * normalization_factor)

        # extend the template0 SED into the IR for this phase
        # using the median of all the lowz templates (or the Hsiao model when
        # templates are not available
        if len(fluxlowzarray)<5:
            print("only %i templates for phase = %.1f. Using Hsiao model" % (
                len(fluxlowzarray), thisphase))
            fluxlowzmedian = modeldict['hsiao'].flux(thisphase, waveir)
        else:
            fluxlowzarray = np.array(fluxlowzarray)
            fluxlowzmedian = np.median(fluxlowzarray, axis=0)

        ijoin0 = np.argmin(abs(wave0 - wavejoin))
        ijoinlowz = np.argmin(abs(waveir - wavejoin))

        if flux0[ijoin0]>0:
            scalefactor = fluxlowzmedian[ijoinlowz]/flux0[ijoin0]
        else:
            scalefactor = 1
        wavenew = np.append(wave0[:ijoin0], waveir.tolist())
        fluxnew = np.append(flux0[:ijoin0], (scalefactor*fluxlowzmedian))

        # append to the list of output data lines
        for j in range( len(wavenew) ) :
            outlines.append( "%6.2f    %12i  %12.7e\n"%(
                    thisphase, wavenew[j], fluxnew[j] ) )

    # write it out to the new template0 .dat file
    fout = open( temp0fileOUT, 'w' )
    fout.writelines( outlines )
    fout.close()

def _round_up_to_odd(f):
    f = int(np.ceil(f))
    return f + 1 if f % 2 == 0 else f

def plot_extended_template0(salt2indir = 'salt2ir/salt2-4', salt2outdir = 'data/IRSalt2'):
    """Plot the extended SALT2 template0 component at 5 phases"""
    pl.clf()
    oldtemp0 = TimeSeriesGrid(os.path.join(salt2indir,'salt2_template_0.dat'))
    newtemp0 = TimeSeriesGrid(os.path.join(salt2outdir,'salt2_template_0.dat'))
    for phase, fluxoffset, color in zip([-10,-5,0,5,25],
                                        [0.14, 0.09, 0.06, 0.03, 0.0],
                                        ['m','b','g','r','k']):
        oldtemp0.plot_at_phase(phase=phase, fluxoffset=fluxoffset,
                               color=color, ls='-',
                               label='phase = %i'%int(phase))
        newtemp0.plot_at_phase(phase=phase, fluxoffset=fluxoffset,
                               color=color, ls='--',
                               label='__nolabel__')
    pl.legend(loc='upper right', ncol=2)
    ax = pl.gca()
    ax.set_xlim(5800, 15000)
    ax.set_ylim(-0.01, 0.31)
    pl.setp(ax.get_yticklabels(), visible=False)
    ax.set_ylabel('SALT2 M0 value (arbitrary units)',fontsize=14)
    ax.set_xlabel(r'Wavelength ($\rm{\AA}$)',fontsize=14)
    ax.set_title('SALT2 and SALT2-IR at Various Phases',fontsize=16)


def extend_template1_ir(modeldict = None, cmin=-0.1, cmax=0.3,
                        x1min=-1, x1max=1,
                        modeldir=__MODELDIR__,
                        salt2indir = 'salt2-4', salt2outdir = 'salt2-ir',
                        wavejoinstart = 8500, wavemax = 24990,
                        wavesmoothwindow=100, timesmoothwindow=5,
                        showplots=False, verbose=False):
    """ extend the M1 component of the SALT2 model into the IR.  This is the
    component that is multiplied by the SALT2 parameter x1, and therefore
    reflects changes in the shape of the light curve.

    :param modeldict:
    :param modeldir:
    :param salt2indir:
    :param salt2outdir:
    :param wavejoinstart:
    :param wavemax:
    :param wavesmoothwindow: smoothing window size in Angstrom
    :param showphase0plots:
    :return:
    """
    medsmooth = lambda f, N : np.array([np.median(
        f[max(0,i-N):min(len(f),max(0,i-N)+2*N)]) for i in range(len(f))])

    if showplots:
        phaselist = [-5, 0, 5, 10, 25]

        fig = pl.gcf()
        fig.clf()
        ax1 = pl.subplot2grid((len(phaselist), 2), (0,0), rowspan=len(phaselist))

    if modeldict == None:
        modeldict = load_sncosmo_models()
    salt2indir = os.path.join(modeldir, salt2indir)
    salt2outdir = os.path.join(modeldir, salt2outdir)

    # read in the extended template0 data (i.e. the M0 model
    # component that has already been extended to IR wavelengths)
    temp0extfile = os.path.join( salt2outdir, 'salt2_template_0.dat' )
    temp0ext = TimeSeriesGrid(temp0extfile)

    temp1fileIN = os.path.join( salt2indir, 'salt2_template_1.dat' )
    temp1 = TimeSeriesGrid(temp1fileIN)

    temp1fileOUT = os.path.join( salt2outdir, 'salt2_template_1.dat' )
    outlines = []

    # define wavelength arrays. The same wavelength arrays are used
    # for every phase
    wavelist_old = temp1.wave[0]
    wavestep = wavelist_old[1] - wavelist_old[0]
    wavemin = wavelist_old[0]
    wavejoinend = wavelist_old[-1]
    wavelist_join = np.arange(wavejoinstart, wavejoinend+wavestep, wavestep)
    wavelist_opt = np.arange(4500, 7000, wavestep)
    wavelist_ir = np.arange(wavejoinstart, wavemax+wavestep, wavestep)
    wavelist_all = np.arange(wavemin, wavemax+wavestep, wavestep)

    # indices to pick out the join sections for the Old and new M1 models
    ijoinold = np.where((wavelist_old >= wavejoinstart) &
                        (wavelist_old <= wavejoinend))[0]
    ijoinnew = np.where((wavelist_ir >= wavejoinstart) &
                        (wavelist_ir <= wavejoinend))[0]

    assert np.all(wavelist_all==temp0ext.wave[0])
    #        "New M1 wavelength array must match extended M0"
    #        " wavelength array exactly.")

    assert np.all(wavelist_all == np.append(wavelist_old[:ijoinold[0]],
                                            wavelist_ir))
    #        "New M1 wavelength array must match the join of old (optical)"
    #        " and new (IR) wavelength arrays exactly.")

    # define weight functions to be applied when combining the new
    #  M1 curve with the old M1 curve.  It increases linearly from 0
    # at the start of the join window (8500 A) to 1 at the end of
    # the join window (9200 A), which is the last wavelength for the
    # old optical SALT2 model
    whtslope = 1.0 / (wavejoinend - wavejoinstart)
    whtintercept = -whtslope * wavejoinstart
    newM1weight = whtslope * wavelist_join + whtintercept
    oldM1weight = 1 - newM1weight

    newM1grid = []
    for iphase1 in range(len(temp1.phase)): # iphaselist:
        thisphase = temp1.phase[iphase1]

        # get the SALT2 template0 and template1 data for this day
        # and define interpolating functions
        iphase0 = np.argmin(np.abs(temp0ext.phase-thisphase))
        assert iphase0 == iphase1

        M0func = scinterp.interp1d(
            temp0ext.wave[iphase0], temp0ext.value[iphase0],
            bounds_error=False, fill_value=0)
        M1func = scinterp.interp1d(
            temp1.wave[iphase1], temp1.value[iphase1],
            bounds_error=False, fill_value=0)

        if verbose:
            print( 'solving for IR tail of template_1 for day : %i'%thisphase )

        M1extlist = []
        for name in modeldict.keys():
            if not name.startswith('sn'):
                continue
            lowzsn = modeldict[name]
            if lowzsn.maxwave()<wavemax:
                continue
            if lowzsn.mintime()>thisphase: continue
            if lowzsn.maxtime()<thisphase: continue
            if lowzsn.salt2c<cmin: continue
            if lowzsn.salt2c>cmax: continue
            if lowzsn.salt2x1<x1min: continue
            if lowzsn.salt2x1>x1max: continue

            fluxlowz_opt = lowzsn.flux(thisphase, wavelist_opt)
            if np.sum(fluxlowz_opt)==0: continue

            # load the original SALT2 colorlaw:
            CL = load_salt2_colorlaw(salt2indir=salt2indir)

            # Solve for x0 as a function of lambda (should be ~constant!)
            # over optical wavelengths. Then define a scalar x0 as the
            # median value over the optical wavelength range
            x0_vs_lambda = (fluxlowz_opt /
                            ((M0func(wavelist_opt) +
                             lowzsn.salt2x1 * M1func(wavelist_opt)) *
                             np.exp(lowzsn.salt2c * CL(wavelist_opt))))
            x0median = np.median(x0_vs_lambda)

            if False:
                pl.clf()
                pl.plot(wavelist_opt, x0_vs_lambda)
                ax = pl.gca()
                ax.axhline(x0median, color='r', ls='-')
                ax.text(0.95, 0.95, '%.3f'%lowzsn.salt2x0,
                        ha='right', va='top', color='r', transform=ax.transAxes)
                pl.show()
                pl.draw()
                raw_input("Showing x0 vs lambda for %s. Return to continue"%name)

            x0 = x0median

            # solve for the M1 array over the NIR wavelength range (8500+ A)
            fluxlowz_ir = lowzsn.flux(thisphase, wavelist_ir)
            M1ext = (1/lowzsn.salt2x1) * ((fluxlowz_ir/x0) - M0func(wavelist_ir))
            M1extlist.append(M1ext)

        if len(M1extlist)<3:
            print("only %i templates for phase = %.1f." % (
                len(M1extlist), thisphase))
            print("Using M1(lambda)=0 for all IR wavelengths")
            newM1median = np.zeros(len(wavelist_ir))
        else:
            newM1median = np.median(M1extlist, axis=0)

        if wavesmoothwindow:
            Nsmooth = int(wavesmoothwindow / wavestep)
            newM1median = medsmooth(newM1median, Nsmooth)

        # apply the predefined weight functions to combine the new M1
        # curve with the old M1 curve in the "join window"
        # and then stitch together the old M1 curve up to the join
        # wavelength, the joining curve through the join window, and the
        # new M1 curve beyond that into the IR wavelengths.
        newM1joinvalues = (oldM1weight * temp1.value[iphase1][ijoinold] +
                           newM1weight * newM1median[ijoinnew])
        newM1values = np.append(
            np.append(temp1.value[iphase1][:ijoinold[0]], newM1joinvalues),
            newM1median[ijoinnew[-1]+1:])
        newM1grid.append(newM1values)

        if len(newM1values) != len(wavelist_all):
            import pdb; pdb.set_trace()

        if showplots and thisphase==0:
            for M1ext in M1extlist:
                ax1.plot(wavelist_ir, M1ext, lw=2, alpha=0.3, color='b',
                         label='__nolabel__')
            #ax1.plot(wavelist_ir, newM1median, lw=3.5, color='0.5',
            #         label='median of IR templates')
            ax1.plot(temp1.wave[0], M1func(temp1.wave[0]), color='0.5', lw=3.5,
                     label='SALT2-4')
            ax1.plot(wavelist_all, newM1values, lw=1.2, color='darkorange',
                     label='SALT2-IR')
            ax1.set_ylim(-0.035, 0.045)
            ax1.set_xlim(2805, 18495)
            ax1.set_xlabel(r'Wavelength ($\rm{\AA}$)')
            ax1.set_ylabel('SALT2 M1 model component value')
            pl.legend(loc='upper right')

    newM1grid = np.array(newM1grid)
    if timesmoothwindow:
        for iwave in range(len(wavelist_all)):
            newM1grid[:,iwave] = savitzky_golay(
                newM1grid[:,iwave], timesmoothwindow)

    for iphase1 in range(len(temp1.phase)): # iphaselist:
        thisphase = temp1.phase[iphase1]
        newM1values = newM1grid[iphase1]

        # append to the list of output data lines
        for j in range(len(wavelist_all)) :
            outlines.append( "%6.2f    %12i  %12.7e\n"%(
                    thisphase, wavelist_all[j], newM1values[j] ) )

    # write it out to the new template sed .dat file
    fout = open( temp1fileOUT, 'w' )
    fout.writelines( outlines )
    fout.close()

    if showplots:
        oldtemp1 = TimeSeriesGrid('salt2-4/salt2_template_1.dat')
        newtemp1 = TimeSeriesGrid('salt2ir/salt2_template_1.dat')

        iax = -1
        for phase in phaselist:
            if phase==0: continue
            iax+=1
            ax = pl.subplot2grid((len(phaselist)-1, 2), (iax, 1), sharex=ax1)
            oldtemp1.plot_at_phase(phase=phase, color='0.5', ls='-', lw=2.5,
                                   label='SALT2-4')
            newtemp1.plot_at_phase(phase=phase, color='darkorange', ls='-', lw=1,
                                   label='SALT2-IR')
            ax.text(0.5,0.9, 'phase = % i'%int(phase), ha='left', va='top',
                    fontsize='large', transform=ax.transAxes)

            if phase==phaselist[-1]:
                ax.set_xlabel(r'Wavelength ($\rm{\AA}$)')
            else:
                pl.setp(ax.get_xticklabels(), visible=False)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
            ax.yaxis.set_minor_locator(ticker.MaxNLocator(6))
            ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
            ax.xaxis.set_minor_locator(ticker.MaxNLocator(10))
            ax.set_xlim(2805, 17995)
            
            pl.setp(ax.get_yticklabels(), fontsize=7)
        ax1.xaxis.set_major_locator(ticker.MaxNLocator(5))
        ax1.xaxis.set_minor_locator(ticker.MaxNLocator(10))
        pl.suptitle('SALT2 and SALT2-IR M1 Model component',fontsize=16)
        fig = pl.gcf()
        fig.subplots_adjust(hspace=0, left=0.1, bottom=0.12, right=0.95)
        pl.draw()
    return



def extend_dispersion_variance_covariance( wref=8500,
        modeldir=__MODELDIR__,
        salt2indir='salt2-4', salt2irsubdir='salt2-ir',
        showplots = False):
    """ extrapolate the *lc* and *spec* .dat files for SALT2
     using a flatline extension of the red tail to 2.5 microns
     """
    salt2modeldir = os.path.join(modeldir, salt2indir)
    salt2irmodeldir = os.path.join(modeldir, salt2irsubdir)

    filelist = ['salt2_lc_dispersion_scaling.dat',
                'salt2_lc_relative_covariance_01.dat',
                'salt2_lc_relative_variance_0.dat',
                'salt2_lc_relative_variance_1.dat',
                'salt2_spec_covariance_01.dat',
                'salt2_spec_variance_0.dat',
                'salt2_spec_variance_1.dat']
    titlelist = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

    if showplots:
        fig = pl.gcf()
        fig.clf()
        iax=0

    for filename, title in zip(filelist, titlelist):
        infile = os.path.join(salt2modeldir, filename)
        outfile = os.path.join(salt2irmodeldir, filename)
        timeseries = TimeSeriesGrid(infile)

        if 'dispersion' in filename:
            timeseries.extrapolate_linearto0(
                outfilename=outfile, refwave=8500,
                maxwave=25000, maxwaveval=1)
        elif 'lc' in filename:
            timeseries.extrapolate_linearto0(
                outfilename=outfile, refwave=wref,
                maxwave=25000, maxwaveval=0)
        else:
            timeseries.extrapolate_flatline(
                outfilename=outfile, refwave=wref,
                maxwave=25000)
        if showplots:
            iax+=1
            if iax==1:
                ax1 = pl.subplot2grid((4,2), (iax-1,0))
                ax = ax1
            elif iax<=4:
                ax = pl.subplot2grid((4,2), (iax-1,0))
            else:
                ax = pl.subplot2grid((4,2), (iax-5,1))
            timeseries.plot_at_phase(0, color='0.5', lw=3.2, ls='-', marker=' ')
            timeseriesNew = TimeSeriesGrid(outfile)
            timeseriesNew.plot_at_phase(0, color='darkorange', lw=1,
                                        ls='--', marker=' ')
            ax.text(0.90, 0.9, title, transform=ax.transAxes,
                    fontsize='large', ha='right', va='top')
            if iax==1:
                ax.set_ylim(-0.2, 6)
            elif iax in [2,3,4]:
                ax.set_ylim(-5e-4,3e-3)
            else:
                ax.set_ylim(-5e-6,5e-5)
            if iax!=4 and iax!=7:
                pl.setp(ax.get_xticklabels(), visible=False)
            pl.setp(ax.get_yticklabels(), visible=False)
            ax.set_xlim(2000,25000)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(10000))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(2000))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
            ax.yaxis.set_minor_locator(ticker.MaxNLocator(6))
            pl.subplots_adjust(hspace=0, wspace=0)

    # read in the salt2_color_correction.dat file
    # and repeat the fit parameters into the .INFO file:
    colorlawfile = os.path.join(salt2irmodeldir,
                                'salt2_color_correction.dat')
    fin = open(colorlawfile, 'r')
    words = fin.read().split()
    fin.close()

    # Get colorlaw coeffecients.
    npoly = int(words[0])
    colorlaw_coeffs = [float(word) for word in words[1: 1 + npoly]]

    # Look for keywords in the rest of the file.
    version = 0
    colorlaw_range = [3000., 7000.]
    for i in range(1+npoly, len(words)):
        if words[i] == 'Salt2ExtinctionLaw.version':
            version = int(words[i+1])
        if words[i] == 'Salt2ExtinctionLaw.min_lambda':
            colorlaw_range[0] = float(words[i+1])
        if words[i] == 'Salt2ExtinctionLaw.max_lambda':
            colorlaw_range[1] = float(words[i+1])
    COLORLAWLINE = 'COLORCOR_PARAMS: %i %i %i ' % (
        colorlaw_range[0], colorlaw_range[1], npoly)
    for i in range(npoly):
        COLORLAWLINE += ' %.8f' % colorlaw_coeffs[i]

    outinfo = os.path.join(salt2irmodeldir, 'SALT2.INFO')
    fout = open(outinfo, 'w')
    print >> fout, """
# open rest-lambda range WAAAY beyond nominal 2900-7000 A range.
RESTLAMBDA_RANGE:  2000. 25000.
COLORLAW_VERSION: 1
""" +\
COLORLAWLINE + \
"""
COLOR_OFFSET:  0.0

MAG_OFFSET: 0.27
SEDFLUX_INTERP_OPT: 1  # 1=>linear,    2=>spline
ERRMAP_INTERP_OPT:  1  # 0=snake off;  1=>linear  2=>spline
ERRMAP_KCOR_OPT:    1  # 1/0 => on/off

MAGERR_FLOOR:   0.005            # don;t allow smaller error than this
MAGERR_LAMOBS:  0.1  2000  4000  # magerr minlam maxlam
MAGERR_LAMREST: 0.1   100   200  # magerr minlam maxlam
"""

def load_salt2_colorlaw(salt2indir='salt2-4'):
    """Read in the polynomial parameters for a salt2 color law and
    return a callable function that provides the colorlaw value CL(l) 
    for any input wavelength or wavelength array. """
    clfilename = os.path.join(salt2indir,'salt2_color_correction.dat')
    fin = open(clfilename, 'r')

    nparam = int(fin.readline().strip())
    params = []
    for iline in range(nparam):
        params.append(float(fin.readline().strip()))
    version = int(fin.readline().split()[1])
    rangemin = float(fin.readline().split()[1])
    rangemax = float(fin.readline().split()[1])
    cl = lambda w: salt2_colorlaw(w, params,
                                  colorlaw_range=[rangemin, rangemax])
    return cl


def salt2_colorlaw(wave, params,
                   colorlaw_range=[2800,7000]):
    """Return the  extinction in magnitudes as a function of wavelength,
    for c=1. This is the version 1 extinction law used in SALT2 2.0
    (SALT2-2-0).  From S. Bailey's "salt2" github. 

    Notes
    -----
    From SALT2 code comments:

    if(l_B<=l<=l_R):
        ext = exp(color * constant *
                  (alpha*l + params(0)*l^2 + params(1)*l^3 + ... ))
            = exp(color * constant * P(l))

        where alpha = 1 - params(0) - params(1) - ...

    if (l > l_R):
        ext = exp(color * constant * (P(l_R) + P'(l_R) * (l-l_R)))
    if (l < l_B):
        ext = exp(color * constant * (P(l_B) + P'(l_B) * (l-l_B)))
    """
    v_minus_b = _V_WAVELENGTH - _B_WAVELENGTH
    l = (wave - _B_WAVELENGTH) / v_minus_b
    l_lo = (colorlaw_range[0] - _B_WAVELENGTH) / v_minus_b
    l_hi = (colorlaw_range[1] - _B_WAVELENGTH) / v_minus_b

    alpha = 1. - sum(params)
    coeffs = [0., alpha]
    coeffs.extend(params)
    coeffs = np.array(coeffs)
    prime_coeffs = (np.arange(len(coeffs)) * coeffs)[1:]

    extinction = np.empty_like(wave)

    # Blue side
    idx_lo = l < l_lo
    p_lo = np.polyval(np.flipud(coeffs), l_lo)
    pprime_lo = np.polyval(np.flipud(prime_coeffs), l_lo)
    extinction[idx_lo] = p_lo + pprime_lo * (l[idx_lo] - l_lo)

    # Red side
    idx_hi = l > l_hi
    p_hi = np.polyval(np.flipud(coeffs), l_hi)
    pprime_hi = np.polyval(np.flipud(prime_coeffs), l_hi)
    extinction[idx_hi] = p_hi + pprime_hi * (l[idx_hi] - l_hi)

    # In between
    idx_between = np.invert(idx_lo | idx_hi)
    extinction[idx_between] = np.polyval(np.flipud(coeffs), l[idx_between])

    return -extinction


def fit_salt2colorlaw_to_ccm(c=0.1, Rv=3.1,
                             wmin=2000, wjoin=6750, wmax=25000,
                             salt2_colorlaw_range = [2800, 7000],
                             modeldir=__MODELDIR__,
                             salt2indir='salt2-4', salt2irsubdir='salt2-ir',
                             uselowzfit=False, showfit=False):
    """ Find the SALT2 color law parameters that will cause the color law to
    approximately match the Cardelli+ 1989 extinction law (as extended
    by O'Donnell 1994) for IR wavelengths.
    :param c: SALT2 color parameter (~E(B-V))
    :return:
    """
    if uselowzfit:
        # Use the function that converts from MLCS A_V to SALT2 c
        # to find the value of E(B-V) for this value of c
        d2x1, av2c = get_mlcs_to_salt2_parameter_conversion_functions()
        avrange = np.arange(0,5,0.01)
        crange = av2c(avrange)
        ithisc = np.argmin(np.abs(crange-c))
        Av = avrange[ithisc]
        EBV = avrange/Rv
    else:
        # To first order (and for Rv=3.1) we can use: E(B-V) = c
        EBV = c

    # Load some functions for computing the CCM extinction law
    alambda = lambda w: ccm_extinction(w, EBV, Rv)
    aB = ccm_extinction(_B_WAVELENGTH, EBV, Rv)

    # define a wavelength grid that extends from the blue end of the SALT2
    # color law range out to the extreme red end that we are extrapolating to
    wave = np.append(np.arange(salt2_colorlaw_range[0], salt2_colorlaw_range[1], 1.0),
                     np.arange(salt2_colorlaw_range[1], wmax, 100.0))

    # The color-law parameters for SALT2-4:
    param0 = [-0.504294, 0.787691, -0.461715, 0.0815619]

    # First-guess parameters for the new fit:
    paramA = np.array([-0.5, 0.8, -0.5, 0.1])

    def chi2_function(params, w):
        salt2ext = c * salt2_colorlaw(
            w, params, colorlaw_range=salt2_colorlaw_range)
        ccmext = np.where(
            w < wjoin, c * salt2_colorlaw(w, param0), alambda(w) - aB)
        chi2 = (salt2ext-ccmext)**2
        return np.sum(chi2)

    fitres = scopt.minimize( chi2_function, paramA, args=(wave,),
                             method='Nelder-Mead')
    paramfit = fitres['x']

    if showfit:
        pl.clf()
        salt2colorlaw0 = salt2_colorlaw(wave, params=param0[:4],
                                        colorlaw_range=[2800,7000])
        salt2colorlawfit = salt2_colorlaw(wave, params=paramfit,
                                          colorlaw_range=salt2_colorlaw_range)

        pl.plot(wave, c * salt2colorlaw0, 'k-', label='SALT2-4 color law')
        pl.plot(wave, c * salt2colorlawfit, 'r--', label='SALT2ir color law')
        pl.plot(wave, alambda(wave) - aB, 'b-.', label='Dust with Rv=3.1')
        ax = pl.gca()
        ax.set_xlim(3000,20000)
        ax.set_ylim(-0.45,0.25)

        ax.text(0.95, 0.6,
                'SALT2-4 colorlaw range = [%i, %i]'%tuple([2800,7000]),
                ha='right', va='bottom', transform=ax.transAxes, color='k')
        ax.text(0.95, 0.5,
                ('SALT2-4 colorlaw parameters =\n'
                 '[%.3f, %.3f, %.3f, %.3f]'%tuple(param0)),
                 ha='right', va='bottom', transform=ax.transAxes, color='k')

        ax.text(0.95, 0.35,
                'SALT2ir colorlaw range = [%i, %i]'%tuple(salt2_colorlaw_range),
                ha='right', va='bottom', transform=ax.transAxes, color='r')
        ax.text(0.95, 0.25,
                ('SALT2ir colorlaw parameters =\n'
                '[%.3f, %.3f, %.3f, %.3f]'%tuple(paramfit)),
                ha='right', va='bottom', transform=ax.transAxes, color='r')
        #ax.text(0.95, 0.2,
        #        '$\lambda_{\\rm join}$= %i'%wjoin,
        #        ha='right', va='bottom', transform=ax.transAxes, color='r')
        ax.legend(loc='upper right')
        ax.set_xlabel('Wavelength ($\AA$)',fontsize=14)
        ax.set_ylabel('A$_{\lambda}$ - A$_{B}$  or  c$\\times$ CL($\lambda$)',fontsize=14)
        ax.set_title('SALT2-IR Color Law Extrapolation',fontsize=16)
        #pl.draw()
        pl.savefig('figs/color.pdf',format='pdf',overwrite=True)

    # write out the revised color law as a file
    outfile= os.path.join(modeldir,
                          salt2irsubdir + '/salt2_color_correction.dat')
    fout = open(outfile, 'w')
    print >> fout, '%i' % len(paramfit)
    for param in paramfit:
        print >> fout, '%.8f' % param
    print >> fout, 'Salt2ExtinctionLaw.version 1'
    print >> fout, 'Salt2ExtinctionLaw.min_lambda %i' % salt2_colorlaw_range[0]
    print >> fout, 'Salt2ExtinctionLaw.max_lambda %i' % salt2_colorlaw_range[1]
    fout.close()
    print("Updated SALT2 color law parameters written to %s" % outfile)
    return


def mk_salt2ir_template0_plots(modeldir=__MODELDIR__):
    salt2outdir='salt2-ir'
    salt2irmodeldir = os.path.join(modeldir, salt2outdir)
    salt2irsource = sncosmo.models.SALT2Source(modeldir=salt2irmodeldir, name='salt2-ir')
    salt2irmodel = sncosmo.Model(source=salt2irsource)
    wir = np.arange(2000, 24000, 10)

    salt2model = sncosmo.Model('salt2')
    wopt = np.arange(2000, 9000, 10)

    iax = 0
    fig = pl.gcf()
    fig.clf()
    for phase  in [-5,0,15,30]:
        iax+=1
        fig.add_subplot(2,2,iax)
        fir = salt2irmodel.flux(phase, wir)
        fopt = salt2model.flux(phase, wopt)

        pl.plot( wir, fir, 'r-', label='salt2-ir')
        pl.plot( wopt, fopt, 'k-', label='SALT2-4')
        ax = pl.gca()
        ax.text(0.5, 0.8, 'Phase=%i' % phase,
                transform=ax.transAxes, fontsize='large')
        ax.set_xlim(2000,16000)

        if iax==2:
            ax.legend(loc='center right')

    ax.set_xlabel('wavelength (Angstrom)')
    pl.draw()


def plot_salt2_component_timeseries(
        component=0, wavemin=10000, wavemax=20000,
        modeldir=__MODELDIR__,
        salt2indir = 'salt2-4', salt2outdir = 'salt2-ir',
        **kwargs):
    salt2indir = os.path.join(modeldir, salt2indir)
    salt2outdir = os.path.join(modeldir, salt2outdir)

    temp0fileOUT = os.path.join( salt2outdir, 'salt2_template_0.dat' )
    temp1fileOUT = os.path.join( salt2outdir, 'salt2_template_1.dat' )

    if component == 1 :
        valuegrid = TimeSeriesGrid(temp1fileOUT)
    else:
        valuegrid = TimeSeriesGrid(temp0fileOUT)

    wavelist = np.linspace(wavemin, wavemax, 5)
    for wave in wavelist:
        phaselist = np.arange(-5,25,1)
        iwave = np.argmin(np.abs(valuegrid.wave[0]-wave))
        iphaselist = np.array([np.where(valuegrid.phase == phase)[0][0]
                               for phase in phaselist])
        valuelist = np.array([valuegrid.value[ip][iwave] for ip in iphaselist])
        pl.plot(phaselist, valuelist, label='%i Ang' % wave, **kwargs)
    ax = pl.gca()
    ax.set_xlabel('time (days)')
    ax.set_ylabel('SALT2 model component M%i value' % component)
    ax.legend(loc='best', ncol=2)



def savitzky_golay(y, window_size=5, order=3, deriv=0):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as pl
    pl.plot(t, y, label='Noisy signal')
    pl.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    pl.plot(t, ysg, 'r', label='Filtered signal')
    pl.legend()
    pl.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m, y, mode='valid')


def extend_dispersion_files(refwaveUV=3200, minwave=1700,
                            refwaveIR=8500, maxwave=25000,
                            modeldir=__MODELDIR__,
                            salt2indir='salt2-4', salt2irsubdir='salt2-extended',
                            extrapolate=True, showplots=True):
    """ extend and/or plot the *dispersion.dat files for SALT2 that define
    the model uncertainty
     -- using a downward linear extension of the IR tail
     -- TODO: using a flatline extension of the UV tail
     """
    #TODO: add a flatline extension on the UV side

    salt2modeldir = os.path.join(modeldir, salt2indir)
    salt2irmodeldir = os.path.join(modeldir, salt2irsubdir)

    filelist = ['salt2_lc_dispersion_scaling.dat',
                'salt2_color_dispersion.dat']
    titlelist = ['Light Curve Dispersion Scaling',
                 'Color Dispersion Scaling']

    if showplots:
        fig = pl.gcf()
        fig.clf()
        iax=0

    for filename, title in zip(filelist, titlelist):
        infile = os.path.join(salt2modeldir, filename)
        outfile = os.path.join(salt2irmodeldir, filename)

        if extrapolate:
            if 'color' in filename:
                datagrid = WavelengthGrid(infile)

                datagrid.extrapolate_linearto0(
                    outfilename=outfile, refwave=refwaveIR,
                    maxwave=maxwave, maxwaveval=0.01)
                datagrid = WavelengthGrid(outfile)

            else:
                datagrid = TimeSeriesGrid(infile)

                # linear extrapolation to zero at 25000 angstroms
                datagrid.extrapolate_linearto0(
                    outfilename=outfile, refwave=refwaveIR,
                    maxwave=maxwave, maxwaveval=1.0)
                datagrid = TimeSeriesGrid(outfile)

        if showplots:
            iax+=1
            if iax==1:
                ax1 = pl.subplot(211)
                ax = ax1
            else:
                ax = pl.subplot(212)

            if 'color' in filename:
                datagrid = WavelengthGrid(infile)
                datagrid.plot(color='0.5', lw=3.2, ls='-', marker=' ')
                datagridNew = WavelengthGrid(outfile)
                datagridNew.plot(color='darkorange', lw=1,
                                            ls='--', marker=' ')
                ax.set_yscale('log')
                ax.set_ylim(0.01, 0.98)
            else:
                datagrid = TimeSeriesGrid(infile)
                for phase, color in zip([-5,0,5,15],
                                        ['r','darkorange','g','b']):
                    datagrid.plot_at_phase(phase, color='0.5', lw=3.2,
                                           ls='-', marker=' ',
                                           label='__nolabel__')
                    datagridNew = TimeSeriesGrid(outfile)
                    datagridNew.plot_at_phase(phase, color=color, lw=1,
                                              ls='--', marker=' ',
                                              label='t=%i'%phase)
                pl.legend(loc='upper right')

            ax.text(0.90, 0.9, title, transform=ax.transAxes,
                    fontsize='large', ha='right', va='top')
            ##if iax==1:
            #    ax.set_ylim(-0.2, 6)
            #elif iax in [2,3,4]:
            #    ax.set_ylim(-5e-4,3e-3)
            #else:
            #    ax.set_ylim(-5e-6,5e-5)
            #if iax!=4 and iax!=7:
            #    pl.setp(ax.get_xticklabels(), visible=False)
            #pl.setp(ax.get_yticklabels(), visible=False)
            ax.set_xlim(2000,25000)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(4000))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(500))
            ax.set_xlabel('Wavelength $({\\rm \\AA})$')
            pl.subplots_adjust(hspace=0, wspace=0)


def extend_variance_covariance(refwaveUV=3200, minwave=1700,
                               refwaveIR=8500, maxwave=25000,
                               modeldir=__MODELDIR__,
                               salt2indir='salt2-4', salt2outdir='salt2-extended',
                               extrapolate=True, showplots=False, plotphase=0):
    """ extrapolate the *lc* and *spec* .dat files for SALT2
     using a flatline extension of the red tail to 2.5 microns
     """
    salt2inmodeldir = os.path.join(modeldir, salt2indir)
    salt2outmodeldir = os.path.join(modeldir, salt2outdir)

    filelist = ['salt2_lc_relative_variance_0.dat',
                'salt2_lc_relative_variance_1.dat',
                'salt2_lc_relative_covariance_01.dat']
    titlelist = ['Var0', 'Var1', 'Covar01']

    iax = 0
    if showplots:
        fig = pl.gcf()
        fig.clf()

    for filename, title in zip(filelist, titlelist):
        infile = os.path.join(salt2inmodeldir, filename)
        outfile = os.path.join(salt2outmodeldir, filename)
        timeseries = TimeSeriesGrid(infile)

        if extrapolate:
            timeseries.extrapolate_flatline(
                outfilename=outfile, refwave=refwaveIR, maxwave=maxwave)
            timeseries2 = TimeSeriesGrid(outfile)
            timeseries2.extrapolate_flatline_uv(outfilename=outfile)

        if showplots:
            iax+=1
            ax = pl.subplot(3, 1, iax)
            timeseriesNew = TimeSeriesGrid(outfile)
            for phase, color in zip([-5, 0, 5, 15],['r','darkorange','g','b']):
                timeseries.plot_at_phase(phase, color='0.5', lw=3.2,
                                         ls='-', marker=' ',
                                         label='__nolabel__')
                timeseriesNew.plot_at_phase(phase, color=color, lw=1,
                                            ls='--', marker=' ',
                                            label='t=%i'%phase)
            #ax.text(0.90, 0.9, title, transform=ax.transAxes,
            #        fontsize='large', ha='right', va='top')
            #if iax==1:
            #    ax.set_ylim(-0.2, 6)
            #elif iax in [2,3,4]:
            #    ax.set_ylim(-5e-4,3e-3)
            #else:
            #    ax.set_ylim(-5e-6,5e-5)
            if iax==1:
                pl.legend(loc='upper right', ncol=2)
            if iax<3:
                pl.setp(ax.get_xticklabels(), visible=False)
                ax.set_yscale('log')
            else:
                ax.set_ylim(-2.2e-3, 1.2e-3)
            #pl.setp(ax.get_yticklabels(), visible=False)
            ax.set_xlim(2000,12500)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(4000))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(500))
            #ax.yaxis.set_major_locator(ticker.MaxNLocator(3))
            #ax.yaxis.set_minor_locator(ticker.MaxNLocator(6))
            ax.text(0.5, 0.95, title,
                    ha='right', va='top', transform=ax.transAxes)
            ax.set_xlabel('Wavelength $(\\rm{\\AA})$')
            pl.subplots_adjust(hspace=0, wspace=0)


def extend_info_file(wref=8500, modeldir=__MODELDIR__,
                    salt2irindir='salt2-4', salt2iroutdir='salt2-extended',
                    extrapolate=False, showplots=True):
    """ Update the SALT2.INFO file to reflect the modified color law
    """
    # TODO : read in the pre-existing SALT2.INFO file, and only update
    # the color law polynomial parameters.
    # salt2modeldir = os.path.join(modeldir, salt2indir)
    salt2irmodeldir = os.path.join(modeldir, salt2irindir)

    # read in the salt2_color_correction.dat file
    # and repeat the fit parameters into the .INFO file:
    colorlawfile = os.path.join(salt2irmodeldir,
                                'salt2_color_correction.dat')
    fin = open(colorlawfile, 'r')
    words = fin.read().split()
    fin.close()

    # Get colorlaw coeffecients.
    npoly = int(words[0])
    colorlaw_coeffs = [float(word) for word in words[1: 1 + npoly]]

    # Look for keywords in the rest of the file.
    version = 0
    colorlaw_range = [3000., 7000.]
    for i in range(1+npoly, len(words)):
        if words[i] == 'Salt2ExtinctionLaw.version':
            version = int(words[i+1])
        if words[i] == 'Salt2ExtinctionLaw.min_lambda':
            colorlaw_range[0] = float(words[i+1])
        if words[i] == 'Salt2ExtinctionLaw.max_lambda':
            colorlaw_range[1] = float(words[i+1])
    COLORLAWLINE = 'COLORCOR_PARAMS: %i %i %i ' % (
        colorlaw_range[0], colorlaw_range[1], npoly)
    for i in range(npoly):
        COLORLAWLINE += ' %.8f' % colorlaw_coeffs[i]

    outinfo = os.path.join(__MODELDIR__,salt2iroutdir, 'SALT2.INFO')
    fout = open(outinfo, 'w')
    print >> fout, """
# ====================================== 
# SALT2.INFO from SALT2.JLA-B14
# 
RESTLAMBDA_RANGE: 2000. 25000.  # extended using snsedextend tools 
COLORLAW_VERSION: 1 
""" +\
COLORLAWLINE + \
"""
COLOR_OFFSET: 0.0 

MAG_OFFSET: 0.27 # to get B-band mag from cosmology fit (Nov 23, 2011) 

SEDFLUX_INTERP_OPT: 2 # 1=>linear, 2=>spline 
ERRMAP_INTERP_OPT: 1  # 0=snake off; 1=>linear 2=>spline 
ERRMAP_KCOR_OPT: 1    # 0=ignore salt2_color_dispersion.dat; 1=> use it to define K-corr. uncertainty vs lambda

MAGERR_FLOOR: 0.005 # don;t allow smaller error than this 
MAGERR_LAMOBS: 0.1 2000 4000 # magerr minlam maxlam 
MAGERR_LAMREST: 0.1 100 200 # magerr minlam maxlam 

SIGMA_INT: 0.106 # used in simulation 
"""
    fout.close()
    print("SALT2.INFO file updated: %s"%outinfo)
    return
