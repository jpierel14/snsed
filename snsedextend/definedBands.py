import numpy as np
import os,sncosmo
from .utils import __dir__

u_wave,u_trans=np.loadtxt(os.path.join(__dir__,'data','bands','uBand','tophatU.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(u_wave,u_trans,name='tophatU'),force=True)

b_wave,b_trans=np.loadtxt(os.path.join(__dir__,'data','bands','bBand','tophatB.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(b_wave,b_trans,name='tophatB'),force=True)

v_wave,v_trans=np.loadtxt(os.path.join(__dir__,'data','bands','vBand','tophatV.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(v_wave,v_trans,name='tophatV'),force=True)

r_wave,r_trans=np.loadtxt(os.path.join(__dir__,'data','bands','rBand','tophatR.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(r_wave,r_trans,name='tophatR'),force=True)

i_wave,i_trans=np.loadtxt(os.path.join(__dir__,'data','bands','iBand','tophatI.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(i_wave,i_trans,name='tophatI'),force=True)

j_wave,j_trans=np.loadtxt(os.path.join(__dir__,'data','bands','jBand','tophatJ.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(j_wave,j_trans,name='tophatJ'),force=True)

h_wave,h_trans=np.loadtxt(os.path.join(__dir__,'data','bands','hBand','tophatH.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(h_wave,h_trans,name='tophatH'),force=True)

k_wave,k_trans=np.loadtxt(os.path.join(__dir__,'data','bands','kBand','tophatK.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(k_wave,k_trans,name='tophatK'),force=True)

for f in ['J','H','Ks']:
    wave,trans=np.loadtxt(os.path.join(__dir__,'data','bands',f[0].lower()+'Band','paritel'+f+'.dat'),unpack=True)
    wave*=10000
    sncosmo.registry.register(sncosmo.Bandpass(wave,trans,name='paritel::'+f.lower()),force=True)
