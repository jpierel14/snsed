import numpy as np
import os,sncosmo
#os.path.dirname(snsedextend)
u_wave,u_trans=np.loadtxt(os.path.join('snsedextend','data','bands','uBand','tophatU.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(u_wave,u_trans,name='tophatU'),force=True)

b_wave,b_trans=np.loadtxt(os.path.join('snsedextend','data','bands','bBand','tophatB.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(b_wave,b_trans,name='tophatB'),force=True)

v_wave,v_trans=np.loadtxt(os.path.join('snsedextend','data','bands','vBand','tophatV.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(v_wave,v_trans,name='tophatV'),force=True)

r_wave,r_trans=np.loadtxt(os.path.join('snsedextend','data','bands','rBand','tophatR.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(r_wave,r_trans,name='tophatR'),force=True)

i_wave,i_trans=np.loadtxt(os.path.join('snsedextend','data','bands','iBand','tophatI.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(i_wave,i_trans,name='tophatI'),force=True)

j_wave,j_trans=np.loadtxt(os.path.join('snsedextend','data','bands','jBand','tophatJ.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(j_wave,j_trans,name='tophatJ'),force=True)

h_wave,h_trans=np.loadtxt(os.path.join('snsedextend','data','bands','hBand','tophatH.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(h_wave,h_trans,name='tophatH'),force=True)

k_wave,k_trans=np.loadtxt(os.path.join('snsedextend','data','bands','kBand','tophatK.dat'),unpack=True)
sncosmo.registry.register(sncosmo.Bandpass(k_wave,k_trans,name='tophatK'),force=True)

