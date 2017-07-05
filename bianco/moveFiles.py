import os
import numpy as np
import sncosmo

lcs=sncosmo.read_lc('../lc.standardsystem.sesn_allphot.dat')
sne=np.unique(lcs['Name'])
for s in sne:
	sncosmo.write_lc(lcs[lcs['Name']==s],'lc_'+s+'.dat')