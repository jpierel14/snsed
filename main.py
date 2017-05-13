import snsedextend
from astropy.table import Table
t=[0,5,10,15,20]
c=[.5,.5,.5,.5,.5]
d=[1,1,1,1,1]
colorTable=Table([t,c,d],names=['mjd','v-h','u-b'],masked=True)
#snsedextend.extendNon1a(colorTable,{'V':'bessellv','B':'f435w'},sedlist='SDSS-019323.SED',verbose=True)
#effective wavelength, blueward of 4000, redward of 9000, right anchor point is going to be right edge of jwst 55000, left anchor point 1200
colorTable=snsedextend.curveToColor('red.lc',['U-V',"U-r"],zpsys='Vega',constants={'z':.015511})
print(colorTable)
snsedextend.extendNon1a(colorTable,{'V':'bessellv','B':'f435w'},sedlist='SDSS-019323.SED',verbose=True)
#make sure to add in dust