import snsedextend
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import glob,os,sys
from astropy.table import vstack
from astropy.io import ascii
import sncosmo
sndataroot = os.environ['SNDATA_ROOT']
dir='hicken'
type=['II','IIP']
filelist=[os.path.basename(file) for file in glob.glob(os.path.join(dir,'*clipped.dat'))]
#filelist=['lc_2002bx.dat','lc_2006ca.dat','lc_2006cd.dat','lc_2006it.dat','lc_2007aa.dat','lc_2007av.dat','lc_2008bj.dat','lc_2008bn','lc_2008in.dat','lc_2009ay.dat','lc_2009kn.dat',]
#IIn=['lc_2008ip.dat','lc_2010bq.dat']
sn,redshift=np.loadtxt(os.path.join(dir,'redshift.ref'),dtype='str',unpack=True)
redshift={sn[i]:redshift[i] for i in range(len(sn))}
sn,colors=np.loadtxt(os.path.join(dir,'color.ref'),dtype='str',unpack=True)
colors={sn[i]:colors[i] for i in range(len(sn))}
sn,dusts=np.loadtxt(os.path.join(dir,'dust.ref'),dtype='str',unpack=True)
dust={sn[i]:dusts[i] for i in range(len(sn))}
sn,types=np.loadtxt(os.path.join(dir,'type.ref'),dtype='str',unpack=True)
typ={sn[i]:types[i] for i in range(len(sn))}
sne,peaks=np.loadtxt(os.path.join(dir,'timeOfPeak.dat'),dtype='str',unpack=True)
peaks={sne[i]:float(peaks[i]) for i in range(len(sne))}

typeColors=None

for k in colors:
    if len(colors[k])>3:
        colors[k]=colors[k].split(',')
    elif not isinstance(colors[k],(list,tuple)):
        colors[k]=[colors[k]]


for i in range(len(filelist)):
    if typ[filelist[i][:-12]] in type and filelist[i] not in ['lc_2005kl_clipped.dat','lc_2005az_clipped.dat','lc_2008aq_clipped.dat','lc_2007av_clipped.dat']:#['lc_2006fo_clipped.dat','lc_2006jc_clipped.dat','lc_2004ao_clipped.dat','lc_2007D_clipped.dat','lc_2005bf_clipped.dat','lc_2005nb_clipped.dat','lc_2006ld_clipped.dat']:

        print(filelist[i])
        if filelist[i] not in ['lc_2008bj_clipped.dat','lc_2008bn_clipped.dat']:#,'lc_2008in_clipped.dat']:
            continue
        else:
            print(colors[filelist[i][:-12]],typ[filelist[i][:-12]])
            colorTable=snsedextend.curveToColor(os.path.join(dir,filelist[i]),colors[filelist[i][:-12]],snType=typ[filelist[i][:-12]],zpsys='Vega',
                                                bounds={'hostebv':(-1,1),'t0':(peaks[filelist[i][:-12]]-5,peaks[filelist[i][:-12]]+5)},
                                                constants={'z':redshift[filelist[i][:-12]],'hostr_v':3.1,'mwr_v':3.1,'mwebv':dust[filelist[i][:-12]]},
                                                dust='CCM89Dust',effect_frames=['rest','obs'],effect_names=['host','mw'])
        #else:
        #    continue
        '''
        for color in colors[filelist[i][:-12]]:
            fig=plt.figure()
            ax=plt.gca()
            ax.errorbar(np.array(colorTable['time'][colorTable[color].mask==False]),np.array(colorTable[color][colorTable[color].mask==False]),
                        yerr=np.array(colorTable[color.replace('-','')+'_err'][colorTable[color].mask==False]),fmt='x')
            ax.set_title(filelist[i][3:])
            ax.invert_yaxis()
            plt.title(color+' Color')
            fig.text(0.5, 0.01, 'Time (Since Peak)', ha='center')
            fig.text(0.01, 0.5, 'Color Magnitude (Vega)', va='center', rotation='vertical')
            plt.savefig(os.path.join(dir,"type"+type[0],"plots",filelist[i][:-12]+"_"+color[0]+color[-1]+".pdf"),format='pdf')
            plt.close()
        
        if typeColors:
            typeSNe=typeSNe+[filelist[i][:-12].replace('lc','SN') for j in range(len(colorTable))]

            typeColors=vstack([typeColors,colorTable])
        else:
            typeSNe=[filelist[i][:-12].replace('lc','SN') for j in range(len(colorTable))]
            typeColors=colorTable
        '''

            #snsedextend.extendNon1a(colorTable,sedlist='SDSS-0 13449.SED',verbose=True)
sys.exit()
typeColors['SN']=typeSNe

from astropy.table import MaskedColumn
#typeColors=ascii.read(os.path.join(dir,'type'+type[0],'tables','all'+type[0]+'Colors.dat'))
#vColors=ascii.read(os.path.join('SEDs','typeII','vColors.dat'))
#temp=MaskedColumn([np.nan for i in range(len(typeColors))],name='V-r',mask=[True for i in range(len(typeColors))])
#temp2=MaskedColumn([np.nan for i in range(len(typeColors))],name='Vr_err',mask=[True for i in range(len(typeColors))])
#typeColors['V-r']=temp
#typeColors['Vr_err']=temp2
#for i in range(len(vColors)):
#    typeColors.add_row(np.append(vColors['time'][i],np.append([np.nan for j in range(len(typeColors.colnames)-3)],np.append(vColors['V-r'][i],vColors['Vr_err'][i]))),mask=np.append([False],np.append([True for k in range(len(typeColors.colnames)-3)],[False,False])))
typeColors.sort('time')
ascii.write(typeColors,os.path.join(dir,'type'+type[0],'tables','all'+type[0]+'Colors.dat'))
sys.exit()
seds=np.loadtxt(os.path.join('/Users','jpierel','rodney','snsedextend','SEDs','NON1A.LIST'),dtype='str',unpack=True)
sedlist=[seds[3][i]+'.SED' for i in range(len(seds[3])) if seds[2][i] in type and seds[3][i]+'.SED' in [os.path.basename(x) for x in glob.glob(os.path.join(sndataroot,'snsed','NON1A','*.SED'))]]

snsedextend.extendNon1a(typeColors,colors=['U-B','r-J','r-H','r-K'],sedlist=sedlist,zpsys='Vega',verbose=True)
    #snsedextend.extendNon1a(colorTable,sedlist='SDSS-0 13449.SED',verbose=True)
#typeColors.sort('time')
#sncosmo.write_lc(typeColors,'lcs_clipped2/tables/allColors.dat')
