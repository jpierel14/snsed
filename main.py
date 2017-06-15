import snsedextend
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import glob,os,sys



#colorTable=snsedextend.curveToColor('lcs/lc_'+filelist[i]+'.dat',['V-H'],bandFit='V',zpsys='Vega',constants={'z':redshift[i],'hostr_v':3.1,'mwr_v':3.1},dust='OD94Dust',effect_frames=['rest','obs'],effect_names=['host','mw'])

#colorTable=snsedextend.curveToColor('lcs/lc_2007av.dat',['V-H'],bandFit='V',zpsys='Vega',constants={'z':.004650},singleBand=True)
'''
import glob
import sncosmo
import os
files=glob.glob('/Users/jpierel/rodney/snsedextend/lcs_clipped/tables/uncorr*')
zpmag=sncosmo.get_magsystem('Vega')
bands={v:k for k,v in snsedextend.extrapolate._filters.iteritems()}
sne=[]
uncorr=[]
corr=[]
for file in files:
    #band,date,mag,=np.loadtxt(file,unpack=True,skiprows=1,usecols=(1,2,3))
    #lc=snsedextend.extrapolate._standardize(sncosmo.read_lc(file))
    color=sncosmo.read_lc(file)
    col=[x for x in color.colnames() if '-' in x][0]
    sne.append('SN'+file[file.rfind('_'):-4])
    uncorr.append(np.mean(color[col]))
    if 'U' not in lc['band'] and 'J' not in lc['band'] and 'H' not in lc['band'] and 'K' not in lc['band']:
        os.remove(file)
    
    try:
        lc.remove_column('zp')
    except:
        pass
    try:
        lc.remove_column('zpsys')
    except:
        pass
    try:
        lc.remove_column('flux')
    except:
        pass
    try:
        lc.remove_column('fluxerr')
    except:
        pass
    
    try:
        band=[bands[b] for b in lc['band']]
    except:
        continue
    lc.remove_column('band')
    lc['band']=band
    
    sncosmo.write_lc(lc,file)
    
    if 'bessellux' in lc['band']:
        blue='U'
        num=0
        for band in [x for x in np.unique(lc['band']) if x in ['bessellb','bessellv','bessellr','sdss::g','sdss::r']]:
            temp=len(lc[lc['band']==band])
            if temp>num:
                maxBand=band
                num=temp
        red=maxBand
        if 'cspk' in lc['band']:
            num=0
            for band in [x for x in np.unique(lc['band']) if x in ['bessellb','bessellv','bessellr','sdss::g','sdss::r']]:
                temp=len(lc[lc['band']==band])
                if temp>num:
                    maxBand=band
                    num=temp
            blue2=maxBand
            num=0
            for band in [x for x in np.unique(lc['band']) if x in ['cspjs','csphs','cspk']]:
                temp=len(lc[lc['band']==band])
                if temp>num:
                    maxBand=band
                    num=temp
            red2=maxBand
            print(os.path.basename(file),blue+'-'+red,blue2+'-'+red2)
            continue
    else:
        try:
            num=0
            for band in [x for x in np.unique(lc['band']) if x in ['bessellb','bessellv','bessellr','sdss::g','sdss::r']]:
                temp=len(lc[lc['band']==band])
                if temp>num:
                    maxBand=band
                    num=temp
            blue=maxBand
            num=0
            for band in [x for x in np.unique(lc['band']) if x in ['cspjs','csphs','cspk']]:
                temp=len(lc[lc['band']==band])
                if temp>num:
                    maxBand2=band
                    num=temp
            red=maxBand2
        except:
            os.remove(file)
            continue
    print(os.path.basename(file),blue+'-'+red)
    
    if lc['band'][0] in snsedextend.extrapolate._filters.values():
        continue
    elif 'U' not in lc['band'] and set(np.unique(lc['band']))==set(['J','H','K']):
        os.remove(file)
        continue

    zp=[]
    for band in lc['band']:
        zp.append(zpmag.band_flux_to_mag(1,snsedextend.extrapolate._filters[band]))
    lc['zp']=zp
    lc['zpsys']=['Vega' for i in range(len(lc))]
    band=[snsedextend.extrapolate._filters[band] for band in lc['band']]
    lc.remove_column('band')
    lc['band']=band
    lc=snsedextend.extrapolate._mag_to_flux(lc)
    if max(lc['time'])-min(lc['time'])<50:
        continue
    else:
        #maxes={(band,lc['time'][lc[lc['band']==band]['flux']==max(lc[lc['band']==band]['flux'])]) for band in lc['band']}
        maxTime=lc['time'][lc['flux']==max(lc['flux'])]
    #for row in lc.rows():
    temp=lc[lc['band']=='bessellux']
    temp2=lc[lc['band']==]
    lc=lc[lc['time']<=maxTime+50]
    lc=lc[lc['time']>=maxTime-50]
    sncosmo.write_lc(lc,file)
    #plt.show()
    

#filelist=['2007rt','2008in','2008ip','2009ay','2010bq','2007av']#,'2007aa']
#redshift=[.022365,.005224,.015124,.022182,.030988,.004650]#.004887]
host_ebvs={'SN2006kd':0.035,'SN2006ca':0.104,'SN2006cd':-0.596,'SN2006it':0.075,'SN2007aa':0.146,'SN2007av':0.086,
           'SN2007rt':0.888,'SN2008aj':-0.736,'SN2008bj':-0.002,'SN2008bn':-0.116,'SN2008gm':-0.128,'SN2008in':-0.005,
           'SN2008ip':-0.793,'SN2010al':-0.944,'SN2010bq':-0.843}
mw_ebvs={'SN2006kd':0.036,'SN2006ca':0.106,'SN2006cd':-0.631,'SN2006it':0.077,'SN2007aa':0.147,'SN2007av':0.086,
         'SN2007rt':-1.000,'SN2008aj':-0.750,'SN2008bj':-0.001,'SN2008bn':-0.121,'SN2008gm':-0.133,'SN2008in':-0.005,
         'SN2008ip':-0.019,'SN2010al':1.000,'SN2010bq':1.000}
filelist=[os.path.basename(file) for file in glob.glob('/Users/jpierel/rodney/snsedextend/lcs_clipped/*.dat')]
#filelist=['lc_2007rt.dat','lc_2010al.dat','lc_2010bq.dat']
sn,redshift=np.loadtxt('lcs_clipped/redshift.ref',dtype='str',unpack=True)
redshift={'lc_'+sn[i].replace('SN',''):redshift[i] for i in range(len(sn))}
sn,colors=np.loadtxt('lcs_clipped/color.ref',dtype='str',unpack=True)
colors={sn[i].replace('SN',''):colors[i] for i in range(len(sn))}
#snsedextend.extendNon1a(colorTable,{'V':'bessellv','B':'f435w'},sedlist='SDSS-019323.SED',verbose=True)
#effective wavelength, blueward of 4000, redward of 9000, right anchor point is going to be right edge of jwst 55000, left anchor point 1200
#fig,axs=plt.subplots(nrows=3,ncols=2)
row=0
col=0
for i in range(len(filelist)):
    print(filelist[i])
    #ax=axs[row,col]
    fig=plt.figure()
    ax=plt.gca()
    colorTable=snsedextend.curveToColor('lcs_clipped/'+filelist[i],[colors[filelist[i][:-4]]],zpsys='Vega',bounds={'hostebv':(-1,1),'mwebv':(-1,1)},constants={'z':redshift[filelist[i][:-4]],'hostr_v':3.1,'mwr_v':3.1},dust='CCM89Dust',effect_frames=['rest','obs'],effect_names=['host','mw'])
    print(colorTable)
    ax.errorbar(np.array(colorTable['time']),np.array(colorTable[colors[filelist[i][:-4]]]),yerr=np.array(colorTable[colors[filelist[i][:-4]].replace('-','')+'_err']),fmt='x')
    ax.set_title(filelist[i][3:])
    if i<2:
        row+=1
    elif i==2:
        row=0
        col+=1
    else:
        row+=1

#fig.suptitle('V-H Colors')
    ax.invert_yaxis()
    plt.title(colors[filelist[i][:-4]]+' Color')
    fig.text(0.5, 0.01, 'Time (Since Peak)', ha='center')
    fig.text(0.01, 0.5, 'Color Magnitude (Vega)', va='center', rotation='vertical')
    plt.savefig("lcs_clipped/plots/"+filelist[i][:-4]+".pdf",format='pdf')

plt.tight_layout()
plt.show() 
'''
#snsedextend.extendNon1a(colorTable,{'V':'bessellv','B':'f435w'},sedlist='SDSS-019323.SED',verbose=True)
#make sure to add in dust

from astropy.table import vstack
from astropy.io import ascii
import sncosmo

dir='modjaz'
type=['II']
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

typeIColors=None

for k in colors:
    if len(colors[k])>3:
        colors[k]=colors[k].split(',')
    elif not isinstance(colors[k],(list,tuple)):
        colors[k]=[colors[k]]
for i in range(len(filelist)):
    if filelist[i]=='lc_2005hg_clipped.dat':#typ[filelist[i][:-12]] in type and filelist[i] not in ['lc_2006fo_clipped.dat','lc_2006jc_clipped.dat','lc_2004ao_clipped.dat','lc_2007D_clipped.dat','lc_2005bf_clipped.dat','lc_2005nb_clipped.dat','lc_2006ld_clipped.dat']:

        print(filelist[i])
        colorTable=snsedextend.curveToColor(os.path.join(dir,filelist[i]),colors[filelist[i][:-12]],snType=typ[filelist[i][:-12]],zpsys='Vega',
                                            bounds={'hostebv':(-1,1),'t0':(peaks[filelist[i][:-12]]-5,peaks[filelist[i][:-12]]+5)},
                                            constants={'z':redshift[filelist[i][:-12]],'hostr_v':3.1,'mwr_v':3.1,'mwebv':dust[filelist[i][:-12]]},
                                            dust='CCM89Dust',effect_frames=['rest','obs'],effect_names=['host','mw'])
        snsedextend.extendNon1a(colorTable,sedlist='SDSS-000020.SED',verbose=True)

        continue
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
        if typeIColors:
            typeIColors=vstack([typeIColors,colorTable])
        else:
            typeIColors=colorTable

    #snsedextend.extendNon1a(colorTable,sedlist='SDSS-0 13449.SED',verbose=True)
typeIColors.sort('time')
#sncosmo.write_lc(allColors,'lcs_clipped2/tables/allColors.dat')
ascii.write(typeIColors,os.path.join(dir,'type'+type[0],'tables','all'+type[0]+'Colors.dat'))
