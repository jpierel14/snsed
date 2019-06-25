from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

#import seaborn as sns
import sys,os,warnings,matplotlib
#from astropy.io import ascii
from contextlib import contextmanager

import pymc3 as pm
import logging

warnings.filterwarnings('ignore')
logger=logging.getLogger('pymc3')
logger.disabled=True

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def create_poly_modelspec(k=1):
    '''
    Convenience function:
    Create a polynomial modelspec string for patsy
    '''
    return ('y ~ 1 + x ' + ' '.join(['+ np.power(x,{})'.format(j)
                                     for j in range(2,k+1)])).strip()




def run_models(df, upper_order=5):
    '''
    Convenience function:
    Fit a range of pymc3 models of increasing polynomial complexity.
    Suggest limit to max order 5 since calculation time is exponential.
    '''

    models, traces = OrderedDict(), OrderedDict()
    for k in range(1,upper_order+1):

        nm = 'k{}'.format(k)
        fml = create_poly_modelspec(k)
        #with suppress_stdout():
        with pm.Model() as models[nm]:
            pm.glm.GLM.from_formula(fml, df, family=pm.glm.families.Normal())

            # For speed, we're using Metropolis here

            traces[nm] = pm.sample(1000, pm.Metropolis(),progressbar=False)[200::5]
    return models, traces

from matplotlib.offsetbox import AnchoredText


def getModels(models, traces, rawdata,table, xlims,
              datamodelnm='linear', modelnms=None,bestModel=None,allDat=None,typ=None,bic=None,color=None,savefig=False):
    if savefig:
        colors=['b','k','y','orange','cyan','violet','r','g']
        import seaborn as sns
        fig=plt.figure()
        ax=fig.gca()


        if typ =='II':
            d='hicken'
        elif typ=='IIn':
            d='hicken/IIn'
        else:
            d='bianco'
        sne,types=np.loadtxt(os.path.join(d,'type.ref'),dtype='str',unpack=True)
        typeDict={'SN_'+sne[i][3:]:types[i] for i in range(len(sne))}
        colorDict={np.unique([typeDict[x] for x in table['SN']])[i]:colors[i] for i in range(len(np.unique([typeDict[x] for x in table['SN']])))}
        table['snTypes']=[typeDict[x] for x in table['SN']]

    for modelnm in modelnms:
        if modelnm==bestModel:
            ## Get traces and calc posterior prediction for npoints in x
            npoints = 100
            mdl = models[modelnm]
            trc = pm.trace_to_dataframe(traces[modelnm][-1000:])
            trc = trc[[str(v) for v in mdl.cont_vars[:-1]]]

            ordr = int(modelnm[-1:])
            x = np.linspace(xlims[0], xlims[1], npoints).reshape((npoints,1))
            pwrs = np.ones((npoints,ordr+1)) * np.arange(ordr+1)
            X = x ** pwrs
            cr = np.dot(X,trc.T)

            ## Calculate credible regions and plot over the datapoints
            dfp = pd.DataFrame(np.percentile(cr,[2.5, 25, 50, 75, 97.5], axis=1).T
                               ,columns=['025','250','500','750','975'])
            dfp['x'] = x
            if savefig:

                '''
                pal = sns.color_palette('Reds')
                ax.fill_between(dfp['x'], dfp['025'], dfp['975'], alpha=0.5
                                ,color=pal[1], label='CR 95%')
                ax.fill_between(dfp['x'], dfp['250'], dfp['750'], alpha=0.5
                                ,color=pal[4], label='CR 50%')
                ax.plot(dfp['x'], dfp['500'], alpha=0.6, color=pal[5], label='Median')
                _ = ax.legend(loc='lower right',fontsize=10)
                _ = ax.set_xlim(xlims)
                #_ = ax1d.errorbar(allDat['x'],allDat['y'],yerr=allDat['error'],fmt='.',color='red')
                types=[]
                for type in np.unique([typeDict[x] for x in rawdata['names']]):
                    types.append(ax.errorbar(table['time'][table['snTypes']==type],table['mag'][table['snTypes']==type],yerr=table['magerr'][table['snTypes']==type],fmt='.',color=colorDict[type]))
                _ = sns.regplot(x='x', y='y', data=rawdata, fit_reg=False
                                ,scatter_kws={'alpha':0.7,'s':100, 'lw':2,'edgecolor':'w'}, ax=ax)
                ax.set_xlabel('Days After Peak',size=14)
                ax.set_ylabel('Color (Magnitude)', size=14)
                plt.figlegend(types,[str(x) for x in np.unique([typeDict[x] for x in rawdata['names']])],bbox_to_anchor=(.7,.23))
                
                plt.savefig('type'+typ+'_'+out.upper()+'_fits.pdf',format='pdf')
                plt.close()
                '''
                if color[0]=='U':
                    out='U'
                    #ax.set_title('Posterior Predictive Fits -- Data: U-B, Type {} -- Best Model: Order {}'.format(
                    #    typ, bestModel[1]), fontsize=12)
                else:
                    out=color[-1]
                    #ax.set_title('Posterior Predictive Fits -- Data: r-{}, Type {} -- Best Model: Order {}'.format(
                    #   color[-1],typ, bestModel[1]), fontsize=12)
                #plt.savefig('type'+typ+'_'+out.upper()+'_fits.pdf',format='pdf')
                #plt.close()
                cmap = matplotlib.cm.get_cmap('jet')

                fig=plt.figure()
                ax=fig.gca()
                sne=np.unique(np.asarray(table['SN']))
                snColorDict={sne[i]:colors[i] for i in range(len(sne))}
                allSne=[]
                #for sn in sne:
                #    #allSne.append(ax.scatter(table['time'][table['SN']==sn],table['mag'][table['SN']==sn],color=snColorDict[sn]))
                #norm = matplotlib.colors.Normalize(vmin=min(table['mag']), vmax=max(table['mag']))
                allVR=[]
                shapes=['.','*','o','^','+','8','s']

                b=0
                for sn in sne:
                    ax.scatter(table['time'][table['SN']==sn],table['mag'][table['SN']==sn],label=sn)
                    b+=1
                    ax.errorbar(table['time'][table['SN']==sn],table['mag'][table['SN']==sn],yerr=table['magerr'][table['SN']==sn],fmt=None, marker=None, mew=0,lw=.5,color='k',alpha=.4,label=None)#,c=cmap,fmt='.')
                    temp=ax.scatter(table['time'][table['SN']==sn],table['mag'][table['SN']==sn],marker=shapes[b],label=sn)
                    b+=1
                    ax.errorbar(table['time'][table['SN']==sn],table['mag'][table['SN']==sn],yerr=table['magerr'][table['SN']==sn],fmt=None, marker=None, mew=0,lw=.5,color='k',alpha=.4,label=None)#,c=cmap,fmt='.')
                plt.colorbar(temp)
                ax.legend(loc='lower right')
                ax.plot(dfp['x'], dfp['500'], color='k', label='Median')
                ax.plot(dfp['x'],dfp['500']+np.std(table['mag']),color='r',linestyle='--')
                ax.plot(dfp['x'],dfp['500']-np.std(table['mag']),color='b',linestyle='--')
                #plt.figlegend(allSne,sne,bbox_to_anchor=(.9,.4))
                ax.set_xlabel('Days After Peak',size=14)
                ax.set_ylabel('Color (Magnitude)', size=14)
                ax.set_title('Color Curve by Supernova: Type {} -- Color={}'.format(typ,color))
                
                plt.savefig('type'+typ+'_'+out.upper()+'_fits_SNe.pdf',format='pdf')
                #sys.exit()
            return(dfp)

def BICrun(table,color=None,type='II',verbose=False,savefig=False):
    #table=table[table['time']<=50]
    modelList=['k1','k2']
    #print(logging.Logger.manager.loggerDict)
    temp=pd.DataFrame({'x':np.array(table['time']),'y':np.array(table['mag']),'error':np.array(table['magerr']),'names':table['SN']})
    temp_xlims = (temp['x'].min() - np.ptp(temp['x'])/10,temp['x'].max() + np.ptp(temp['x'])/10)
    models_lin,traces_lin=run_models(temp,upper_order=2)
    dfdic = pd.DataFrame(index=modelList, columns=['dic','waic'])
    dfdic.index.name = 'model'
    for nm in dfdic.index:
        #dfdic.loc[nm, 'dic'] = pm.stats.dic(traces_lin[nm], models_lin[nm])
        dfdic.loc[nm, 'waic'] = pm.stats.waic(traces_lin[nm], models_lin[nm])[0]

    dfwaic = pd.DataFrame(index=modelList, columns=['lin'])
    dfwaic.index.name = 'model'

    for nm in dfwaic.index:
        dfwaic.loc[nm, 'lin'] = pm.stats.waic(traces_lin[nm], models_lin[nm])[0]

    best=dfwaic[dfwaic['lin']==np.min(dfwaic['lin'])].index[0]
    dfwaic = pd.melt(dfwaic.reset_index(), id_vars=['model'], value_name='waic')

    return(getModels(models_lin,traces_lin,temp,table,temp_xlims, bestModel=best,modelnms=modelList,bic=dfwaic,typ=type,color=color,savefig=savefig))

















