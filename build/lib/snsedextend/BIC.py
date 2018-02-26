from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sns
import sys,os,warnings
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

def plot_posterior_cr(models, traces, rawdata, xlims,
                      datamodelnm='linear', modelnms=None,bestModel=None,allDat=None,typ=None,bic=None):
    '''
    Convenience function:
    Plot posterior predictions with credible regions shown as filled areas.
    '''
    ax={}
    #f, ((ax[1],ax[2]),(ax[3],ax[4])) = plt.subplots(nrows=2,ncols=2,sharex=True,sharey=True, figsize=(16,12))
    f, (ax[1],ax[2]) = plt.subplots(nrows=2,ncols=1,sharex=True,sharey=True, figsize=(16,12))
    i=1
    f.text(0.5, 0.05, 'Days After Peak', ha='center', va='center',size=22)
    f.text(0.05, 0.5, 'Color (Magnitude)', ha='center', va='center', rotation='vertical',size=22)
    if datamodelnm.upper()=='UV':
        f.suptitle('Posterior Predictive Fits -- Data: U-Optical, Type {} -- Best Model: Order {}'.format(
            typ, bestModel[1]), fontsize=24)
    else:
        f.suptitle('Posterior Predictive Fits -- Data: Optical-{}, Type {} -- Best Model: Order {}'.format(
            datamodelnm.upper(),typ, bestModel[1]), fontsize=24)
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
            return(dfp)

def getModels(models, traces, rawdata, xlims,
            datamodelnm='linear', modelnms=None,bestModel=None,allDat=None,typ=None,bic=None):

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
            return(dfp)

def BICrun(table,type='II',verbose=False):

    modelList=['k1','k2']
    #print(logging.Logger.manager.loggerDict)
    temp=pd.DataFrame({'x':np.array(table['time']),'y':np.array(table['mag']),'error':np.array(table['magerr'])})
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


    return(getModels(models_lin,traces_lin,temp,temp_xlims, bestModel=best,modelnms=modelList,bic=dfwaic))

















