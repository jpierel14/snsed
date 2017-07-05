from collections import OrderedDict
from time import time

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.optimize import fmin_powell
from scipy import integrate

import pymc3 as pm
import theano as thno
import theano.tensor as T


# configure some basic options
sns.set(style="darkgrid", palette="muted")
pd.set_option('display.notebook_repr_html', True)
plt.rcParams['figure.figsize'] = 12, 8
rndst = np.random.RandomState(0)

def generate_data(n=20, p=0, a=1, b=1, c=0, latent_sigma_y=20):
    '''
    Create a toy dataset based on a very simple model that we might
    imagine is a noisy physical process:
        1. random x values within a range
        2. latent error aka inherent noise in y
        3. optionally create labelled outliers with larger noise

    Model form: y ~ a + bx + cx^2 + e

    NOTE: latent_sigma_y is used to create a normally distributed,
    'latent error' aka 'inherent noise' in the 'physical process'
    generating thses values, rather than experimental measurement error.
    Please don't use the returned `latent_error` values in inferential
    models, it's returned in e dataframe for interest only.
    '''

    df = pd.DataFrame({'x':rndst.choice(np.arange(100),n,replace=False)})

    ## create linear or quadratic model
    df['y'] = a + b*(df['x']) + c*(df['x'])**2

    ## create latent noise and marked outliers
    df['latent_error'] = rndst.normal(0,latent_sigma_y,n)
    df['outlier_error'] = rndst.normal(0,latent_sigma_y*10,n)
    df['outlier'] = rndst.binomial(1,p,n)

    ## add noise, with extreme noise for marked outliers
    df['y'] += ((1-df['outlier']) * df['latent_error'])
    df['y'] += (df['outlier'] * df['outlier_error'])

    ## round
    for col in ['y','latent_error','outlier_error','x']:
        df[col] = np.round(df[col],3)

    ## add label
    df['source'] = 'linear' if c == 0 else 'quadratic'

    ## create simple linspace for plotting true model
    plotx = np.linspace(df['x'].min() - np.ptp(df['x'])*.1
                        ,df['x'].max() + np.ptp(df['x'])*.1, 100)
    ploty = a + b*plotx + c*plotx**2
    dfp = pd.DataFrame({'x':plotx, 'y':ploty})

    return df, dfp


def interact_dataset(n=20, p=0, a=-30, b=5, c=0, latent_sigma_y=20):
    '''
    Convenience function:
    Interactively generate dataset and plot
    '''

    df, dfp = generate_data(n, p, a, b, c, latent_sigma_y)

    g = sns.FacetGrid(df, size=8, hue='outlier', hue_order=[True,False]
                      ,palette=sns.color_palette('Set1'), legend_out=False)

    _ = g.map(plt.errorbar, 'x', 'y', 'latent_error', marker="o"
              ,ms=10, mec='w', mew=2, ls='', elinewidth=0.7).add_legend()

    _ = plt.plot(dfp['x'], dfp['y'], '--', alpha=0.8)

    plt.subplots_adjust(top=0.92)
    _ = g.fig.suptitle('Sketch of Data Generation ({})'.format(df['source'][0])
                       ,fontsize=16)


def plot_datasets(df_lin, df_quad, dfp_lin, dfp_quad):
    '''
    Convenience function:
    Plot the two generated datasets in facets with generative model
    '''

    df = pd.concat((df_lin, df_quad), axis=0)
    dfp_lin, dfp_quad

    g = sns.FacetGrid(col='source', hue='source', data=df, size=6
                      ,sharey=False, legend_out=False)

    _ = g.map(plt.scatter, 'x', 'y', alpha=0.7, s=100, lw=2, edgecolor='w')

    _ = g.axes[0][0].plot(dfp_lin['x'], dfp_lin['y'], '--', alpha=0.6)
    _ = g.axes[0][1].plot(dfp_quad['x'], dfp_quad['y'], '--', alpha=0.6)


def plot_traces(traces, retain=1000):
    '''
    Convenience function:
    Plot traces with overlaid means and values
    '''

    ax = pm.traceplot(traces[-retain:], figsize=(12,len(traces.varnames)*1.5),
                      lines={k: v['mean'] for k, v in pm.df_summary(traces[-retain:]).iterrows()})

    for i, mn in enumerate(pm.df_summary(traces[-retain:])['mean']):
        ax[i,0].annotate('{:.2f}'.format(mn), xy=(mn,0), xycoords='data'
                         ,xytext=(5,10), textcoords='offset points', rotation=90
                         ,va='bottom', fontsize='large', color='#AA0022')


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

        with pm.Model() as models[nm]:

            print('\nRunning: {}'.format(nm))
            pm.glm.glm(fml, df, family=pm.glm.families.Normal())

            # For speed, we're using Metropolis here
            traces[nm] = pm.sample(5000, pm.Metropolis())[1000::5]

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
def run():
    n = 12
    import sys,os,inspect
    from astropy.io import ascii
    modelList=['k1','k2']#,'k3','k4']
    t='Ic'
    res=dict([])
    for ir in ['uv','J','H','K']:
        table=ascii.read(os.path.join('modjaz','type'+t,'tables',ir+'.dat'))
        allData=ascii.read(os.path.join('modjaz','type'+t,'tables',ir+'all.dat'))
        temp=pd.DataFrame({'x':np.array(table['time']),'y':np.array(table['mag']),'error':np.array(table['magerr'])})
        alltemp=pd.DataFrame({'x':np.array(allData['time']),'y':np.array(allData['mag'])*(-1),'error':np.array(allData['magerr'])})
        temp_xlims = (temp['x'].min() - np.ptp(temp['x'])/10,temp['x'].max() + np.ptp(temp['x'])/10)
        models_lin,traces_lin=run_models(temp,2)
        dfdic = pd.DataFrame(index=modelList, columns=['dic','waic'])
        dfdic.index.name = 'model'

        for nm in dfdic.index:
            dfdic.loc[nm, 'dic'] = pm.stats.dic(traces_lin[nm], models_lin[nm])
            dfdic.loc[nm, 'waic'] = pm.stats.waic(traces_lin[nm], models_lin[nm])[0]

        dfdic = pd.melt(dfdic.reset_index(), id_vars=['model'], var_name='poly', value_name='Information Criterion')

        #g = sns.factorplot(x='model', y='Information Criterion', col='poly', hue='poly', data=dfdic, kind='bar', size=6)
        #plt.show()
        dfwaic = pd.DataFrame(index=modelList, columns=['lin'])
        dfwaic.index.name = 'model'

        for nm in dfwaic.index:
            dfwaic.loc[nm, 'lin'] = pm.stats.waic(traces_lin[nm], models_lin[nm])[0]

        best=dfwaic[dfwaic['lin']==np.min(dfwaic['lin'])].index[0]
        dfwaic = pd.melt(dfwaic.reset_index(), id_vars=['model'], var_name=ir.upper(), value_name='waic')
        g = sns.factorplot(x='model', y='waic', col=ir.upper(), hue=ir.upper(), data=dfwaic, kind='bar', size=6)
        #plt.savefig(os.path.join('type'+t,'plots',ir.upper()+'_waic.pdf'),fmt='pdf')

        res[ir]=plot_posterior_cr(models_lin,traces_lin,temp,temp_xlims,datamodelnm=ir, bestModel=best,modelnms=modelList,allDat=alltemp,typ=t,bic=dfwaic)

    return(res)

















