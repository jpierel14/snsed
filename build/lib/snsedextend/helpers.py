import os,sncosmo
from .utils import _props,_filters
import numpy as np

__all__=['_findFile','_get_default_prop_name','_bandCheck','bandRegister','_standardize']

def _findFile(filename):
    """Helper function that does a quick recurisive directory seach for a transmission file.
    :param filename: Transmission file to find
    """
    for root, subFolders, files in os.walk('.'):
        if filename in files:
            return(os.path.abspath(os.path.join(root,filename)))

def _get_default_prop_name(prop):
    for key,value in _props.items():
        if {prop.lower()} & value:
            return key
    return prop



def _bandCheck(bandDict,band):
    """
    (Private)
     Helper function to check a band dictionary to make sure it contains all sncosmo.Bandpass objects. If it cannot find
     the band in the sncosmo registry, it will attempt to find a file with the name of that band.
    """
    try:
        bandDict[band]=sncosmo.get_bandpass(bandDict[band])
    except:
        try:
            wave,trans=np.loadtxt(_findFile(bandDict[band]),unpack=True)
            bandDict[band]=bandRegister(wave,trans,band,os.path.basename(bandDict[band]))
        except:
            raise RuntimeError('Band "%s" listed in bandDict but not in sncosmo registry and file not found.'%band)
    return(bandDict)

def bandRegister(wave,trans,band,bandName):
    sncosmo.registry.register(sncosmo.Bandpass(wave,trans,name=bandName))
    _filters[bandName]=sncosmo.get_bandpass(bandName)

def _standardize(table):
    for col in table.colnames:
        if '-' not in col and '_' not in col:
            if col != _get_default_prop_name(col.lower()):
                table.rename_column(col, _get_default_prop_name(col.lower()))
    return(table)