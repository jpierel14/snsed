from collections import OrderedDict as odict
import os

__dir__=os.path.abspath(os.path.dirname(__file__))


__all__=['_filters','_opticalBands','_props','_fittingParams','_needs_bounds','_UVrightBound','_UVleftBound','_IRleftBound','_IRrightBound','_zp','__dir__']


#default transmission files, user can define their own
_filters={
    'U':'bessellux',
    'B':'bessellb',
    'V':'bessellv',
    'R':'bessellr',
    'I':'besselli',
    'J':'paritel::j',
    'H':'paritel::h',
    'K':'paritel::ks',
    'Ks':'paritel::ks',
    'u':'sdss::u',
    'g':'sdss::g',
    'r':'sdss::r',
    'i':'sdss::i',
    'Z':'sdss::z'
}

_opticalBands=['B','V','R','g','r']

_props=odict([
    ('time',{'mjd', 'mjdobs', 'jd', 'time', 'date', 'mjd_obs','mhjd','jds','mjds'}),
    ('band',{'filter', 'band', 'flt', 'bandpass'}),
    ('flux',{'flux', 'f','fluxes'}),
    ('fluxerr',{'flux_error', 'fluxerr', 'fluxerror', 'fe', 'flux_err','fluxerrs'}),
    ('zp',{'zero_point','zp', 'zpt', 'zeropoint'}),
    ('zpsys',{'zpsys', 'magsys', 'zpmagsys'}),
    ('mag',{'mag','magnitude','mags'}),
    ('magerr',{'magerr','magerror','magnitudeerror','magnitudeerr','magerrs','dmag'})
])

_fittingParams=odict([
    ('curve',None),
    ('method','minuit'),
    ('params',None),
    ('bounds',None),
    ('ignore',None),
    ('constants',None),
    ('dust',None),
    ('effect_names',None),
    ('effect_frames',None)
])

#dictionary of zero-points (calculated later based on transmission files)
_zp={}

_needs_bounds={'z'}


_UVrightBound=4000
_UVleftBound=1200
_IRleftBound=9000
_IRrightBound=55000