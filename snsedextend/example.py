import sncosmo,os,glob
from .helpers import *
from .utils import __dir__
__all__=['load_example_lc','load_example_lcs','example_sed','example_seds']

example_sed=[os.path.join(__dir__,'example_setup','SEDs','SDSS-004012.SED')]
example_seds=glob.glob(os.path.join(__dir__,'example_setup','SEDs','*.SED'))

def load_example_lc():
	return(_standardize(sncosmo.read_lc(os.path.join(__dir__,'example_setup','lcs','example_2006aj.lc'))))

def load_example_lcs():
	return([_standardize(sncosmo.read_lc(f))for f in glob.glob(os.path.join(__dir__,'example_setup','lcs','*.lc'))])


