import sncosmo,os
from .helpers import *
from .utils import __dir__
__all__=['load_example_data','example_sed']

example_sed=os.path.join(__dir__,'example_setup','SEDs','SDSS-004012.SED')

def load_example_data():
	return(_standardize(sncosmo.read_lc(os.path.join(__dir__,'example_setup','example.lc'))))
