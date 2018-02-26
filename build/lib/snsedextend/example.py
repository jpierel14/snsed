import sncosmo,os
from .helpers import *

__all__=['load_example_data','example_sed']

example_sed=os.path.join('example_setup','SDSS-004012.SED')

def load_example_data():
	return(_standardize(sncosmo.read_lc(os.path.join('example_setup','example.lc'))))
