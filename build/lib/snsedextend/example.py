import sncosmo,os
from .helpers import *
dir=os.path.abspath(os.path.dirname(__file__))
__all__=['load_example_data','example_sed']

example_sed=os.path.join(dir,'example_setup','SEDs','SDSS-004012.SED')

def load_example_data():
	return(_standardize(sncosmo.read_lc(os.path.join(dir,'example_setup','example.lc'))))
