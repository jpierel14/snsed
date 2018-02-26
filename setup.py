from setuptools import setup
import os

os.system('pip install git+https://github.com/sncosmo/sncosmo')
setup(
	name='snsedextend',
	install_requires=['astropy'],
	author='J. R. Pierel',
	author_email='jr23@email.sc.edu',
	version='0.0.1'
)