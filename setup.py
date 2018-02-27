from setuptools import setup
import os,glob,warnings,sys,fnmatch

def recursive_glob(basedir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(basedir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

PACKAGENAME='snsedextend'
# Add the project-global data
pkgdatadir = os.path.join(PACKAGENAME, 'data')
exdatadir = os.path.join(PACKAGENAME, 'example_setup')
data_files = []
data_files.extend(recursive_glob(pkgdatadir, '*'))
data_files.extend(recursive_glob(exdatadir, '*'))
data_files = [f[len(PACKAGENAME)+1:] for f in data_files]

#warnings.warn("I'm about to install sncosmo from github, the development version",RuntimeWarning)
#answer=raw_input('Proceed (y/n)? ')
#if answer.lower() != 'y':
#    print('Exiting...')
#    sys.exit()
os.system('pip install git+https://github.com/sncosmo/sncosmo')
setup(
	name='snsedextend',
	install_requires=['cython','numpy','scipy',
                          'astropy','matplotlib','pandas','pymc3>=3.0'],
        packages=['snsedextend'],
	author='J. R. Pierel',
	author_email='jr23@email.sc.edu',
	version='0.2.0',
        package_data={'snsedextend': data_files}
)
