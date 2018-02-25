import warnings
warnings.filterwarnings('ignore',category=FutureWarning)
#warnings.filterwarnings('ignore',category=ShimWarning)
import sncosmo, definedBands,BIC,os
from .colorCalc import *


#Automatically sets environment to your SNDATA_ROOT file (Assumes you've set this environment variable,otherwise will use the builtin version)
try:
    sndataroot = os.environ['SNDATA_ROOT']
except:
    os.environ['SNDATA_ROOT']=os.path.join(os.path.dirname(snsedextend),'SNDATA_ROOT')
    sndataroot=os.environ['SNDATA_ROOT']
