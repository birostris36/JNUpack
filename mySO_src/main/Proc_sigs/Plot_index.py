import matplotlib.pyplot as plt
import xarray as xr
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import matplotlib as mpl
from scipy.interpolate import griddata 
import warnings
import os
import numpy as np
sys.path.append('C:/Users/shjo/Bridge/JNUpack/mySO_src/libs/')
warnings.filterwarnings('ignore')
from myTools import myInfo
from myTrend import myfitting2d_sttcs,myRegress3d_sttcs
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib as mpl
mpl.use('agg')
### ======================================================================
npth='J:/Reanalysis/myOISST_198001_202212_sst.nc'


myindx='SAM'
varnm='ice'

t_rng=['1993-01', '2017-12']
lat_rng=[-75,-30]; lon_rng=[0,360]

wpth='C:/Users/shjo/OneDrive/mySO/Regression_H/'+myindx+'/'
mySig_npth='C:/Users/shjo/OneDrive/mySO/mySignals/my'+myindx+'.pkl'


### load Index ==================================================================
with open(mySig_npth, 'rb') as f:
    data = pickle.load(f)
myIdx=data.loc['1993-01':'2017-12']
myIdx=myIdx.rolling(window=12,center=True).mean()[6:-5].values.reshape(-1)


















