
import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/SO/libs/')
import matplotlib as mpl
mpl.use('agg')
import pickle
import pandas as pd
import numpy as np
import os
import numpy as np
import xarray as xr
from myTrend import myRegress3d_sttcs
from myPlot import  figmaster,myClrbr
from myTools import myInfo
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Index path
spth='C:/Users/shjo/OneDrive/mySO/mySignals/SAM_index.pickle'
#Data Path
pthMO='J:/MDLS_OBS_OHC/'

wpth='C:/Users/shjo/OneDrive/mySO/Regress/SAM_OHC/'
t_rng=[1980, 2017]
vrnm='OHC700'

### Preparation ============================================================
wpth=wpth+vrnm+str(t_rng[0])+'_'+str(t_rng[-1])+'/'
try :
    os.mkdir(wpth)
    loc=sys._getframe().f_code.co_filename
    myInfo(loc,wpth)
except:
    raise

myMDOB=[pthMO+i for i in os.listdir(pthMO) if i.endswith('.nc')]

myDATA=myMDOB

if t_rng[0]<1992:
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myECCO')]

### Read index =============================================================
with open(spth, 'rb') as f:
    data = pickle.load(f)
SAM_9320_=data['1970-01':'2020-12'].rolling(window=12,center=True).mean()
SAM_9320=SAM_9320_['1980-01':'2017-12']
# SAM_9320=data['1980-01':'2017-12']
### Read myDATA =============================================================
print('!!! Open files !!!')
lat_rng=[-80,-30] ; time_rng=[str(t_rng[0])+'-01',str(t_rng[-1])+'-12']
for i in myDATA: 
    print('!!! Open: '+i+' !!!')
    tmp=xr.open_dataset(i)

    mySST = tmp[vrnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1])\
        ,time=slice(time_rng[0],time_rng[-1]))]

    mySST=mySST.where(mySST<10**30)
    lonR,latR=mySST.lon.values,mySST.lat.values
    lonR_m,latR_m=np.meshgrid(lonR,latR)
    time=mySST.time.values
    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+vrnm+' Regression (SAM) '+\
        str(time[0])[:4]+' '+str(time[-1])[:4]

    ### SST Coef =============================================================
    print('!!! linregress !!!')
    slope,intercept,r_value,p_value,std_err,smask=myRegress3d_sttcs(SAM_9320.values.reshape(-1),mySST,threshold=0.05)
    CoefD=slope*10**-9 # Decadal^-1
    ### Figure configs =======================================================
    myN=16
    sstTlim=[-1.,1.]
    CMAP,mylevel=myClrbr('myblc2',sstTlim,myN)

    CoefD[CoefD<sstTlim[0]]=sstTlim[0]
    CoefD[CoefD>sstTlim[-1]]=sstTlim[-1]

    mySetting={'title_loc':'right',
            'wpth':wpth}

    mySetting={
        'figsize': '',
        'mylabel': '',
        'Label_size':18,
        'title_loc':'right',
        'fontParams':'Arial',
        'wpth':wpth}

    ### Figure ==============================================================
    print('!!!   figure   !!!')
    F=figmaster(mySetting)
    F.myCrtpy_sph2(latR_m,lonR_m,CoefD,smask,CMAP,mylevel,dta_nm.replace('(','').replace(')',''))
    
