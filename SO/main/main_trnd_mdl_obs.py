import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/SO/libs/')
import matplotlib as mpl
mpl.use('agg')
import os
import numpy as np
import xarray as xr
from myTrend import myfitting2d_sttcs
from myPlot import  figmaster,myClrbr
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

myttl='trend 1980 2020'

# pthrn='J:/Reanalysis/'
pthmd='J:/Models/'
pthob='J:/Obs/'

wpth='J:/mytest/'

# myRnly=[pthrn+i for i in os.listdir(pthrn) if i.endswith('.nc')]
myMdls=[pthmd+i for i in os.listdir(pthmd) if i.endswith('.nc')]
myObsv=[pthob+i for i in os.listdir(pthob) if i.endswith('.nc')]

myDTA=myMdls+myObsv

### Reanalysis SST trend =================================================
print('!!! Open files !!!')
lat_rng=[-80,-30] ; time_rng=['1982-01','2020-12']
for i in myDTA: 
    print('!!! Open: '+i+' !!!')
    mySST = xr.open_dataset(i).temp.loc[dict(depth=slice(0,10),lat=slice(lat_rng[0],lat_rng[-1])\
        ,time=slice(time_rng[0],time_rng[-1]))].mean(dim='depth')
    mySST=mySST.where(mySST<1000)
    lonR,latR=mySST.lon.values,mySST.lat.values
    lonR_m,latR_m=np.meshgrid(lonR,latR)

    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+myttl

    ### SST Coef =============================================================
    print('!!! linregress !!!')
    slope,intercept,r_value,p_value,std_err,smask=myfitting2d_sttcs(mySST,threshold=0.05)
    CoefD=slope*12 # year^-1

    ### Figure configs =======================================================
    myN=16
    sstTlim=[-0.03,0.03]
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
    print('!!! figure !!!')
    F=figmaster(mySetting)
    F.myCrtpy_sph2(latR_m,lonR_m,CoefD,smask,CMAP,mylevel,dta_nm)
