import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/SO/libs/')
import matplotlib as mpl
# mpl.use('agg')
import os
import numpy as np
import pickle
import xarray as xr
from myTrend import myfitting2d_sttcs,myfitting1d_sttcs
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata 
import warnings
warnings.filterwarnings('ignore')


pthMO='J:/MDLS_OBS_OHC/Area_data/'
wpth='C:/Users/shjo/OneDrive/mySO/zonal_sum/tmp/'

t_rng=[1993, 2020]
varnm='OHC7002000'
fig_bool=1

myName='Zonally integrated '+varnm+' trend'


### Preparation ============================================================
wpth=wpth+varnm+'_'+str(t_rng[0])+'_'+str(t_rng[-1])+'/'
try :
    os.mkdir(wpth)
except:
    raise

myMDOB=[pthMO+i for i in os.listdir(pthMO) if i.endswith('.nc')]

myDATA=myMDOB
if t_rng[0]<1992:
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myECCO')]
print(myDATA)

### Read myDATA =============================================================


print('!!! Open files !!!')
myCoefs,myNm,myLat=[],[],[]
lat_rng=[-80,-30]; lon_rng=[0,360]

time_rng=[str(t_rng[0])+'-01',str(t_rng[-1])+'-12']


for i in myDATA: 
            
    print('!!! Open: '+i+' !!!')
    tmp=xr.open_dataset(i)
    Area=tmp['Area'].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]))]
    mySST = tmp[varnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1])\
        ,time=slice(time_rng[0],time_rng[-1]))]
    mySST=mySST.where(mySST<10**30)
    # mySST=mySST.sum(dim='lon',skipna=True)
    time,latR,lonR=mySST.time.values,mySST.lat.values,mySST.lon.values
    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+varnm+' trend '+\
        str(time[0])[:4]+' '+str(time[-1])[:4]
    ### SST Coef =============================================================
    print('!!! linregress !!!')
    slope,intercept,r_value,p_value,std_err,smask=myfitting2d_sttcs(mySST,threshold=0.05)
    CoefD=slope*12*10**-8*Area # Decadal^-1
    
    CoefD_P=CoefD.loc[dict(lon=slice(153,290))]
    CoefD_I=CoefD.loc[dict(lon=slice(20,117))]
    CoefD_S=CoefD.loc[dict(lon=slice(180,240),lat=slice(-60,-50))]
    CoefD_A=xr.concat([CoefD.loc[dict(lon=slice(-70+360,360))],CoefD.loc[dict(lon=slice(0,20))]],dim='lon')

    ### Figure configs =======================================================
    CoefD_posi=np.nansum(CoefD.where(CoefD>0))
    CoefD_nega=np.nansum(CoefD.where(CoefD<0))

    CoefD_P_posi=np.nansum(CoefD_P.where(CoefD_P>0))
    CoefD_P_nega=np.nansum(CoefD_P.where(CoefD_P<0))

    CoefD_I_posi=np.nansum(CoefD_I.where(CoefD_I>0))
    CoefD_I_nega=np.nansum(CoefD_I.where(CoefD_I<0))
    
    CoefD_A_posi=np.nansum(CoefD_A.where(CoefD_A>0))
    CoefD_A_nega=np.nansum(CoefD_A.where(CoefD_A<0))

    CoefD_S_posi=np.nansum(CoefD_S.where(CoefD_S>0))
    CoefD_S_nega=np.nansum(CoefD_S.where(CoefD_S<0))
    
    print('!!! ========== !!!')
    print(i.split('/')[-1])
    print('SCP: ',np.nansum(CoefD_S_nega)/CoefD_nega*100,'%')
    print('Pacific: ',np.nansum(CoefD_P_nega)/CoefD_nega*100,'%')
    print('Atlantic: ',np.nansum(CoefD_A_nega)/CoefD_nega*100,'%')
    print('Indian: ',np.nansum(CoefD_I_nega)/CoefD_nega*100,'%')
    print('!!! ========== !!!')
