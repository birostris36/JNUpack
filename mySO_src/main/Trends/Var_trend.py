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
from myTrend import myfitting2d_sttcs,myRegress3d_sttcs,myXfitting2d
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib as mpl
mpl.use('agg')
### ======================================================================
pthMO='J:/tmp_proc/Models/'
pthrn='J:/tmp_proc/Obs/'

varnm='temp'

t_rng=['2011-01', '2022-12']
d_rng=[700,2000]
lat_rng=[-75,-30]; lon_rng=[0,360]

wpth='C:/Users/shjo/OneDrive/mySO/Trend/Depth_avg/'


# print(mySAM)

### Preparation ============================================================
wpth=wpth+varnm+'_trend_'+str(t_rng[0][:4])+'_'+str(t_rng[-1][:4])+'_'+\
    str(lat_rng[0])+'S'+str(lat_rng[-1])+'S'+'_'+\
        'E'+str(lon_rng[-1])+'E_'+str(d_rng[0])+'m'+str(d_rng[-1])+'m/'
                
wpth=wpth.replace('-','')

try :
    os.mkdir(wpth)
except:
    pass
loc=sys._getframe().f_code.co_filename
myInfo(loc,wpth)

myRnly=[pthrn+i for i in os.listdir(pthrn) if i.endswith('.nc')]
myMDOB=[pthMO+i for i in os.listdir(pthMO) if i.endswith('.nc')]                     
myDATA=myMDOB+myRnly
if int(t_rng[0].split('-')[0])<1992:
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myECCO')]
myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myISA')]

### Read Data ==============================================================
for i in myDATA: 
    print('!!! Open: '+i+' !!!')
    tmp=xr.open_dataset(i)
    if len(tmp.coords)==4:
        mydata = tmp[varnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1])\
            ,time=slice(t_rng[0],t_rng[-1]),depth=slice(d_rng[0],d_rng[-1]))].mean(dim='depth')
    
    mydata=mydata.where(mydata<1000)
    mydata=mydata.where(mydata>-1000)
    
    mydata=mydata.fillna(0)
    mydata=mydata.rolling(time=12,center=True).mean()[6:-5]
    
    lonR,latR=mydata.lon.values,mydata.lat.values
    lonR_m,latR_m=np.meshgrid(lonR,latR)
    time=mydata.time.values
    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+varnm+' '+\
        str(d_rng[0])+'~'+str(d_rng[-1])+'m'+' trend\n'+\
        str(lon_rng[0])+'~'+str(lon_rng[-1])+'E '+\
            str(time[0])[:4]+' '+str(time[-1])[:4]
                
    dta_snm=i.split('/')[-1][2:-3].split('_')[0]+' '+varnm+'_trend'+\
        str(lat_rng[0])+'S'+str(lat_rng[-1])+'S'+' '+str(lon_rng[0])+'E'+str(lon_rng[-1])+'E_'+\
            str(time[0])[:4]+' '+str(time[-1])[:4]+str(d_rng[0])+'~'+str(d_rng[-1])+'m'
    
    dta_snm=dta_snm.replace(' ','_').replace('salt','salinity').replace('-','')
    dta_nm=dta_nm.replace('salt','salinity').replace('-','')
    
    ### SST Coef =============================================================
    print('!!! linregress !!!')
    # plt.figure()
    # plt.scatter(range(300),mydata[:,30,20])
    # plt.show()
    # raise
    # slope,intercept,r_value,p_value,std_err,smask=myRegress3d_sttcs(myIdx,mydata,threshold=0.05)
    # CoefD=slope*10**10 # Decadal^-1
    
    _,slope=myXfitting2d(mydata,'time')
    
    CoefD=slope.values*12*10 # Decadal^-1
    
    # smask[CoefD==0]=np.nan
    ### Figure configs =======================================================
    # myCoefs.append(CoefD); myNm.append(dta_nm); myLat.append(latR)
    # raise
    
    myN=16
    mylim=[-.5,.5]
    CMAP,mylevel=myClrbr('myblc2',mylim,myN)
    CMAP_salt,mylevel_salt=myClrbr('salt',mylim,myN)
    CMAP_temp,mylevel_temp=myClrbr('balance',mylim,myN)

    CoefD[CoefD<mylim[0]]=mylim[0]
    CoefD[CoefD>mylim[-1]]=mylim[-1]
    
    mySetting={
        'figsize': '',
        'mylabel': '',
        'Label_size':18,
        'title_loc':'right',
        'fontParams':'Arial',
        'wpth':wpth}
    
    lat_rng_,lon_rng_=[-60,-53],[200,250]
    lat_rng_,lon_rng_=[10,11],[200,250]

    smask=np.zeros_like(CoefD); CoefD[0]=-1
    F=figmaster(mySetting)
    print(dta_snm)
    F.myCrtpy_sph4_X(latR_m,lonR_m,CoefD,CMAP,mylevel,dta_nm,dta_snm,lat_rng_,lon_rng_)

