import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/mySO_src/libs/')
import matplotlib as mpl
mpl.use('agg')
from myPlot import  figmaster,myClrbr
from myTools import myInfo
import matplotlib.path as mpath
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
from eofs.xarray import Eof
import numpy as np
from scipy import io
import xarray as xr
import pickle
from myTrend import myfitting2d_sttcs,myfitting1d_sttcs
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata 
import warnings
warnings.filterwarnings('ignore')

pthmd='J:/tmp_proc/Models/'
pthob='J:/tmp_proc/Obs/'

wpth='C:/Users/shjo/OneDrive/mySO/REOF/Verical_temp/'

lat_rng=[-75,-30]; lon_rng=[230,250]
t_rng=[1993, 2020]
d_rng=[0,2000]
vrnm='temp'
    
### Preparation ============================================================
time_rng=[str(t_rng[0])+'-01',str(t_rng[-1])+'-12']
tmp_sv_nm=str(lon_rng[0])+'E'+str(lon_rng[-1])+'E'+'_'+str(lat_rng[0])+'S'+str(lat_rng[-1])+'S'
tmp_sv_nm=tmp_sv_nm.replace('-','')
# wpth=wpth+varnm+'_'+tmp_sv_nm+'/'

wpth_re=wpth+vrnm+'_'+str(t_rng[0])+'_'+str(t_rng[-1])+'_'+tmp_sv_nm+'/'
try :
    os.mkdir(wpth_re)
except:
    pass
loc=sys._getframe().f_code.co_filename
myInfo(loc,wpth)
# myRnly=[pthrn+i for i in os.listdir(pthrn) if i.endswith('.nc')]
myMdls=[pthmd+i for i in os.listdir(pthmd) if i.endswith('.nc')]
myObsv=[pthob+i for i in os.listdir(pthob) if i.endswith('.nc')]

myDATA=myMdls+myObsv

if t_rng[0]<1992:
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myECCO')]
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myISAS')]
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myARMOR')]

plt.rcParams["font.family"] = 'Arial'

### Read myDATA =============================================================
# str(lon_rng[0])+'E'+str(lon_rng[-1])+'E'+'_'+str(lat_rng[0])+'S'+str(lat_rng[-1])+'S'

print('!!! Open files !!!')
myEofs,myNm,myLat,myLon=[],[],[],[]
myPcs,myVar,myVar2=[],[],[]
for i in myDATA: 
    print('!!! Open: '+i+' !!!')
    tmp=xr.open_dataset(i)
    
    mydata = tmp.loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1])\
        ,time=slice(time_rng[0],time_rng[-1]),depth=slice(d_rng[0],d_rng[-1]))]
    
    mydata=mydata.where(mydata>-10**26)
    mydata=mydata.where(mydata< 10**26)
    
    mydata=mydata.mean(dim='lon',skipna=True)

    myvrnm=mydata[vrnm]
    time,latR,depthR=mydata.time.values,mydata.lat.values,mydata.depth.values
    
    # print(myvrnm)
    myvrnm_1Y=myvrnm.rolling(time=12,center=True).mean()[6:-5]
    
    # mydata=mydata.values

    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+vrnm+' '+\
        str(time[0])[:4]+' '+str(time[-1])[:4]+' '+tmp_sv_nm
    mySgnl={'sgnl':myvrnm_1Y.values}
    io.savemat(wpth_re+dta_nm.replace(' ','_')+\
    '.mat', mySgnl)
    



