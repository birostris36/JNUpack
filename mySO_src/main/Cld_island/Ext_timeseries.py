import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/mySO_src/libs/')
import matplotlib as mpl
# mpl.use('agg')
import os
import numpy as np
import xarray as xr
from myTrend import myfitting2d_sttcs
from myPlot import  figmaster,myClrbr,dta_colr
from myTools import myInfo
from scipy import io
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

pthmd='J:/tmp_proc/Models/'
pthob='J:/tmp_proc/Obs/'

wpth='C:/Users/shjo/OneDrive/mySO/Cld_island/TimeSeries/'

lat_rng=[-60,-53]; lon_rng=[200,220]
t_rng=[1960, 2020]; d_rng=[0,100]
varnm='temp'

### Preparation ============================================================
time_rng=[str(t_rng[0])+'-01',str(t_rng[-1])+'-12']
tmp_sv_nm=str(lon_rng[0])+'E'+str(lon_rng[-1])+'E'+'_'+str(lat_rng[0])+'S'+str(lat_rng[-1])+'S'
tmp_sv_nm=tmp_sv_nm.replace('-','')
# wpth=wpth+varnm+'_'+tmp_sv_nm+'/'

my4DNM=''

wpth=wpth+varnm+'_'+str(t_rng[0])+'_'+str(t_rng[-1])+'_'+tmp_sv_nm+'/'
try :
    os.mkdir(wpth)
except:
    pass
loc=sys._getframe().f_code.co_filename
myInfo(loc,wpth)

# myRnly=[pthrn+i for i in os.listdir(pthrn) if i.endswith('.nc')]
myMdls=[pthmd+i for i in os.listdir(pthmd) if i.endswith('.nc')]
myObsv=[pthob+i for i in os.listdir(pthob) if i.endswith('.nc')]

myDATA=myMdls+myObsv

# if t_rng[0]<1992:
#     myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myECCO')]
#     myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myISAS')]
#     myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myARMOR')]

plt.rcParams["font.family"] = 'Arial'

### Read myDATA =============================================================
# str(lon_rng[0])+'E'+str(lon_rng[-1])+'E'+'_'+str(lat_rng[0])+'S'+str(lat_rng[-1])+'S'
print('!!! Open files !!!')
myEofs,myNm,myLat,myLon=[],[],[],[]
myPcs,myVar,myVar2=[],[],[]
for i in myDATA: 
    print('!!! Open: '+i+' !!!')
     
    tmp=xr.open_dataset(i)
    
    if len(tmp.coords)==4:
        mydata = tmp[varnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]),\
            time=slice(time_rng[0],time_rng[-1]),depth=slice(d_rng[0],d_rng[-1]))].mean(dim='depth')
        my4DNM='_'+str(d_rng[0])+'m'+str(d_rng[-1])+'m'
    elif len(tmp.coords)==3:
        mydata = tmp[varnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]),\
            time=slice(time_rng[0],time_rng[-1]))]

    mydata=mydata.where(mydata>-10**26)
    mydata=mydata.where(mydata<10**26)
    
    mydata=mydata.mean(dim=['lat','lon'])
    
    time=mydata.time.values
    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+varnm+' '+\
        str(time[0])[:4]+' '+str(time[-1])[:4]+' '+tmp_sv_nm  
    
    myNM=i.split('/')[-1][2:-3].split('_')[0]
    
    mytt=pd.DatetimeIndex([str(tt)[:7] for tt in mydata.time.values])

    if i==myDATA[0]:
        myPD=pd.DataFrame({myNM:mydata.values},index=mytt)
    else:
        myPD=pd.concat([myPD, pd.DataFrame({myNM:mydata.values},index=mytt)] ,axis=1)
            # np.save(wpth+dta_nm.replace(' ','_'),mydata)
    # print(myPD)
    
with open(wpth+'/'+wpth.split('/')[-2]+my4DNM+'.pickle', 'wb') as f:
    pickle.dump(myPD, f, pickle.HIGHEST_PROTOCOL)



