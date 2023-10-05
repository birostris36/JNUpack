import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/SO/libs/')
import matplotlib as mpl
# mpl.use('agg')
import os
import numpy as np
import pickle
import pandas as pd
import xarray as xr
from myTrend import myfitting2d_sttcs,myfitting1d_sttcs
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata 
import warnings
warnings.filterwarnings('ignore')

pthMO='J:/MDLS_OBS_OHC/Area_data/'
wpth='C:/Users/shjo/OneDrive/mySO/SCP_signal/'

t_rng=[1993, 2017]
varnm='OHC700'
fig_bool=1

myName='SCP areal mean'

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
lat_rng=[-60,-50]; lon_rng=[180,240]

time_rng=[str(t_rng[0])+'-01',str(t_rng[-1])+'-12']

mySig=[]
for i in myDATA: 
            
    print('!!! Open: '+i+' !!!')
    tmp=xr.open_dataset(i)
    Area=tmp['Area'].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]))]
    mySST = tmp[varnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1])\
        ,time=slice(time_rng[0],time_rng[-1]))]
    mySST=mySST.where(mySST<10**30)
    mySST=mySST.mean(dim=['lat','lon'])
    mySig.append(mySST)
time=mySST.time.values

xtick_location = time[::12*5]
xtick_labels = time[::12*5]

### Figure mySig =============================================================
Label_size = 18
fig, axs = plt.subplots(1,1,figsize=(10,4.5),constrained_layout = True,dpi=200)
for i in mySig:
    axs.plot(time,i,linewidth=1,color='grey',zorder=0,alpha=0.7)
axs.plot(time,np.nanmean(np.array(mySig)),label='Ensembled',linewidth=2,color='r',zorder=0)
axs.axhline(y=-55,color='k',linestyle='-.')
axs.axhline(y=-60,color='k',linestyle='-.')
axs.axvline(x=0,color='k',linestyle='-')
axs.set_xticks(ticks=xtick_location)
axs.set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
axs.set_title(myName+'\nIndian '+str(t_rng[0])+' '+str(t_rng[-1]),loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'normal'})
axs.tick_params(axis='both', labelsize=Label_size)
axs.grid(axis='x',linestyle='-.')
# axs.set_xlim([-0.35,0.35])
axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size-3, labelcolor='k', top=True,width=1.)
axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
plt.legend(fontsize='13')
plt.tight_layout()
if fig_bool:
    myName.replace(' ','_')
    plt.savefig(wpth+'/'+myName.replace(' ','_')+'_Indian'+'_ensembled_'+str(t_rng[0])+'_'+str(t_rng[-1]),bbox_inches='tight')
plt.show()






