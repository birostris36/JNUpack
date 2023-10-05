import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/SO/libs/')
import matplotlib as mpl
mpl.use('agg')
import os
import numpy as np
import xarray as xr
import pickle
from myTrend import myfitting2d_sttcs,myfitting1d_sttcs
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata 
import warnings
warnings.filterwarnings('ignore')

pthMO='J:/MDLS_OBS_OHC/'
pthrn='J:/Reanalysis/'

wpth='C:/Users/shjo/OneDrive/mySO/zonal_mean_Indian/'

t_rng=[2005, 2020]
varnm='sst'
fig_bool=1

myName='Zonal mean '+varnm+' trend'
    

### Preparation ============================================================
wpth=wpth+varnm+'_'+str(t_rng[0])+'_'+str(t_rng[-1])+'/'
try :
    os.mkdir(wpth)
except:
    raise

myRnly=[pthrn+i for i in os.listdir(pthrn) if i.endswith('.nc')]
myMDOB=[pthMO+i for i in os.listdir(pthMO) if i.endswith('.nc')]

myDATA=myRnly+myMDOB

if t_rng[0]<1992:
    myDATA=[i for i in myDATA if not i.split('/')[-1].startswith('myECCO')]

plt.rcParams["font.family"] = 'Arial'

### Read myDATA =============================================================

print('!!! Open files !!!')
myCoefs,myNm,myLat=[],[],[]
lat_rng=[-75,-30]; lon_rng=[20,117]; time_rng=[str(t_rng[0])+'-01',str(t_rng[-1])+'-12']
for i in myDATA: 
    print('!!! Open: '+i+' !!!')
    tmp=xr.open_dataset(i)
    if len(tmp.coords)==3:
        mySST = tmp[varnm].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]),\
            time=slice(time_rng[0],time_rng[-1]))]
    else:
        mySST = tmp['temp'].loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1])\
    ,time=slice(time_rng[0],time_rng[-1]),depth=slice(0,10))].mean(dim='depth')
        
    mySST=mySST.where(mySST<10**30)
    mySST=mySST.mean(dim='lon',skipna=True)
    
    time,latR=mySST.time.values,mySST.lat.values
    dta_nm=i.split('/')[-1][2:-3].split('_')[0]+' '+varnm+' trend '+\
        str(time[0])[:4]+' '+str(time[-1])[:4]

    ### SST Coef =============================================================
    print('!!! linregress !!!')
    slope,intercept,r_value,p_value,std_err,smask=myfitting1d_sttcs(mySST,threshold=0.05)
    CoefD=slope*12*10 # Decadal^-1
    ### Figure configs =======================================================
    myCoefs.append(CoefD); myNm.append(dta_nm); myLat.append(latR)


myCoefs_interp=myCoefs[0].reshape([1,len(myCoefs[0])])
for i,j in zip(myCoefs[1:],myLat[1:]):
    myCoefs_interp=np.concatenate([myCoefs_interp,griddata(j,i,myLat[0]).reshape([1,len(myCoefs[0])])],axis=0)
Ensemble=np.nanmean(myCoefs_interp,axis=0)

### Save myTrend ====================================================================
mydata={}
for i,j,nm in zip(myCoefs_interp,myLat,myNm):
    mydata[nm.replace(' ','_')]=i
    mydata[nm.replace(' ','_')+'_lat']=j
mydata['Ensemble']=Ensemble
mydata['Ensemble_lat']=myLat[0]

with open(wpth+myName.replace(' ','_')+'_'+str(t_rng[0])+'_'+str(t_rng[-1])+'.pickle', 'wb') as f:
    pickle.dump(mydata, f, pickle.HIGHEST_PROTOCOL)

### Figure Ensembled ==============================================================
Label_size = 18
fig, axs = plt.subplots(1,1,figsize=(4.5,6.5),constrained_layout = True,dpi=200)
for i,j,k in zip(myCoefs,myNm,myLat):
    axs.plot(i,k,linewidth=1,color='grey',zorder=0,alpha=0.7)
axs.plot(Ensemble,myLat[0],label='Ensembled',linewidth=2,color='r',zorder=0)
axs.axhline(y=-55,color='k',linestyle='-.')
axs.axhline(y=-60,color='k',linestyle='-.')
axs.axvline(x=0,color='k',linestyle='-')
axs.set_title(myName+'\nIndian '+str(t_rng[0])+' '+str(t_rng[-1]),loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'normal'})
axs.tick_params(axis='both', labelsize=Label_size)
axs.grid(axis='x',linestyle='-.')
axs.set_xlim([-0.35,0.35])
axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size-3, labelcolor='k', top=True,width=1.)
axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
plt.legend(fontsize='13')
plt.tight_layout()
if fig_bool:
    myName.replace(' ','_')
    plt.savefig(wpth+'/'+myName.replace(' ','_')+'_Indian'+'_ensembled_'+str(t_rng[0])+'_'+str(t_rng[-1]),bbox_inches='tight')

plt.show()

myNm=[i.split(' ')[0] for i in myNm]
### Figure ==============================================================
Label_size = 18
fig, axs = plt.subplots(1,1,figsize=(4.5,6.5),constrained_layout = True,dpi=200)
for i,j,k in zip(myCoefs,myNm,myLat):
    axs.plot(i,k,label=j,linewidth=2,zorder=0,alpha=1.,color=dta_colr(j))
axs.axhline(y=-55,color='k',linestyle='-.')
axs.axhline(y=-60,color='k',linestyle='-.')
axs.axvline(x=0,color='k',linestyle='-')
axs.set_title(myName+'\nIndian '+str(t_rng[0])+' '+str(t_rng[-1]),loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'normal'})
axs.tick_params(axis='both', labelsize=Label_size)
axs.grid(axis='x',linestyle='-.')
axs.set_xlim([-0.35,0.35])
axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size-3, labelcolor='k', top=True,width=1.)
axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
plt.legend(fontsize="11",loc='lower right')
plt.tight_layout()
if fig_bool:
    myName.replace(' ','_')
    plt.savefig(wpth+'/'+myName.replace(' ','_')+'_Indian_'+str(t_rng[0])+'_'+str(t_rng[-1]),bbox_inches='tight')
plt.show()
















