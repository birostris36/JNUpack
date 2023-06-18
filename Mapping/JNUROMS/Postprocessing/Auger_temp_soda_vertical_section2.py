# -*- coding: utf-8 -*-
"""
Created on Tue May  2 06:22:17 2023

@author: shjo9
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May  2 06:03:35 2023

@author: shjo9
"""


import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import xarray as xr
import dask
import matplotlib as mpl
import cmocean
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp2d, griddata
import os
import datetime as dt
from matplotlib import colors,colorbar
from tqdm import tqdm
import pandas as pd

# =============================================================================
#  PLOT 1
# =============================================================================
data_path='G:/SODA/'
t_rng=['1980-01','2017-12']
w_path_sig='D:/OneDrive/base142/Model_results/TK0525ED_CLM/vertical_temp_60_soda/'
save_name='Vertical_SODA_temp_lon60'
Tcline=300
fig_bool=0

DATA_path=[data_path+i for i in os.listdir(data_path) if i.endswith('.nc')]

Sample_data=Dataset(DATA_path[0])
LAT,LON=Sample_data['yt_ocean'][:],Sample_data['xt_ocean'][:]
DEPTH=Sample_data['st_ocean'][:]
Sample_data.close()

lon_co=np.where( (LON>=139.8)&(LON<=140.3) )[0]
lat_co=np.where( (LAT>=-74)&(LAT<=-24) )[0]

# Process Times
DATA_TIMES=MFDataset(data_path+'*.nc')['time']
TIME_UNIT=DATA_TIMES.units
DATA_times=num2date(DATA_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),28)
TIMES_co=np.where( (DATA_times>=Tst)&(DATA_times<=Ted) )[0]
# =============================================================================

Time_array=np.array([[i.year, i.month, i.day] for i in DATA_times])
DATA_times[0].strftime('%Y-%m')

VAR1=[]
for i in tqdm(range(len(DATA_path))):
    ncA=np.squeeze(Dataset(DATA_path[i])['temp'][:,:,lat_co,lon_co].data)
    VAR1.append(np.mean(ncA,axis=0));
VAR1=np.array(VAR1)
VAR1[VAR1<-100]=np.nan

X,Z=np.meshgrid(LAT[lat_co],-DEPTH)

My_date=pd.date_range('1980-01','2018-01',freq='Y').strftime('%Y')


# plt.pcolormesh(X,Y,DATA[0])

plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams["font.family"] = 'Arial'
cmap_list=['#333333','#7F00BF','#9200E4','#AD07FF','#C23EFF','#DA87FF','#BC0D0F','#B30000','#C91F25',\
   '#D9343E','#E84A56','#F65E67','#FF6E6E','#FF8051','#FF8B1B','#FF9C00','#FFAA09','#FFBC1E',\
       '#FFD039','#FFE256','#FFF26F','#008000','#00A400','#00D500','#1EF31E','#69FC69','#0077B3',\
           '#008DDE','#07ABFF','#3EC1FF','#87D9FF','#EEEEEE']

MyCmap=ListedColormap(cmap_list).reversed()

NN=15
data_lim=[-2.1,15]

My_levels=np.linspace(data_lim[0],data_lim[-1],NN)
MyCmap = ListedColormap(MyCmap(np.linspace(0, .95, len(My_levels)+1,endpoint=True)))
MyCmap = ListedColormap(cmocean.cm.thermal(np.linspace(0, 1, len(My_levels)+1,endpoint=True)))

Label_size=18
t=0
for i,j in zip(VAR1,My_date):
    t+=1
    Title_name='Years: '+j+f' (+{t-1:02d}~{t:02d})'
    # Figures
    fig, axs = plt.subplots(2,1,figsize=(6,4),constrained_layout = True,
                            sharex=True,gridspec_kw={'height_ratios': [1, 1.3]},dpi=200)
    axs[0].set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
    im0=axs[0].contour(X,Z,i,colors='k',levels=[1,3])
    im0.collections[1].set_linestyle('dashed')
    im1=axs[0].pcolor(X,Z,i,cmap=MyCmap,vmin=data_lim[0],vmax=data_lim[-1])
    axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
    axs[0].set_ylim(-Tcline,0)
    axs[0].set_xlim(-66.5,-42.5)
    im3=axs[1].contour(X,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
    im3.collections[1].set_linestyle('dashed')
    # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
    im2=axs[1].pcolor(X,Z,i,cmap=MyCmap,vmin=-2,vmax=24)
    axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
    axs[1].set_ylim(-5000,-Tcline)
    axs[1].set_xlim(-66.5,-42.5)
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("bottom", size="7%", pad=.35)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    if fig_bool:
        plt.savefig(w_path_sig+'ppt/'+save_name+'_'+j.replace('-','_'),
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path_sig+save_name+'_'+j.replace('-','_'))
    plt.show()










