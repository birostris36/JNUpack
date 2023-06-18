#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:31:35 2023

@author: shjo
"""


import numpy as np
from netCDF4 import Dataset,MFDataset, num2date, date2num
import os
import xarray as xr
import dask
import matplotlib as mpl
import cmocean
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import datetime as dt
from tqdm import tqdm

AVG_path='/data4/base158/CaseSODA/TK25EM/Outputs/Avg/'
# save_dir1='/data4/base158/Warehouse01/MyPro01/TK0525EM/TK25EM/ZETA/'
fig_bool=0
latrng=[-80,-24]
latrng=[-80,-24]


AVGS=[AVG_path+'/'+i for i in os.listdir(AVG_path) if i.endswith('.nc')]

Sample_nc=Dataset(AVGS[0])
latA,lonA=Sample_nc['lat_rho'][:], Sample_nc['lon_rho'][:]

tmp_time_var='ocean_time'
t_rng=['2011-01','2014-12']
My_time_ref='days since 1950-1-1 00:00:00'
AVGS_TIMES=MFDataset(AVG_path+'*.nc')[tmp_time_var]
TIME_UNIT=AVGS_TIMES.units
AVGS_times=num2date(AVGS_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),28)
TIMES_co=np.where( (AVGS_times>=Tst)&(AVGS_times<=Ted) )[0]
# =============================================================================

Zeta = xr.open_dataset(AVGS[0]).zeta.loc[dict(xi_rho=slice(int(350/2),int(370/2)))]

from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import colors,colorbar

cmap_list=['w','w','w','w','w','w',]
MyCmap=ListedColormap(cmap_list) #.reversed()



plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True

plt.rcParams["font.family"] = 'Helvetica'
mpl.rcParams['axes.unicode_minus'] = False

zeta_CMAP=cmocean.cm.balance
zeta_lim=[-2,1.5]
spstere_size=(11,11)
merc_size=(11,4)
FS=28
LS=20
levels1=7
levels2=16
levels_sss1=7
levels_sss2=16


def plot_model_data(lon_m,lat_m,data,datalim,CLABEL,title,w_path,save_name):
    fig, ax = plt.subplots(figsize=spstere_size,linewidth=1,dpi=200)
    ax = plt.gca()
    m = Basemap(projection='spaeqd',
                llcrnrlat=latrng[0],urcrnrlat=latrng[-1],boundinglat=10,lon_0=270,\
                llcrnrlon=170,urcrnrlon=190,resolution='c')
    plt.title(title+'\n',
          loc='center',pad=50, fontsize=FS+2,fontweight='bold')
    m.fillcontinents(color=[.75,.75,.75],lake_color='grey')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,80.,20),labels=[True,False,False,False],
                    dashes=[2,2],fontsize=FS,fontweight='bold',color='grey')
    m.drawmeridians(np.arange(0.,359.,60.),labels=[True,True,True,True],
                    dashes=[2,2],fontsize=FS,fontweight='bold',color='grey')
    lon_mr,lat_mr=m(lon_m,lat_m)
    # cs1 = m.contour(lon_mr,lat_mr,data,colors='k',levels=levels1,linestyles='-.',alpha=1.)
    # plt.clabel(cs1, inline=1, fontsize=26,fmt=r'%1.1f',colors='k')
    # cs2 = m.contourf(lon_mr,lat_mr,data,cmap=CLABEL,levels=30)
    cs2 = m.pcolor(lon_mr,lat_mr,data,cmap=CLABEL,edgecolors='k', linewidths=.5)#,shading='gouraud')
    plt.clim(datalim[0],datalim[-1])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=1.)
    cax.tick_params(labelsize=FS)
    cax.set_ylabel('',{'fontsize':FS,'fontweight':'bold','style':'italic'})
    h = plt.colorbar(label='',cax=cax)
    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+save_name)
    plt.show()

for i in TIMES_co[0]:
    tmp_data=ZETA[i].data
    Title_name=AVGS_times[i].strftime('%Y-%m-%d')
    
    plot_model_data(Zeta.lon_rho,Zeta.lat_rho,Zeta.squeeze(),zeta_lim,MyCmap,'Zeta '+\
                    Title_name,save_dir1,'zeta_'+Title_name.replace('-',''))

plt.pcolor(Zeta.lon_rho,Zeta.lat_rho,Zeta.squeeze(),cmap=zeta_CMAP,edgecolors='k', linewidths=.5) #,shading='gouraud')
plt.clim(zeta_lim)

# plt.scatter(lonA,latA)

# Title_name=AVGS_times[371].strftime('%Y-%m-%d')

tmp_data[tmp_data==tmp_data]=np.nan





