#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 16:35:10 2023

@author: shjo
"""

import numpy as np
from netCDF4 import Dataset,MFDataset
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
from scipy.interpolate import interp2d, griddata

save_dir1='/data4/base158/Warehouse01/MyPro01/aice_2015/'

INI_path = '/data4/DATA/GLORYS/'
# INI_path2= '/data4/DATA/GLORYS/'

INI_list= [INI_path+i for i in os.listdir(INI_path) if i.startswith('mercatorglorys12v1_gl12_mean_2015')]
# GLORYS_list= [INI_path2+i for i in os.listdir(INI_path2) if i.startswith('mercatorglorys12v1_gl12_mean_1993')]

fig_bool=1


import xarray as xr
ncI=xr.open_mfdataset(INI_list).siconc


# Variables2=ncI2.variables.keys()

latrng=[-80,-24]
lonrng=[0,360]

LAT=ncI['latitude'][:]
LON=ncI['longitude'][:]

posi_co=np.where(LON>=0)[0]
nega_co=np.where(LON<0)[0]
LON_new=np.concatenate([LON[posi_co],360+LON[nega_co]])

AICE=np.concatenate([ncI[:,:,posi_co].values,ncI[:,:,nega_co].values],2)

AICE[AICE!=AICE]=0

lat_co=np.where((LAT>=latrng[0])&(LAT<=latrng[-1]))[0]
lon_co=np.where((LON_new>=lonrng[0])&(LON_new<=lonrng[-1]))[0]

lat,lon=LAT[lat_co],LON_new[lon_co]


# =============================================================================
# 
# =============================================================================
lon_m,lat_m=np.meshgrid(lon,lat)

INI_aice=AICE[:,lat_co[0]:lat_co[-1],lon_co[0]:lon_co[-1]]



from scipy import io
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
NN=13
cmap_dir='/data4/base158/Factory/'
Colorbars=io.loadmat(cmap_dir+'color_cbrewer2')['spectral']
cmap_sst=ListedColormap(Colorbars[np.arange(0,256,16),:]).reversed()


zeta_CMAP=cmocean.cm.balance
cmap_zeta=cmocean.cm.balance
latrng=[-80,-24]

sst_lim=[-3,27]
hice_lim=[0,1]
sss_lim=[32,36.5]
u_lim=[-10,10]
zeta_lim=[-2,1.5]
aice_lim=[0,1]
u_CMAP=cmap_zeta
sst_CMAP=ListedColormap(Colorbars[np.arange(0,256,16),:]).reversed()
sss_CMAP=plt.get_cmap('RdYlBu_r',15)
aice_CMAP=cmocean.cm.ice

spstere_size=(11,11)
merc_size=(11,4)
FS=28
LS=20
levels1=7
levels2=16
levels_sss1=7
levels_sss2=16

# =============================================================================
# 
# =============================================================================

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


def plot_model_data(lon_m,lat_m,data,datalim,CLABEL,title,w_path,save_name):
    fig, ax = plt.subplots(figsize=spstere_size,linewidth=1,dpi=200)
    ax = plt.gca()
    m = Basemap(projection='spstere',boundinglat=-23,lon_0=180,
                llcrnrlat=latrng[0],urcrnrlat=latrng[-1],\
                llcrnrlon=0,urcrnrlon=360,resolution='c',round=True)
    plt.title(title+' month: '+str(i+1),
          loc='left',pad=50, fontsize=FS+2,fontweight='bold')
    m.fillcontinents(color=[.75,.75,.75],lake_color='black')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,80.,20),labels=[True,False,False,False],
                    dashes=[2,2],fontsize=FS,fontweight='bold',color='grey')
    m.drawmeridians(np.arange(0.,359.,60.),labels=[True,True,True,True],
                    dashes=[2,2],fontsize=FS,fontweight='bold',color='grey')
    lon_mr,lat_mr=m(lon_m,lat_m)
    # cs1 = m.contour(lon_mr,lat_mr,data,colors='k',levels=levels1,linestyles='-.',alpha=1.)
    # plt.clabel(cs1, inline=1, fontsize=26,fmt=r'%1.1f',colors='k')
    # cs2 = m.contourf(lon_mr,lat_mr,data,cmap=CLABEL,levels=30)
    cs2 = m.pcolormesh(lon_mr,lat_mr,data,cmap=CLABEL)#,shading='gouraud')
    plt.clim(datalim[0],datalim[-1])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=1.)
    cax.tick_params(labelsize=FS)
    cax.set_ylabel('',{'fontsize':FS,'fontweight':'bold','style':'italic'})
    h = plt.colorbar(label='',cax=cax)
    plt.tight_layout()
    if fig_bool:
        plt.savefig(save_dir1+'ppt/'+save_name+'{0:02d}'.format(i+1),
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(save_dir1+save_name+'{0:02d}'.format(i+1))
    plt.show()
 
def plot_model_data_merc(lon_m,lat_m,data,datalim,CLABEL,title,w_path,save_name):
    fig, ax = plt.subplots(figsize=merc_size,linewidth=1,dpi=200)
    ax = plt.gca()
    m = Basemap(projection='merc',
                llcrnrlat=latrng[0]+.3,urcrnrlat=latrng[-1],\
                llcrnrlon=0,urcrnrlon=360,resolution='c')
    plt.title(title+' month: '+str(i+1),
          loc='left',pad=15, fontsize=FS+2,fontweight='bold')
    m.fillcontinents(color=[.75,.75,.75],lake_color='black')
    m.drawcoastlines()
    m.drawparallels(np.arange(-80.,80.,20),labels=[True,False,False,False],
                    dashes=[2,2],fontsize=FS-2,fontweight='bold',color='grey')
    m.drawmeridians(np.arange(0.,359.,60.),labels=[True,True,False,True],
                    dashes=[2,2],fontsize=FS-2,fontweight='bold',color='grey')
    lon_mr,lat_mr=m(lon_m,lat_m)
    # cs1 = m.contour(lon_mr,lat_mr,data,colors='k',levels=levels1,linestyles='-.',alpha=1.)
    # plt.clabel(cs1, inline=1, fontsize=26,fmt=r'%1.1f',colors='k')
    # cs2 = m.contourf(lon_mr,lat_mr,data,cmap=CLABEL,levels=30)
    cs2 = m.pcolormesh(lon_mr,lat_mr,data,cmap=CLABEL)#,shading='gouraud')
    plt.clim(datalim[0],datalim[-1])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=.1)
    cax.tick_params(labelsize=FS)
    cax.set_ylabel('',{'fontsize':FS,'fontweight':'bold','style':'italic'})
    h = plt.colorbar(label='',cax=cax)
    
    plt.tight_layout()
    if fig_bool:
        plt.savefig(save_dir1+'ppt/'+save_name+'{0:02d}'.format(i+1)+'merc',
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(save_dir1+save_name+'{0:02d}'.format(i+1)+'merc')
    plt.show()


for i in range(12):
    aice=INI_aice[i,:,:]

    plot_model_data(lon_m,lat_m,aice,aice_lim,aice_CMAP,'AICE',save_dir1,'aice')
    plot_model_data_merc(lon_m,lat_m,aice,aice_lim,aice_CMAP,'AICE',save_dir1,'aice')

    
    










































