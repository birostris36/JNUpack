#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 21:00:40 2023

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

save_dir1='/data4/base158/Warehouse01/INI_1993_ECCO2/'

INI_path = '/data1/ECCO2/quart_90S_90N/'

latrng=[-80,-24]

fig_bool=1

# INI_path2= '/data4/DATA/GLORYS/'

SSH_list= [INI_path+'SSH_monthly.nc/'+i for i in os.listdir(INI_path+'SSH_monthly.nc/')\
           if i.startswith('SSH.1440x720.1993')]
TEMP_list= [INI_path+'THETA_monthly.nc/'+i for i in os.listdir(INI_path+'THETA_monthly.nc/')\
           if i.startswith('THETA.1440x720x50.1993')]
SALT_list= [INI_path+'SALT_monthly.nc/'+i for i in os.listdir(INI_path+'SALT_monthly.nc/')\
           if i.startswith('SALT.1440x720x50.1993')]
U_list= [INI_path+'UVEL_monthly.nc/'+i for i in os.listdir(INI_path+'UVEL_monthly.nc/')\
           if i.startswith('UVEL.1440x720x50.1993')]
V_list= [INI_path+'VVEL_monthly.nc/'+i for i in os.listdir(INI_path+'VVEL_monthly.nc/')\
           if i.startswith('VVEL.1440x720x50.1993')]


# xarray(open_mfDataset), netCDF4(MFDataset) --> Dimension aggregation problem 

SampleNC=Dataset(SSH_list[0])
Sample_SSH_=SampleNC['SSH'][:].squeeze()

LON,LAT=SampleNC['LONGITUDE_T'][:],SampleNC['LATITUDE_T'][:]
lat=LAT[(LAT>=latrng[0])&(LAT<=latrng[-1])]
lat_co=np.where((LAT>=latrng[0])&(LAT<=latrng[-1]))[0]

Sample_SSH=Sample_SSH_[lat_co,:]

at,on=Sample_SSH.shape
# Read SSH
SSH_1993=np.zeros([12,at,on])
ii=0;
for i in np.sort(SSH_list):
    SSH_1993[ii]=Dataset(i)['SSH'][0,lat_co,:]
    ii+=1
# Read SST
SST_1993=np.zeros([12,at,on])
ii=0;
for i in np.sort(TEMP_list):
    SST_1993[ii]=Dataset(i)['THETA'][0,0,lat_co,:]
    ii+=1
SST_1993[SST_1993<-100]=np.nan
# Read SSS
SSS_1993=np.zeros([12,at,on])
ii=0;
for i in np.sort(SALT_list):
    SSS_1993[ii]=Dataset(i)['SALT'][0,0,lat_co,:]
    ii+=1
SSS_1993[SSS_1993<-100]=np.nan
# Read UVEL
UVEL_1993=np.zeros([12,at,on])
ii=0;
for i in np.sort(U_list):
    UVEL_1993[ii]=Dataset(i)['UVEL'][0,0,lat_co,:]
    ii+=1
UVEL_1993[UVEL_1993<-100]=np.nan
# Read UVEL
VVEL_1993=np.zeros([12,at,on])
ii=0;
for i in np.sort(V_list):
    VVEL_1993[ii]=Dataset(i)['VVEL'][0,0,lat_co,:]
    ii+=1
VVEL_1993[VVEL_1993<-100]=np.nan

# np.save('/data2/ECCO2_SST_1993.npy',SST_1993)
# np.save('/data2/ECCO2_LAT_1993.npy',lat)
# np.save('/data2/ECCO2_LON_1993.npy',LON)



# =============================================================================
# 
# =============================================================================
lon_m,lat_m=np.meshgrid(LON,lat)




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
    u=UVEL_1993[i,:,:]
    v=VVEL_1993[i,:,:]
    sst=SST_1993[i,:,:]
    sss=SSS_1993[i,:,:]
    zeta=SSH_1993[i,:,:]
    
    plot_model_data(lon_m,lat_m,u*10,u_lim,u_CMAP,'Surface U',save_dir1,'surface_u')
    plot_model_data(lon_m,lat_m,v*10,u_lim,u_CMAP,'Surface V',save_dir1,'surface_v')

    plot_model_data(lon_m,lat_m,sst,sst_lim,sst_CMAP,'SST',save_dir1,'sst')
    plot_model_data(lon_m,lat_m,zeta,zeta_lim,zeta_CMAP,'Zeta',save_dir1,'zeta')

    plot_model_data(lon_m,lat_m,sss,sss_lim,sss_CMAP,'SSS',save_dir1,'sss')
    
    plot_model_data_merc(lon_m,lat_m,u*10,u_lim,u_CMAP,'Surface U',save_dir1,'surface_u')
    plot_model_data_merc(lon_m,lat_m,v*10,u_lim,u_CMAP,'Surface V',save_dir1,'surface_v')

    plot_model_data_merc(lon_m,lat_m,sst,sst_lim,sst_CMAP,'SST',save_dir1,'sst')
    plot_model_data_merc(lon_m,lat_m,zeta,zeta_lim,zeta_CMAP,'Zeta',save_dir1,'zeta')

    plot_model_data_merc(lon_m,lat_m,sss,sss_lim,sss_CMAP,'SSS',save_dir1,'sss')













































