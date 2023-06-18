# -*- coding: utf-8 -*-
"""
Created on Tue May 16 20:24:42 2023

@author: shjo9
"""

import sys 
sys.path.append('D:/OneDrive/JNUpack/JNUROMS/')
import Tools.JNUROMS as jr
from Mapping.Tools import d_modules as mm
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
import xarray as xr
import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from scipy.interpolate import interp2d, griddata
import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime as dt
import cmocean
from Tools.dianose_module import *
Avg_pth='G:/Models/TK0525ED_CLM/'
Log_pth='G:/TEST/LogTK25EM_DQD.txt'
save_pth='G:/TEST/'
SODA_pth='G:/SODA/'
Grd_npth='G:/MODEL_DATA/Grd/Grd_SO_05d_sponge.nc'
MyCMAP=''
fig_bool=0

plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams["font.family"] = 'Arial'
mpl.rcParams['axes.unicode_minus'] = False
# Weighted monthly mean
def season_mean(ds, calendar="standard"):
    
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = ds['ocean_time'].dt.days_in_month

    # Calculate the weights by grouping by 'time.season'
    weights = (
        month_length.groupby("ocean_time.season") / month_length.groupby("ocean_time.season").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby("ocean_time.season").sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby("ocean_time.season").sum(dim="ocean_time")
def yearly_mean(ds, calendar="standard"):
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = ds['ocean_time'].dt.days_in_month

    # Calculate the weights by grouping by 'time.season'
    weights = (
        month_length.groupby("ocean_time.year") / month_length.groupby("ocean_time.year").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    # np.testing.assert_allclose(weights.groupby("ocean_time.year").sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby("ocean_time.year").sum(dim="ocean_time")



def zonal_data_diff_Soda(data_nm,cmap,cmap1,data_lim,data_lim1,**kargs):
    global Avg_pth,Grd_npth, save_pth,fig_bool, SODA_pth
    
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams["font.weight"] = "regular"


    t_rng=[kargs['st'],kargs['ed']]

    
    if data_nm=='zeta':
        date_soda_nm='ssh'
    elif data_nm=='u_eastward':
        data_soda_nm='u'
        # data_lim1=[-0.3,0.3]
    elif data_nm=='v_northward':
        data_soda_nm='v'
        # data_lim1=[-0.1,0.1]
    elif data_nm=='temp':
        data_soda_nm='temp'
        # data_lim1=[-3.5,3.5]
    elif data_nm=='salt':
        data_soda_nm='salt'
        # data_lim1=[-1,1.]


    # Proceses SODA
    SODA=xr.open_mfdataset(SODA_pth+'*.nc')[data_soda_nm].\
        loc[dict(yt_ocean=slice(-80,-23.5),time=slice(t_rng[0],t_rng[-1]))].rename({'time':'ocean_time'})
    zonal_soda_m=SODA.mean(dim='xt_ocean')
    SODA_lat=SODA.yt_ocean.values
    SODA_Z=SODA.st_ocean.values
    SODA_lat_m,SODA_Z_m=np.meshgrid(SODA_lat,SODA_Z)
    


    save_name='Zonal_temp_average_section'
    
    
    # Read Grd
    TOPO=xr.open_dataset(Grd_npth).h.mean(dim='xi_rho')
    
    
    AVGS=np.sort([Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')])
    NC=xr.open_mfdataset(AVGS)
    DATA=NC[data_nm].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]
    zeta=NC['zeta'].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]

    zonal_m=DATA.mean(dim='xi_rho')
    zonal_zeta_m=zeta.mean(dim='xi_rho')
    
    lat=NC.lat_rho.values[:,0]
    
    if kargs['mean']=='ann':
        data=yearly_mean(zonal_m).rename({'year':'ocean_time'})
        data_soda=yearly_mean(zonal_soda_m).rename({'year':'ocean_time'})
        d_zeta=yearly_mean(zonal_zeta_m).rename({'year':'ocean_time'})
    elif kargs['mean']=='season':
        data=season_mean(zonal_m).rename({'season':'ocean_time'})
        data_soda=season_mean(zonal_soda_m).rename({'season':'ocean_time'})
        d_zeta=season_mean(zonal_zeta_m).rename({'season':'ocean_time'})
    elif kargs['mean']=='monthly':
        data=zonal_m.resample(ocean_time='1MS').mean()
        data_soda=zonal_soda_m.resample(ocean_time='1MS').mean()
        d_zeta=zonal_zeta_m.resample(ocean_time='1MS').mean()
    elif kargs['mean']=='monthly_clm':
        data=zonal_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        data_soda=zonal_soda_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        d_zeta=zonal_zeta_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
    elif kargs['mean']=='total':
        data=zonal_m.mean(dim='ocean_time',keepdims=True)
        data_soda=zonal_soda_m.mean(dim='ocean_time',keepdims=True)
        d_zeta=zonal_zeta_m.mean(dim='ocean_time',keepdims=True)

    for i in d_zeta.values:
        Z=jr.zlevs(NC['Vtransform'].values[0], NC['Vstretching'].values[0],NC['theta_s'].values[0],\
               NC['theta_b'].values[0], NC['Tcline'].values[0], NC.s_rho.shape[0],1, TOPO.values, i)
    
    Label_size=12
    xtick_location = np.linspace(lat[0], lat[-1],6)
    xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]
    

    lat_m,z_m=np.meshgrid(lat,Z[:,0])
    
    
    # Interpolates SODA grid to model grid
    from scipy.interpolate import interp2d, griddata

 
    for i,j in zip(data,data_soda):
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)
        
        
        tmp_soda_=griddata( (SODA_lat_m.flatten(),-SODA_Z_m.flatten()),j.values.flatten(),
                (lat_m.flatten(),Z.flatten() ),
            method='linear',fill_value=np.nan)
        
        soda_re = tmp_soda_.reshape(lat_m.shape)

        model_soda=i-soda_re
        
        s_name_S1='Zonal_mean_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+'Model_SODA'+'_'+kargs['ed'].replace('-','')
        s_name_S='Zonal_mean_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+'SODA'+'_'+kargs['ed'].replace('-','')

        # Figure 1) Model - SODA  
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        # fig.subplots_adjust(wspace=0, hspace=0)
        axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
        # im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[1,3])
        im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(lat_m,Z,model_soda,cmap=cmap1,vmin=data_lim1[0],vmax=data_lim1[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-NC['Tcline'].values[0],0)
        axs[0].set_xlim(-80,-23.5)
        # im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
        # im3.collections[1].set_linestyle('dashed')
        axs[0].set_xticks(ticks=xtick_location)
        axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(lat_m,Z,model_soda,cmap=cmap1,vmin=data_lim1[0],vmax=data_lim1[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-5000,-NC['Tcline'].values[0])
        axs[1].set_xlim(-80,-23.5)
        axs[1].set_xticks(ticks=xtick_location)
        axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if fig_bool:
            plt.savefig(save_pth+'ppt/'+s_name_S,
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(save_pth+s_name_S,bbox_inches='tight')
        plt.show()
        
        # Figure 2) SODA  
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        # fig.subplots_adjust(wspace=0, hspace=0)
        axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
        # im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[1,3])
        im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(SODA_lat_m,-SODA_Z_m,j,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-NC['Tcline'].values[0],0)
        axs[0].set_xlim(-80,-23.5)
        # im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
        # im3.collections[1].set_linestyle('dashed')
        axs[0].set_xticks(ticks=xtick_location)
        axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(SODA_lat_m,-SODA_Z_m,j,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-6000,-NC['Tcline'].values[0])
        axs[1].set_xlim(-80,-23.5)
        axs[1].set_xticks(ticks=xtick_location)
        axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if fig_bool:
            plt.savefig(save_pth+'ppt/'+s_name_S,
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(save_pth+s_name_S,bbox_inches='tight')
        plt.show()      
        
        
      
        
      
        
temp_NN=15
data_lim=[-2.5,28]
data_lim1=[-3.5,3.5]

temp_levels=np.arange(data_lim[0],data_lim[-1]+2/2,2.)
temp_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')(np.linspace(0, 1, temp_NN,endpoint=True)))
            

zonal_data_diff_Soda('temp',temp_CMAP,temp_CMAP,data_lim,data_lim1,\
                 mean='monthly_clm',st='2010-02',ed='2010-12')

    
    
def data_drift(data_nm,lat_rng,My_levels,cmap,data_lim,**kargs):
    global Avg_pth, Log_pth, save_pth,fig_bool

    AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]
    
    Sample_Data=Dataset(AVGS[0])
    lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
    lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
        
    Sample=xr.open_dataset(AVGS[0])
    Sample.s_rho.values
    
    if [i for i in Sample[data_nm].coords].count('s_rho'):
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(s_rho=Sample.s_rho.values[-1],\
                                                        eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]
    else:
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]

    if kargs['mean']=='ann':
        data=yearly_mean(data_).rename({'year':'ocean_time'})
    elif kargs['mean']=='season':
        data=season_mean(data_).rename({'season':'ocean_time'})
    elif kargs['mean']=='monthly':
        data=data_.resample(ocean_time='1MS').mean()
    elif kargs['mean']=='monthly_clm':
        data=data_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
    elif kargs['mean']=='total':
        data=data_.mean(dim='ocean_time',keepdims=True)
        
    for i in data:
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)
        s_name_S='Spherical_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
        
        s_name_M=s_name_S.replace('Spherical','Merc')
        
        mm.Plot_SO_Spherical2(data.lon_rho.values,data.lat_rho.values,i.squeeze().values,\
                              t_name,My_levels,cmap,data_lim,save_pth,s_name_S,fig_bool)
            
        mm.Plot_SO_Merc2(lon_rho,lat_rho,i,t_name,My_levels,cmap,data_lim,\
                          save_pth,s_name_M,fig_bool)
    
    
    
    
    
        
        
        
        