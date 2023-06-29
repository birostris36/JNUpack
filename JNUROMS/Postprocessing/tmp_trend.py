# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 17:10:21 2023

@author: shjo9
"""


A=xr.open_mfdataset('G:/Models/V205701_MW/Avg/*.nc').temp.loc[dict(ocean_time=slice('1990-01','1990-12'))]


# AB = A.zeta.polyfit(dim='ocean_time',deg=1)

AB = A.mean(dim='xi_rho').polyfit(dim='ocean_time',deg=1,skipna=True)



AB.polyfit_coefficients[0].plot()


# =============================================================================
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import os
# =============================================================================
import cmocean
CMAP=cmocean.cm.balance

zeta=xr.open_mfdataset('G:/Models/V205701_MW/Avg/*.nc').zeta.loc[dict(ocean_time=slice('1980-02','1990-12'))]
zeta= zeta.assign_coords({'TT':('ocean_time',range(len(zeta.ocean_time)))})
zeta=zeta.swap_dims({"ocean_time":"TT"})
z_s=zeta.polyfit(dim='TT',deg=1,skipna=True)
# z_s.polyfit_coefficients[0].plot(cmap=CMAP)





CMAP=plt.get_cmap('seismic',15)



def zonal_data_trend(data_nm,lat_rng,cmap,**kargs):
    global Avg_pth,Grd_npth, save_pth, fig_bool

    AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]
    
    
    TOPO=xr.open_dataset(Grd_npth).h.mean(dim='xi_rho').values

    
    Sample_Data=Dataset(AVGS[0])
    lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
    lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
        
    Sample=xr.open_dataset(AVGS[0])
    
    NC=xr.open_mfdataset(AVGS[0])
    data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,ocean_time=slice(kargs['st'],kargs['ed']))].mean(dim='xi_rho')
    
    data=data_.assign_coords({'TT':('ocean_time',range(len(data_.ocean_time)))})
    data=data.swap_dims({"ocean_time":"TT"})
    data_s=data.polyfit(dim='TT',deg=1,skipna=True)
    Coef=data_s.polyfit_coefficients[0]
    Coef_var=Coef.values*12 # (m/year)
    
    My_lim=np.nanmean(Coef_var)+np.nanstd(Coef_var)*3
    my_lim=[-My_lim,My_lim]
    
    Z=jr.zlevs(NC['Vtransform'].values, NC['Vstretching'].values,NC['theta_s'].values,\
           NC['theta_b'].values, NC['Tcline'].values,  data_.s_rho.shape[0] ,1, TOPO, np.zeros_like(TOPO))


    lat=NC.lat_rho.values[:,0]
    lat_m,z_m=np.meshgrid(lat,Z[:,0])

    t_name=''   
    s_name_S=''

    Label_size=12
    xtick_location = np.linspace(lat[0], lat[-1],6)
    xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]
        
    fig, axs = plt.subplots(2,1,figsize=(6,4),
                            sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
    # fig.subplots_adjust(wspace=0, hspace=0)
    axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
  #  im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
  #  axs[0].clabel(im0, inline=1, fontsize=14)
  #  im0.collections[1].set_linestyle('dashed')
    im1=axs[0].pcolor(lat_m,Z,Coef_var,cmap=cmap,vmin=my_lim[0],vmax=my_lim[-1])
    axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
    axs[0].set_ylim(-NC['Tcline'].values,0)
    axs[0].set_xlim(-80,-23.5)
  #  im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
#    axs[1].clabel(im0, inline=1, fontsize=14)
    axs[0].set_xticks(ticks=xtick_location)
    axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
    
    # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
    im2=axs[1].pcolor(lat_m,Z,Coef_var,cmap=cmap,vmin=my_lim[0],vmax=my_lim[-1])
    axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
    axs[1].set_ylim(-5000,-NC['Tcline'].values)
    axs[1].set_xlim(-80,-23.5)
    axs[1].set_xticks(ticks=xtick_location)
    axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("bottom", size="7%", pad=.35)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    if fig_bool:
        plt.savefig(save_pth+'Zonal_trend_'+'_'+data_nm+'/ppt/'+s_name_S,
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(save_pth+'Zonal_trend_'+'_'+data_nm+'/'+s_name_S,bbox_inches='tight')
    plt.show()


def Surface_data_trend(data_nm,lat_rng,cmap,**kargs):
    global Avg_pth, save_pth, fig_bool
    
    AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]
    
    Sample_Data=Dataset(AVGS[0])
    lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
    lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
        
    Sample=xr.open_dataset(AVGS[0])
    
    if [i for i in Sample[data_nm].coords].count('s_rho'):
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(s_rho=Sample.s_rho.values[-1],\
                                                        eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]     
    else:
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]
    # Linear trend
    data=data_.assign_coords({'TT':('ocean_time',range(len(data_.ocean_time)))})
    data=data.swap_dims({"ocean_time":"TT"})
    data_s=data.polyfit(dim='TT',deg=1,skipna=True)
    Coef=data_s.polyfit_coefficients[0]
    Coef_var=Coef.values*12 # (m/year)
    
    My_lim=np.nanmean(Coef_var)+np.nanstd(Coef_var)*1/2
    my_lim=[-My_lim,My_lim]
    
    t_name=''
    s_name_S=''
    save_pth=''
    Plot_SO_Spherical3(data.lon_rho,data.lat_rho,Coef,\
                          t_name,CMAP,my_lim,\
                          save_pth+'Surface_mean_'+'_'+data_nm,s_name_S)
        
    # Plot_SO_Merc3(lon_rho,lat_rho,i,t_name,My_levels,cmap,data_lim,\
    #                   save_pth+'Surface_mean_'+kargs['mean']+'_'+data_nm\
    #                   ,s_name_M,fig_bool)
    
Surface_data_trend('zeta',[-80,-24],CMAP,st='0001-01',ed='0120-12')





