# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 16:29:41 2023

@author: shjo9
"""
import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime as dt
import cmocean
import sys 
sys.path.append('D:/JNUpack/')
sys.path.append('D:/JNUpack/JNUROMS')
import Tools.JNUROMS as jr
import Tools.Inputs as ti
# from Mapping.Tools import d_modules as mm
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
# from pptx import Presentation # 라이브러리 
# from pptx.util import Inches,Cm, Pt # 사진, 표등을 그리기 위해
# from pptx.enum.text import PP_ALIGN
import cartopy
import matplotlib.path as mpath

Drake=[]
Agulas=[]

def Plot_SO_Merc2(w_path,save_name,fig_bool=False):
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    Cyll=ccrs.PlateCarree(central_longitude=180.0,globe=None)
    
    a,b=[295, 295],[-55,-65]

    # Now we will create axes object having specific projection 
    fig, ax = plt.subplots(1, 1, figsize=(12.5,4),
                       subplot_kw={'projection': MERC})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.')
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 24}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    # ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular'})
    # plt.plot([T_point1[0],T_point1[-1]],[T_point2[0], T_point2[-1]])
    ax.plot(a,b,transform=PC,color='k',linestyle='--')

    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
    ax.set_extent([0, 360, -80, -24], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)
    plt.tight_layout()
    if 0:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()
    

lon_rng= 120    #295
lat_rng= [-65,-30]       # [-65,-55]

def Transport_SV_NS(lon_rng,lat_rng):
    global Avg_pth,Grd_npth, save_pth, fig_bool

    AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]

    Sample_Data=Dataset(AVGS[0])
    lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
    lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    
    dl_x=np.min(np.diff(lon_rho))

    lon_co=np.where((lon_rho[0,:]>lon_rng-dl_x)&(lon_rho[0,:]<lon_rng+dl_x))[0]
    
    if len(lon_co) !=1:
        print('!!!!!!!!!!!!!!!!!!')
        raise
    
    lon_rho,lat_rho=np.meshgrid(lon_rho[0,lon_co],lat_rho[lat_co,0])
    NC=xr.open_mfdataset(AVGS[0])

    U=xr.open_mfdataset(AVGS)['u_eastward'].loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze()
    ZETA=xr.open_mfdataset(AVGS)['zeta'].loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze().values
    TOPO=xr.open_dataset(Grd_npth).h.loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze().values
    PN=xr.open_dataset(Grd_npth).pn.loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze().values

    T,th,at=U.shape

    SV=[]

    for i,j in zip(ZETA,U.values):
        Z=jr.zlevs(NC['Vtransform'].values, NC['Vstretching'].values,NC['theta_s'].values,\
               NC['theta_b'].values, NC['Tcline'].values,  U.s_rho.shape[0] ,5, TOPO, ZETA[-1])
        
        Z_=jr.zlevs(NC['Vtransform'].values, NC['Vstretching'].values,NC['theta_s'].values,\
               NC['theta_b'].values, NC['Tcline'].values,  U.s_rho.shape[0] ,1, TOPO, ZETA[-1])
            
        lat_m,_=np.meshgrid(lat_rho,Z_[:,0])
    
        X_dist=np.tile(1/PN,th).reshape([th,at])
        delta_D=Z[1:,:]-Z[:-1,:]
    
        CELL=delta_D*X_dist #(Unit: m*m=m**2)
        CELL=CELL*10**-6
    
        SV_sum=np.nansum(CELL*j)
    
        SV.append(SV_sum)
    SV=np.array(SV)







    plt.figure()
    plt.pcolor(lat_m,Z_,SV,cmap=plt.get_cmap('seismic'), edgecolor='k',vmin=-2,vmax=2)
    plt.ylim(-100,0)
    plt.colorbar()


    plt.figure()
    plt.pcolor(lat_m,Z_,U[0],cmap=plt.get_cmap('seismic'),vmin=-0.3,vmax=0.3)
    plt.colorbar()














