#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 13:01:33 2023

@author: shjo
"""
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date, date2num,MFDataset
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
sys.path.append('D:/OneDrive/JNUpack/')
from Mapping.Tools import d_modules as mm

AVGS_path='G:/Models/TK0525ED_CLM/'
save_path='D:/OneDrive/JNUpack/Mapping/Samples_figs/'

lat_rng=[-80,-24]

AVGS=[AVGS_path+i for i in os.listdir(AVGS_path) if i.endswith('.nc')]
Sample_Data=Dataset(AVGS[0])
lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])

TEMP=Sample_Data['temp'][9,-1,lat_co,:]

plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "regular"
plt.rcParams['axes.linewidth'] = 1
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True

def Plot_SO_Spherical2(lonA,latA,MyDATA,My_levels,CMAP,Mylim,w_path,save_name,fig_bool=False):

    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                       subplot_kw={'projection': Spheric})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.xlabels_top,gl.ylabels_right = True,True
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)

    # plt.contour(lonA,latA,MyDATA,levels=My_levels_dw,transform=PC,colors='w')
    # plt.contour(lonA,latA,MyDATA,levels=My_levels_up,transform=PC,colors='k')
    M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)

    #plt.pcolormesh(lonA, latA, MyDATA,
    #              transform=PC,cmap=CMAP)
    plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
    ax.set_extent([0, 360, -80, -24], crs=PC)
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
    cb.set_label(label='m', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+save_name)
    plt.show()
    
# =============================================================================
# Cartopy Mercator
# =============================================================================

import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt
import cmocean

data_lim=[-2,28]

NN=15
My_levels_dw=np.arange(data_lim[0],0,0.2)
My_levels_up=np.arange(0,data_lim[-1],2)

My_levels=np.concatenate((My_levels_dw, My_levels_up),axis=0)

# My_levels=np.linspace(data_lim[0],data_lim[-1],len(My_levels_dw)+len(My_levels_up))


# My_levels=np.arange(data_lim[0],data_lim[-1]+NN/2,NN)
zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, 2*len(My_levels),endpoint=True)))
# zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, len(My_levels),endpoint=True)))

zeta_CMAP=zeta_CMAP =  cmocean.cm.balance


Plot_SO_Spherical2(lon_rho,lat_rho,TEMP,My_levels,zeta_CMAP,data_lim,save_path,'Spherical_sst_sample',fig_bool=True)
# mm.Plot_SO_Merc2(lon_rho,lat_rho,TEMP,My_levels,zeta_CMAP,data_lim,save_path,'Merc_sst_sample',fig_bool=True)



















