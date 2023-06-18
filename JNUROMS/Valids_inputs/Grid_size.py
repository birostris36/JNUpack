# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 11:24:28 2023

@author: shjo9
"""

import cartopy
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime as dt
import cmocean
import xarray as xr
import numpy as np

lon_rng,lat_rng=[23,24],[-75,-60]

TKGRD=xr.open_dataset('G:/MODEL_DATA/TK_Inputs/roms_grd_sao025.nc')
MyGRD=xr.open_dataset('G:/MODEL_DATA/Grd/Grd_Q0_Rtopo30S_Smooth_sponge.nc')

TK_lon,TK_lat = TKGRD.lon_rho.values[0,:], TKGRD.lat_rho.values[:,0]
MY_lon,MY_lat = MyGRD.lon_rho.values[0,:], MyGRD.lat_rho.values[:,0]

TK_lon_co=np.where( (TK_lon>=lon_rng[0]) & (TK_lon<=lon_rng[-1]))[0]
TK_lat_co=np.where( (TK_lat>=lat_rng[0]) & (TK_lat<=lat_rng[-1]))[0]

MY_lon_co=np.where( (MY_lon>=lon_rng[0]) & (MY_lon<=lon_rng[-1]))[0]
MY_lat_co=np.where( (MY_lat>=lat_rng[0]) & (MY_lat<=lat_rng[-1]))[0]

My_H=MyGRD.h.loc[dict(xi_rho=MY_lon_co,eta_rho=MY_lat_co)]
TK_H=TKGRD.h.loc[dict(xi_rho=TK_lon_co,eta_rho=TK_lat_co)]

# My_H.plot()
# TK_H.plot()

topo_NN=16
topo_lim=[-91,0]
topo_levels=np.arange(topo_lim[0],topo_lim[-1]+100/2,100)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
topo_CMAP = ListedColormap(cmocean.cm.topo(np.linspace(0, 1, len(topo_levels)+1,endpoint=True)))

topo_CMAP = ListedColormap(['w','w','w','w','w','w','w','w']*32)


MY_lon_m,MY_lat_m=np.meshgrid(MY_lon[MY_lon_co],MY_lat[MY_lat_co])
TK_lon_m,TK_lat_m=np.meshgrid(TK_lon[TK_lon_co],TK_lat[TK_lat_co])


# Plot_SO_Merc2(MY_lon_m,MY_lat_m,-My_H,'',topo_CMAP,topo_lim,'','',False)
# Plot_SO_Spherical2(TK_lon_m,TK_lat_m,-TK_H,'',topo_CMAP,'','',False)



Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)

# Now we will create axes object having specific projection 
############################################################
fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                   subplot_kw={'projection': Spheric},dpi=200)
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
# To plot borders and coastlines, we can use cartopy feature
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
# ax.set_title(t_name,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                  linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                      )
gl.rotate_labels=False
gl.xlabels_top,gl.ylabels_right = True,True
gl.xlabel_style = gl.ylabel_style = {"size" : 26}
 
# M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
M=plt.pcolor(MY_lon_m, MY_lat_m, My_H,
              transform=PC,edgecolors='k',cmap=topo_CMAP)
# plt.clim(Mylim[0],Mylim[-1])

# crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
ax.set_extent([lon_rng[0]-1, lon_rng[-1]+1, lat_rng[0]-1, lat_rng[-1]+1], crs=PC)

ax.tick_params(axis='both', which='major', labelsize=28)

plt.tight_layout()
if 0:
    plt.savefig(w_path+'/ppt/'+save_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+'/'+save_name)
plt.show()
############################################################


############################################################
fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                   subplot_kw={'projection': Spheric},dpi=200)
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
# To plot borders and coastlines, we can use cartopy feature
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
# ax.set_title(t_name,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                  linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                      )
gl.rotate_labels=False
gl.xlabels_top,gl.ylabels_right = True,True
gl.xlabel_style = gl.ylabel_style = {"size" : 26}
 
# M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
M=plt.pcolor(TK_lon_m, TK_lat_m, TK_H,
              transform=PC,edgecolors='k',cmap=topo_CMAP)
# plt.clim(Mylim[0],Mylim[-1])

# crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
ax.set_extent([lon_rng[0]-1, lon_rng[-1]+1, lat_rng[0]-1, lat_rng[-1]+1], crs=PC)

ax.tick_params(axis='both', which='major', labelsize=28)

plt.tight_layout()
if 0:
    plt.savefig(w_path+'/ppt/'+save_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+'/'+save_name)
plt.show()
############################################################
