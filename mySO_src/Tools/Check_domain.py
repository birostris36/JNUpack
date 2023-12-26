# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:19:32 2023

@author: shjo9
"""
import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime as dt
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from scipy.interpolate import interp2d, griddata
import cartopy
import matplotlib.path as mpath
from copy import deepcopy

# lat_rng,lon_rng=[18,25.5],[127.,140]
# lat_rng,lon_rng=[32.18810,53.56719],[127.2177,143.9147]

lat_rng,lon_rng=[-80,52.],[115.,150]
# lat_rng,lon_rng=[32.18810,53.56719],[127.2177,143.9147]

a,b=[lon_rng[0], lon_rng[0]],[lat_rng[0],lat_rng[-1]]
c,d=[lon_rng[-1], lon_rng[-1]],[lat_rng[0],lat_rng[-1]]
e,f=[lon_rng[0], lon_rng[-1]],[lat_rng[0],lat_rng[0]]
g,h=[lon_rng[0], lon_rng[-1]],[lat_rng[-1],lat_rng[-1]]

Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)


fs=18
 # Now we will create axes object having specific projection 
fig, ax = plt.subplots(1, 1, figsize=(8,6),
                   subplot_kw={'projection': MERC})

gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                  linewidth=.6, color='k', alpha=0.45, linestyle='-.')
gl.xlabels_top,gl.ylabels_right = False,False
gl.xlabel_style = gl.ylabel_style = {"size" : fs}

# To plot borders and coastlines, we can use cartopy feature
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
# ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular'})
# plt.plot([T_point1[0],T_point1[-1]],[T_point2[0], T_point2[-1]])
ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)

# crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
ax.set_extent([lon_rng[0]-8, lon_rng[-1]+8, lat_rng[0]-8, lat_rng[-1]+8], crs=PC)
ax.tick_params(axis='both', which='major', labelsize=fs)
plt.tight_layout()
if 0:
    plt.savefig(w_path+'/ppt/'+save_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+'/'+save_name)
plt.show()