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
import numpy as np

lat_rng=[-75,-30.] 
lon_rng=[153,290] # Pacific
# lon_rng=[20,117]  # Indian
# lon_rng=[-70,20]  # Atlantic 

lon_rng2=[200,220]
lat_rng2=[-60,-53]

a,b=[lon_rng[0], lon_rng[0]],[lat_rng[0],lat_rng[-1]]
c,d=[lon_rng[-1], lon_rng[-1]],[lat_rng[0],lat_rng[-1]]
e,f=[lon_rng[0], lon_rng[-1]],[lat_rng[0],lat_rng[0]]
g,h=[lon_rng[0], lon_rng[-1]],[lat_rng[-1],lat_rng[-1]]

a1,b1=[lon_rng2[0], lon_rng2[0]],[lat_rng2[0],lat_rng2[-1]]
c1,d1=[lon_rng2[-1], lon_rng2[-1]],[lat_rng2[0],lat_rng2[-1]]
e1,f1=[lon_rng2[0], lon_rng2[-1]],[lat_rng2[0],lat_rng2[0]]
g1,h1=[lon_rng2[0], lon_rng2[-1]],[lat_rng2[-1],lat_rng2[-1]]

Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
fig, ax = plt.subplots(1, 1, figsize=(7,6),
                subplot_kw={'projection': Spheric})
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)

gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                linewidth=.6, color='k', alpha=0.45, linestyle='-.')
gl.rotate_labels=False
gl.xlabels_top,gl.ylabels_right = True,True
gl.xlabel_style = gl.ylabel_style = {"size" : 26}

ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)

ax.plot(a1,b1,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(c1,d1,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(e1,f1,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
ax.plot(g1,h1,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)

ax.set_extent([0, 360, -80, -20], crs=PC)

ax.tick_params(axis='both', which='major', labelsize=28)
plt.tight_layout()
plt.show()