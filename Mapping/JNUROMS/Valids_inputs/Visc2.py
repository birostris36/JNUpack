# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:16:18 2023

@author: shjo9
"""

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

GRD=xr.open_dataset('G:/MODEL_DATA/TK_Inputs/roms_grd_sao025.nc')

pn,pm=GRD.pn,GRD.pm

om_r,on_r=1/pm,1/pn

grdscl=np.sqrt(om_r*on_r)

viscosity2=100

grdmax=np.max([np.max(om_r),np.max(on_r)])

cff=viscosity2/grdmax
visc2_r=cff*grdscl


NN=15
data_lim=[20,101]
data_levels=np.arange(data_lim[0],data_lim[-1]+5/2,5)
data_CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, len(data_levels),endpoint=True)))



Plot_SO_Spherical2(GRD.lon_rho,GRD.lat_rho,visc2_r,'visc2',data_levels,data_CMAP,data_lim,'','',fig_bool=False)

def Plot_SO_Spherical2(lonA,latA,MyDATA,t_name,My_levels,CMAP,Mylim,w_path,save_name,fig_bool=False):

    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                       subplot_kw={'projection': Spheric})
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.rotate_labels=False
    gl.xlabels_top,gl.ylabels_right = True,True
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
     
    M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    #M=plt.pcolormesh(lonA, latA, MyDATA,
    #              transform=PC,cmap=CMAP,vmin=Mylim[0],vmax=Mylim[-1])
    plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
    ax.set_extent([0, 360, -80, -24], crs=PC)
    
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()