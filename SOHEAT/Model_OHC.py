# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 17:01:09 2023

@author: shjo9
"""

import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import cartopy
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
npth='D:/OneDrive/base142/Warehouse03/H0701/OHC_H0701.nc'

OCH=xr.open_dataset(npth).loc[dict(ocean_time=slice('2013-01','2016-12'))]


OCH=OCH.assign_coords({'TT':('ocean_time',range(len(OCH.ocean_time)))})
OCH=OCH.swap_dims({"ocean_time":"TT"})
OCH_s=OCH.polyfit(dim='TT',deg=1,skipna=True)

Coef=OCH_s.OHC_polyfit_coefficients[0]
Coef_var=Coef.values*12*10 # (m/year)

lon_m,lat_m=np.meshgrid(OCH.longitude,OCH.latitude)

Plot_SO_Merc3(lon_m,lat_m,Coef_var,'Coef_var',plt.get_cmap('jet',15),[-5*10**9,5*10**9],\
              'w_path','save_name',fig_bool=False)


    
npth='G:/tmp/OHC_TEST_01.nc'

OCH=xr.open_dataset(npth).loc[dict(ocean_time=slice('2003-01','2012-01'))]


OCH=OCH.assign_coords({'TT':('ocean_time',range(len(OCH.ocean_time)))})
OCH=OCH.swap_dims({"ocean_time":"TT"})
OCH_s=OCH.polyfit(dim='TT',deg=1,skipna=True)

Coef=OCH_s.polyfit_coefficients[0]
Coef_var=Coef.values*12 # (m/year)









  
def Plot_SO_Merc3(lonA,latA,MyDATA,t_name,CMAP,Mylim,w_path,save_name,fig_bool=False):
    
    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,4),
                       subplot_kw={'projection': MERC})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 24}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular'})

    # M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    M=plt.pcolormesh(lonA, latA, MyDATA,
                  transform=PC,cmap=CMAP)
    plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
 #   ax.set_extent([0, 360, -80, -24], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()