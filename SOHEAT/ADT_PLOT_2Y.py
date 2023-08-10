# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 09:42:19 2023

@author: shjo9
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.path as mpath
from mpl_toolkits.axes_grid1 import make_axes_locatable

pth='G:/AVISO/madt_h/'


Sample=xr.open_dataset(pth+'dt_global_allsat_madt_h_y2002_m04.nc')

NC=xr.open_mfdataset(pth+'*.nc').loc[dict(latitude=slice(-65,40),nv=0)].adt


NC=NC.assign_coords({'TT':('time',range(len(NC.time)))})
NC=NC.swap_dims({"time":"TT"})
NC_s=NC.polyfit(dim='TT',deg=1,skipna=True)
fit = xr.polyval(NC.TT, NC_s.polyfit_coefficients)
Coef=NC_s.polyfit_coefficients[0]
Coef_var=Coef.values*12 # (m/year)
adt_dt=NC-fit
adt_dt=adt_dt.swap_dims({"TT":"time"})

adt_dt_2Y=adt_dt.rolling(time=24,center=True).mean()

adt_dt_2Y_=adt_dt_2Y.values

# adt_dt_2Y[10].plot()


NAME=[str(i)[:7] for i in adt_dt_2Y.time.values]

lon_m,lat_m=np.meshgrid(adt_dt_2Y.longitude.values,adt_dt_2Y.latitude.values)


def Plot_SO_Merc3(lonA,latA,MyDATA,t_name,CMAP,Mylim,w_path,save_name,fig_bool=False):
    
    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(16.5,6),
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
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True,dpi=200)
        plt.savefig(w_path+'/'+save_name,bbox_inches='tight',dpi=200)
    plt.show()

for i,j in zip(adt_dt_2Y_,NAME):
    t_name='adt dt2Y '+j
    s_name='adt_dt2Y_'+j.replace('-','_')

    Plot_SO_Merc3(lon_m,lat_m,i,t_name,plt.get_cmap('RdBu_r',15),[-0.1,0.1],\
                  'D:/HEAT/tmp/',s_name,fig_bool=False)




NC=xr.open_mfdataset(pth+'*.nc').loc[dict(latitude=slice(-40,-20),longitude=slice(180,240),nv=0)]\
    .adt.mean(dim=['latitude','longitude'])

NC.rolling(time=24,center=True).mean().plot()
adt_dt_2Y.plot()




