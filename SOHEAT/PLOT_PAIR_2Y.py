


# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:47:00 2023

@author: shjo9
"""
import matplotlib
matplotlib.use('agg')
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.path as mpath
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from tqdm import tqdm


pth='G:/ERA5_monthly_85/'

NC=xr.open_dataset(pth+'ERA5_monthly_pair_1980_2021.nc')\
    .loc[dict(lat=slice(-65,40))].Pair

NC_a = NC - NC.mean(dim='pair_time')


NC_2Y=NC_a.rolling(pair_time=24,center=True).mean()

NC_2Y_=NC_2Y.values

NAME=[str(i)[:7] for i in NC_2Y.pair_time.values]

lon_m,lat_m=np.meshgrid(NC_2Y.lon.values,NC_2Y.lat.values)

def Plot_SO_Merc3(lonA,latA,MyDATA,t_name,CMAP,Mylim,My_levels,w_path,save_name,fig_bool=False):
    
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

    M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    # M=plt.pcolormesh(lonA, latA, MyDATA,
    #               transform=PC,cmap=CMAP)
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
        plt.savefig(w_path+'/'+save_name,dpi=200)
    plt.show()

OHC_cmap=plt.get_cmap('RdYlBu_r',15)

NN=15
OHC_lim=[-5,5]
# OHC_levels=np.arange(OHC_lim[0],OHC_lim[-1])
OHC_levels=np.linspace(OHC_lim[0],OHC_lim[-1],NN)

OHC_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')
                           (np.linspace(0, 1, NN,endpoint=True)) )


for i,j in zip(NC_2Y_,NAME):
    t_name='Air pressure 10m anomaly 2Y '+j
    s_name='pair_a2Y_'+j.replace('-','_')

    i[i>OHC_lim[-1]]=OHC_lim[-1]
    i[i<OHC_lim[0]]=OHC_lim[0]
    
    Plot_SO_Merc3(lon_m,lat_m,i,t_name,plt.get_cmap('RdBu_r',15),OHC_lim,OHC_levels,\
                  'D:/HEAT/Pair/',s_name,fig_bool=True)





'''
NC=xr.open_dataset(pth+'EN4_OHC_GLOBAL_c14_700m_1980_2023.nc')\
    .loc[dict(lat=slice(-40,-20),lon=slice(180,240))].OHC.mean(dim=['lat','lon'])

NC.rolling(time=24,center=True).mean().plot()

'''









