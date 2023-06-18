# -*- coding: utf-8 -*-
"""
Created on Tue May  2 00:08:32 2023

@author: shjo9
"""

import sys
sys.path.append('D:/OneDrive/JNUpack/')
sys.path.append('D:/OneDrive/JNUpack/JNUROMS/')
from Tools.Manta_WindStress import ra_windstrcurl,ra_windstr
import xarray as xr
import numpy as np
from netCDF4 import Dataset, num2date, date2num,MFDataset
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Mapping.Tools import d_modules as mm
from scipy.interpolate import interp2d, griddata
import matplotlib.ticker as mticker
import numpy.ma as ma
from tqdm import tqdm
from matplotlib import ticker
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pandas as pd
import cmocean

fig_bool=0
w_path_sig='D:/WIND_/monthly/'
save_name='UV_WSC_M'

Mytime=['2000-01','2000-12']


Wind=xr.open_dataset('G:/ERA5_monthly_wind_1980_2021.nc').loc[dict(wind_time=slice(Mytime[0],Mytime[-1]),\
                                                                   lat=slice(-80,-10))]

YM_U=Wind.Uwind
YM_V=Wind.Vwind

U_lon=Wind.lon.values
U_lat=Wind.lat.values

U_lon_m,U_lat_m=np.meshgrid(U_lon,U_lat)

lat_new = np.arange(U_lat[0],U_lat[-1], 5)
lon_new = np.arange(U_lon[0],U_lon[-1], 5)
lon_new_m,lat_new_m=np.meshgrid(lon_new,lat_new)

U_re,V_re,WSC=[],[],[]

for i,j in tqdm(zip(YM_U,YM_V)):
    U_=griddata((U_lon_m.flatten(),U_lat_m.flatten()),i.values.flatten(),
                  (lon_new_m.flatten(),lat_new_m.flatten()),
               method='linear',fill_value=np.nan)
    U_re_ = U_.reshape(lon_new_m.shape)

    V_=griddata((U_lon_m.flatten(),U_lat_m.flatten()),j.values.flatten(),
                  (lon_new_m.flatten(),lat_new_m.flatten()),
               method='linear',fill_value=np.nan)
    V_re_ = V_.reshape(lon_new_m.shape)

    WSC_=ra_windstrcurl(U_lat,U_lon,i.values,j.values)
    
    U_re.append(U_re_); V_re.append(V_re_); WSC.append(WSC_);
U_re=np.array(U_re); V_re=np.array(V_re); WSC=np.array(WSC); 


# =============================================================================
# CMAP
# =============================================================================
NN=16

data_lim=[-10**-7,10**-7]

My_levels=np.arange(data_lim[0],data_lim[-1]+10**-8/2,10**-8)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))

CMAP=ListedColormap(cmocean.cm.curl(np.linspace(0, 1, len(My_levels)+1,endpoint=True)))

My_date=pd.date_range(Mytime[0],Mytime[-1],freq='m').strftime('%Y-%m')
Label_size=12
# =============================================================================
# # Figure 1 
# =============================================================================
Spheric=ccrs.SouthPolarStereo(central_longitude=180.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=180.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)

for i, j, k, tt in zip(U_re, V_re, WSC, My_date):
    
    Title_name='Date: '+tt+'   '

    fig, ax = plt.subplots(1, 1, figsize=(8,2),
                       subplot_kw={'projection': PC},dpi=200)
    ax.set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.')
        
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 14}
        
    # gl.xlocator = mticker.FixedLocator( np.arange(0,359,20) )
    # gl.xformatter = LONGITUDE_FORMATTER
    
    # ax.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree(180))
    
    ax.add_feature(cf.COASTLINE.with_scale("50m"), lw=1,zorder=110)
    ax.add_feature(cf.LAND.with_scale("50m"),color=[.75,.75,.75],zorder=100)
    
    # q1 = plt.quiver(lon_new_m,lat_new_m,i,j,
    #     scale=300,headwidth=8.,headaxislength=10,headlength=13,color='k',
    #     minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
    #     pivot='mid',angles='xy')
    
    q1 = plt.quiver(lon_new_m,lat_new_m,i,j,
        scale=300, width=0.0025, headwidth=4.5, headaxislength=4,color='k',
        edgecolor='k',alpha=.75,transform=ccrs.PlateCarree(),zorder=1,minlength=1,
        pivot='mid',angles='xy')
    
    plt.quiverkey(q1,-168.,-4.,10,"$1\cdot10^{1} m/s$",coordinates='data',color='r',
                labelpos='E',alpha=1.,labelcolor='k',fontproperties={'size':12},
                labelsep=0.13,transform=ccrs.PlateCarree(),zorder=3)
    
    C1=plt.contourf(U_lon_m,U_lat_m,k,cmap=CMAP,transform=ccrs.PlateCarree(),zorder=0,levels=My_levels)

    plt.clim(data_lim[0],data_lim[-1])
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
    ax.set_extent([0, 359, -80, -10], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)
    
    fig.add_axes(ax_cb)
    cb=plt.colorbar(C1,extend='both',pad=0.01,cax=ax_cb)
    if fig_bool:
        plt.savefig(w_path_sig+'ppt/'+save_name+'_'+tt.replace('-',''),
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path_sig+save_name+'_'+tt.replace('-',''),bbox_inches='tight')
    plt.show()
    

