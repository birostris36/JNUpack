# -*- coding: utf-8 -*-
"""
Created on Tue May  2 13:36:19 2023

@author: shjo9
"""


import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date, date2num,MFDataset
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
sys.path.append('D:/OneDrive/JNUpack/')
from Mapping.Tools import d_modules as mm

AVGS_path1='G:/SODA/'
AVGS_path1='G:/Models/TK0525ED_CLM/'

save_path='D:/OneDrive/JNUpack/Mapping/Samples_figs/'

lat_rng=[-80,-24]

AVGS1=[AVGS_path1+i for i in os.listdir(AVGS_path1) if i.endswith('.nc')]
# AVGS2=[AVGS_path2+i for i in os.listdir(AVGS_path2) if i.endswith('.nc')]

Sample1=xr.open_dataset(AVGS1[0])
# Sample1.st_ocean.values
AICE=xr.open_mfdataset(AVGS1).zeta.loc[dict(ocean_time=slice('1980-02','2016-12'))]
# A_CMm=AICE.groupby('ocean_time.month').mean()

# Sample2=xr.open_dataset(AVGS2[0])
# Sample2.s_rho.values
# SALT2=xr.open_mfdataset(AVGS2).salt.loc[dict(s_rho=Sample2.s_rho.values[-1],ocean_time=slice('2000-01','2016-01'))].mean(dim='ocean_time')


# SALT2_re_=griddata((SALT2.lon_rho.values.flatten(),SALT2.lat_rho.values.flatten()),SALT2.squeeze().values.flatten(),
#               (lon1_m.flatten(),lat1_m.flatten()),
#            method='linear',fill_value=np.nan)
# SALT2_re = SALT2_re_.reshape(lon1_m.shape)


# diff_salt=SALT2_re-SALT1.squeeze().values

plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "regular"
plt.rcParams['axes.linewidth'] = 1
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True

SIZE=(11,4)
FS=14

# =============================================================================
# Cartopy Mercator
# =============================================================================

# Scalar={'u':U,'lon':lon_rho,'lat':lat_rho}
import pandas as pd
My_date=pd.date_range('1980-01','2018-01',freq='Y').strftime('%Y')

Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)

fig_bool=1
w_path='D:/OneDrive/EMER/ZETA/'
save_name='ZETA'

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature
from mpl_toolkits.axes_grid1 import make_axes_locatable


import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt
import cmocean

data_lim=[-2.0,1.5]
Label_size=28
NN=15
#My_levels=np.linspace(data_lim[0],data_lim[-1],NN+1,endpoint=True)
t=0
My_levels=np.arange(data_lim[0],data_lim[-1]+0.2/2,0.1)
CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
# CMAP = ListedColormap(plt.get_cmap('seismic')(np.linspace(0, 1, len(My_levels),endpoint=True)))
for i,tt in zip(AICE.values,My_date):
    # mm.Plot_SO_Spherical2(AICE.lon_rho.values,AICE.lat_rho.values,i,My_levels,CMAP,data_lim,save_path,'Spherical_zeta_sample',fig_bool=False)


    
    # Now we will create axes object having specific projection 
    t+=1
    Title_name='Years: '+tt+f' (+{t-1:02d}~{t:02d})'


    # Now we will create axes object having specific projection 
    
    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                       subplot_kw={'projection': Spheric})
    ax.set_title(Title_name,fontdict={'fontsize':Label_size,'fontweight':'regular'},)

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.xlabels_top,gl.ylabels_right = True,True
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    
    # M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    M=plt.pcolormesh(AICE.lon_rho.values, AICE.lat_rho.values, i,
                  transform=PC,cmap=CMAP)
    plt.clim(data_lim[0],data_lim[-1])
    
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
        plt.savefig(w_path+'ppt/'+save_name+'_'+tt,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+save_name+'_'+tt)
    plt.show()
    
    



    













