# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:00:18 2023

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
AVGS_path2='G:/Models/TK0525ED_CLM/'

save_path='D:/OneDrive/JNUpack/Mapping/Samples_figs/'

lat_rng=[-80,-24]

AVGS1=[AVGS_path1+i for i in os.listdir(AVGS_path1) if i.endswith('.nc')]
AVGS2=[AVGS_path2+i for i in os.listdir(AVGS_path2) if i.endswith('.nc')]

Sample1=xr.open_dataset(AVGS1[0])
Sample1.st_ocean.values
SALT1=xr.open_mfdataset(AVGS1).temp.loc[dict(st_ocean=Sample1.st_ocean.values[0],time=slice('1990-01','2016-12'))].groupby('time.season').mean()


Sample2=xr.open_dataset(AVGS2[0])
Sample2.s_rho.values
SALT2=xr.open_mfdataset(AVGS2).temp.loc[dict(s_rho=Sample2.s_rho.values[-1],ocean_time=slice('1990-01','2016-12'))].groupby('ocean_time.season').mean()

lon1=Sample1.xt_ocean.values
lat1=Sample1.yt_ocean.values

lon1_m,lat1_m=np.meshgrid(lon1,lat1)

SALT2_DJF=SALT2.loc[dict(season='DJF')].squeeze().values
SALT1_DJF=SALT1.loc[dict(season='DJF')].squeeze().values

SALT2_JJA=SALT2.loc[dict(season='JJA')].squeeze().values
SALT1_JJA=SALT1.loc[dict(season='JJA')].squeeze().values


SALT2_DJF_re_=griddata((SALT2.lon_rho.values.flatten(),SALT2.lat_rho.values.flatten()),SALT2_DJF.flatten(),
              (lon1_m.flatten(),lat1_m.flatten()),
           method='linear',fill_value=np.nan)
SALT2_DJF_re = SALT2_DJF_re_.reshape(lon1_m.shape)

SALT2_JJA_re_=griddata((SALT2.lon_rho.values.flatten(),SALT2.lat_rho.values.flatten()),SALT2_JJA.flatten(),
              (lon1_m.flatten(),lat1_m.flatten()),
           method='linear',fill_value=np.nan)
SALT2_JJA_re = SALT2_JJA_re_.reshape(lon1_m.shape)

diff_DJF=SALT2_DJF_re-SALT1_DJF
diff_JJA=SALT2_JJA_re-SALT1_JJA



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

import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt
import cmocean

data_lim=[-1,24]

NN=15
#My_levels=np.linspace(data_lim[0],data_lim[-1],NN+1,endpoint=True)

My_levels=np.arange(data_lim[0],data_lim[-1]+0.2/2,1)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
CMAP = ListedColormap(plt.get_cmap('jet')(np.linspace(0, 1, len(My_levels),endpoint=True)))
CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, len(My_levels),endpoint=True)))

mm.Plot_SO_Spherical2(lon1_m,lat1_m,diff_JJA,My_levels,CMAP,data_lim,save_path,'Spherical_zeta_sample',fig_bool=False)


mm.Plot_SO_Spherical2(SALT2.lon_rho.values,SALT2.lat_rho.values,SALT2_JJA,My_levels,CMAP,data_lim,save_path,'Spherical_zeta_sample',fig_bool=False)


mm.Plot_SO_Merc2(lon_rho,lat_rho,ZETA,My_levels,CMAP,data_lim,save_path,'Merc_zeta_sample',fig_bool=False)

