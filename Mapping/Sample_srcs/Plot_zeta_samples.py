#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 13:01:33 2023

@author: shjo
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

AVGS_path='G:/Models/TK0525ED_CLM/'
save_path='D:/OneDrive/JNUpack/Mapping/Samples_figs/'

lat_rng=[-80,-24]

AVGS=[AVGS_path+i for i in os.listdir(AVGS_path) if i.endswith('.nc')]

Sample_Data=Dataset(AVGS[0])
lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])

ZETA=Sample_Data['zeta'][0,lat_co,:]

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

data_lim=[-2,1.5]

NN=16
#My_levels=np.linspace(data_lim[0],data_lim[-1],NN+1,endpoint=True)

My_levels=np.arange(data_lim[0],data_lim[-1]+0.2/2,0.2)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 0.9, len(My_levels)+1,endpoint=True)))

mm.Plot_SO_Spherical2(lon_rho,lat_rho,ZETA,'22',My_levels,zeta_CMAP,data_lim,save_path,'Spherical_zeta_sample',fig_bool=False)
mm.Plot_SO_Merc2(lon_rho,lat_rho,ZETA,My_levels,zeta_CMAP,data_lim,save_path,'Merc_zeta_sample',fig_bool=False)

