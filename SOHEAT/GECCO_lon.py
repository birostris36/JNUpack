# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 12:39:02 2023

@author: shjo9
"""
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import xarray as xr

w_pth='E:/_data/GECCO/'

DATA=Dataset('E:/_data/GECCO/GECCO_SALT_194801_201812_ori.nc')


data=DATA['salt'][:]
lon=DATA['lon'][:]

posi_co=np.where(lon>=0)[0]
nega_co=np.where(lon<=0)[0] 

lon_new=np.zeros_like(lon)

lon_new[:180]=lon[posi_co]
lon_new[180:]=lon[nega_co]+360

data_new=np.zeros_like(data)

data_new[:,:,:,:180]=data[:,:,:,posi_co]
data_new[:,:,:,180:]=data[:,:,:,nega_co]

# lon_m,lat_m=np.meshgrid(lon_new,DATA['lat'][:])
# plt.pcolor(lon_m,lat_m,data_new[0,0,:,:])

DATA_xr=xr.open_dataset('E:/_data/GECCO/GECCO_SALT_194801_201812_ori.nc')

DATA_xr=DATA_xr.assign({'salt': (('time','Depth','lat','lon'),data_new) })

DATA_xr=DATA_xr.assign_coords( {'lon':('lon',lon_new) } )

DATA_xr.to_netcdf(w_pth+'GECCO_SALT_194801_201812.nc')


# A=xr.open_dataset(w_pth+'GECCO_SALT_194801_201812.nc')

