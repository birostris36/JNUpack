# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 09:47:22 2023

@author: shjo9
"""

import matplotlib.pyplot as plt
import gsw.density as gsw_d
import gsw.conversions as gsw_c
from gsw._wrapped_ufuncs import cp_t_exact
import xarray as xr
import os
import numpy as np

nc_pth='E:/_data/GECCO/'

TEMP_nc=xr.open_dataset(nc_pth+'GECCO_TEMP_194801_201812.nc')\
    .loc[dict(Depth=slice(0,700),time=slice('1980-01','2023-12'),lat=slice(-80,-10),lon=slice(300,360))]
SALT_nc=xr.open_dataset(nc_pth+'GECCO_SALT_194801_201812.nc')\
    .loc[dict(Depth=slice(0,700),time=slice('1980-01','2023-12'),lat=slice(-80,-10),lon=slice(300,360))]

LON,LAT,DEPTH=TEMP_nc.lon,TEMP_nc.lat,TEMP_nc.temp.Depth

TEMP,SA=TEMP_nc.temp, SALT_nc.salt

#Practical Salinity --> Absolute salinity
# SA = gsw_c.SA_from_SP(SALT,DEPTH,LON,LAT) # [SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat)

# Potential temp --> Conservative temp
CT=gsw_c.CT_from_pt(SA,TEMP) #CT = gsw_CT_from_pt(SA,pt)

rho = gsw_d.rho(SA,CT,DEPTH)

# Potential temperature --> In-situ temperature
t=gsw_c.t_from_CT(SA,CT,DEPTH)

# Heat capacity
CP = cp_t_exact(SA,t,DEPTH)

DZ=np.concatenate( ([6], np.diff(DEPTH)) ,axis=0)

dz=np.tile(np.tile(np.tile(DZ.reshape([DEPTH.shape[0],1]), 70 )\
                   .reshape([DEPTH.shape[0],70,-1]), 360), (TEMP.shape[0],1,1,1) )

OHC=CP*CT*rho

# Integrates from ref depth (2000m) 
# OHC=OHC_.sum(dim='Depth',skipna=False)
OHC=OHC.mean(dim='lon')
OHC=OHC.rename('OHC')

OHC.to_netcdf('E:/_data/MyOHC/V_GECCO_OHC_SO_c14_ATL_700m_1980_2023.nc',format='netcdf4')


# OHC[-1].plot(cmap=plt.get_1cmap('jet',15),vmin=-1*10**10,vmax=12*10**10,)

# xr.open_dataset('D:/HEAT/EN4_OHC_GLOBAL_c14_700m_1980_2023.nc')



