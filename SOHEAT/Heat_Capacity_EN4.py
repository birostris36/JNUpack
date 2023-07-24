# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 14:54:44 2023

@author: shjo9
"""

import matplotlib.pyplot as plt
import gsw.density as gsw_d
import gsw.conversions as gsw_c
from gsw._wrapped_ufuncs import cp_t_exact
import xarray as xr
import os
import numpy as np

nc_pth='G:/EN4_22_anal/'

NC = [nc_pth+i for i in os.listdir(nc_pth) if i.endswith('.nc')]

EN4=xr.open_mfdataset(NC).loc[dict(depth=slice(0,700),time=slice('1980-01','2023-12'),lat=slice(-80,-10))]

LON,LAT,DEPTH=EN4.lon,EN4.lat,EN4.depth
db=EN4.depth_bnds

TEMP,SALT=EN4.temperature-273.15, EN4.salinity

#Practical Salinity --> Absolute salinity
SA = gsw_c.SA_from_SP(SALT,DEPTH,LON,LAT) # [SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat)

# Potential temp --> Conservative temp
CT=gsw_c.CT_from_pt(SA,TEMP) #CT = gsw_CT_from_pt(SA,pt)

rho = gsw_d.rho(SA,CT,DEPTH)

# Potential temperature --> In-situ temperature
t=gsw_c.t_from_CT(SA,CT,DEPTH)

# Heat capacity
CP = cp_t_exact(SA,t,DEPTH)

dz = db[:,:,1]-db[:,:,0]

OHC_=CP*CT*rho*dz

# Integrates from ref depth (2000m) 
OHC=OHC_.sum(dim='depth',skipna=False)

OHC=OHC.rename('OHC')

OHC.to_netcdf('D:/HEAT/EN4_OHC_SO_c14_700m2000m_1980_2023.nc',format='netcdf4')


# OHC[-1].plot(cmap=plt.get_1cmap('jet',15),vmin=-1*10**10,vmax=12*10**10,)

# xr.open_dataset('D:/HEAT/EN4_OHC_GLOBAL_c14_700m_1980_2023.nc')







