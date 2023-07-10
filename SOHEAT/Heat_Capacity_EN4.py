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

nc_pth='G:/EN4/'

NC = [nc_pth+i for i in os.listdir(nc_pth) if i.endswith('.nc')]

EN4=xr.open_mfdataset(NC).loc[dict(depth=slice(0,2000))]

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

OHC[0].plot(cmap=plt.get_cmap('jet',15),vmin=-1*10**10,vmax=12*10**10)













