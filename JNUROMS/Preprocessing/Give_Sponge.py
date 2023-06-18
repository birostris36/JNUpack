# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 12:37:11 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


Grd_name='G:/MODEL_DATA/Grd/Grd_Q0_Rtopo30S_Smooth_sponge.nc'

Factor=10
Nfilter=24

ncG=Dataset(Grd_name,'r+')

lon_rho,lat_rho=ncG['lon_rho'][:],ncG['lat_rho'][:]
mask_rho=ncG['mask_rho']
at,on=mask_rho.shape
plt.pcolormesh(mask_rho)
Sponge_co=np.ones_like(mask_rho)
Sponge_co[-1,:]=Factor
Fac_space=np.linspace(1,Factor,Nfilter)
for i in range(on):
    Sponge_co[-Nfilter-1:-1,i]=Fac_space

VISC=ncG.createVariable('visc_factor', 'f4',('eta_rho','xi_rho'))
VISC.long_name='horizontal viscosity sponge factor'
VISC.valid_min=0.0
VISC.coordinates='lon_rho lat_rho'

DIFF=ncG.createVariable('diff_factor', 'f4',('eta_rho','xi_rho'))
DIFF.long_name='horizontal diffusivity sponge factor'
DIFF.valid_min=0.0
DIFF.coordinates='lon_rho lat_rho'

ncG['visc_factor'][:]=Sponge_co
ncG['diff_factor'][:]=Sponge_co
ncG.close()




# ncG=Dataset(Grd_name,'r')

# ncG['diff_factor']
# tmp['diff_factor']


# plt.pcolor(ncG['diff_factor'][:])
# plt.pcolor(tmp['diff_factor'][:])


# tmp=Dataset('D:/Working_hub/OneDrive/base142/Warehouse01/Grd_SO_05d_sponge.nc')

# len(tmp.variables.keys())
# len(ncG.variables.keys())









