# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 10:16:46 2022

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""
import sys
sys.path.append('/data4/base158/Factory/ERA5_daily_shjo/')
from netCDF4 import Dataset
import xarray as xr
import numpy as np
import create_tools as ct

r_path='/data4/base158/Warehouse01/'


if 1:
    nc_=xr.open_dataset(r_path+'ERA5_1993_hourly.nc')
    nc_s=nc_.loc[dict(latitude=slice(10,-80))]
    Daily=nc_s.groupby('time.dayofyear').mean()
    Daily= Daily.rename({'dayofyear':'time'})
    # Daily.to_netcdf('/data4/base158/Warehouse01/daily/EAR5_daily_1993.nc')


# A=nc.sst[1,:,:]-273.15

LON=Daily.longitude.values
LAT=np.flip(Daily.latitude.values,axis=0)
TIME=Daily.time.values

srf_values = np.flip(Daily.ssr.values,axis=1)/3600
lrf_values = np.flip(Daily.str.values,axis=1)/3600
lrf_down_values=np.flip(Daily.strd.values,axis=1)/3600
# latent_values=Daily.slhf.values/86400
u_values=np.flip(Daily.u10.values,axis=1)
v_values=np.flip(Daily.v10.values,axis=1)

cloud_values=np.flip(Daily.tcc.values,axis=1)
rain_values=np.flip(Daily.tp.values,axis=1)*1000/3600

sst_values=np.flip(Daily.sst.values,axis=1) - 273.15
sst_values[sst_values!=sst_values]=0
T2_values=np.flip(Daily.t2m.values,axis=1) - 273.15
pair_values=np.flip(Daily.msl.values,axis=1) / 100

D2_values=np.flip(Daily.d2m.values,axis=1) - 273.15
TK_values=np.flip(Daily.t2m.values,axis=1)

qair_values = 100 * ( np.exp( (17.625 * D2_values ) /  (243.04 + D2_values ) ) / 
                     np.exp( ( 17.625 * T2_values ) / (243.04 + T2_values ) ) )

dqdsst_values = ct.get_dqdsst(sst_values, TK_values, 1.2, (u_values**2+v_values**2)**(1/2), qair_values/1000)

TIME=np.arange(1,366)


ct.create_wind_nc(r_path+'daily/'+'ERA5_daily_wind_1993.nc',LON,LAT,TIME,365,u_values,v_values)

ct.create_srf_nc(r_path+'daily/'+'ERA5_daily_srf_1993.nc',LON,LAT,TIME,365,srf_values)
ct.create_lwrad_nc(r_path+'daily/'+'ERA5_daily_lrf_1993.nc',LON,LAT,TIME,365,lrf_values)
ct.create_lwrad_down_nc(r_path+'daily/'+'ERA5_daily_lrf_down_1993.nc',LON,LAT,TIME,365,lrf_down_values)


ct.create_sst_nc(r_path+'daily/'+'ERA5_daily_sst_1993.nc',LON,LAT,TIME,365,sst_values)
ct.create_pair_nc(r_path+'daily/'+'ERA5_daily_pair_1993.nc',LON,LAT,TIME,365,pair_values)
ct.create_tair_nc(r_path+'daily/'+'ERA5_daily_tair_1993.nc',LON,LAT,TIME,365,T2_values)
ct.create_qair_nc(r_path+'daily/'+'ERA5_daily_qair_1993.nc',LON,LAT,TIME,365,qair_values)
ct.create_dqdsst_nc(r_path+'daily/'+'ERA5_daily_dqdsst_1993.nc',LON,LAT,TIME,365,dqdsst_values)
ct.create_rain_nc(r_path+'daily/'+'ERA5_daily_rain_1993.nc',LON,LAT,TIME,365,rain_values)
ct.create_cloud_nc(r_path+'daily/'+'ERA5_daily_cloud_1993.nc',LON,LAT,TIME,365,cloud_values)















