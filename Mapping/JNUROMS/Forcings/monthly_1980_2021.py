# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:03:43 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""


import sys
sys.path.append('D:/Working_hub/Onedrive/base142/Factory/ERA5_daily_shjo/')
from netCDF4 import Dataset, date2num, num2date
import xarray as xr
import numpy as np
import create_tools2 as ct

r_path='C:/Users/shjo9/Downloads/'
w_path='D:/Working_hub/DATA_dep/Frc_tmp/'

if 1:
    nc_=xr.open_dataset(r_path+'ERA5_monthly_1980_2022_raw.nc',decode_times=False)
    nc_s=nc_.loc[dict(latitude=slice(10,-80))]
    # Daily=nc_s.groupby('time.dayofyear').mean()
    # Daily= Daily.rename({'dayofyear':'time'})
    # Daily.to_netcdf('/data4/base158/Warehouse01/daily/EAR5_daily_1993.nc')


# A=nc.sst[1,:,:]-273.15

LON=nc_s.longitude.values
LAT=np.flip(nc_s.latitude.values,axis=0)
TIME=nc_s.time.values/24-18262

# import datetime as dt 
# date2num(dt.datetime(1950, 1, 1),'days since 1900-01-01')
# num2date(TIME-18262,'days since 1950-01-01')
# Dataset(r_path+'ERA5_monthly_1980_2022_raw.nc')['time']


srf_values = np.flip(nc_s.ssr[:,:,:].values,axis=1)/86400
lrf_values = np.flip(nc_s.str[:,:,:].values,axis=1)/86400
lrf_down_values=np.flip(nc_s.strd[:,:,:].values,axis=1)/86400
# latent_values=Daily.slhf.values/86400
u_values=np.flip(nc_s.u10[:,:,:].values,axis=1)
v_values=np.flip(nc_s.v10[:,:,:].values,axis=1)

cloud_values=np.flip(nc_s.tcc[:,:,:].values,axis=1)
rain_values=np.flip(nc_s.tp[:,:,:].values,axis=1)*1000/86400

sst_values=np.flip(nc_s.sst[:,:,:].values,axis=1) - 273.15
sst_values[sst_values!=sst_values]=0
T2_values=np.flip(nc_s.t2m[:,:,:].values,axis=1) - 273.15
pair_values=np.flip(nc_s.msl[:,:,:].values,axis=1) / 100

D2_values=np.flip(nc_s.d2m[:,:,:].values,axis=1) - 273.15
TK_values=np.flip(nc_s.t2m[:,:,:].values,axis=1)

qair_values = 100 * ( np.exp( (17.625 * D2_values ) /  (243.04 + D2_values ) ) / 
                     np.exp( ( 17.625 * T2_values ) / (243.04 + T2_values ) ) )

dqdsst_values = ct.get_dqdsst(sst_values, TK_values, 1.2, (u_values**2+v_values**2)**(1/2), qair_values/1000)

# TIME=np.arange(1,366)
Ref_time='days since 1950-01-01'

ct.create_wind_nc(w_path+'ERA5_monthly_wind_1980_2021.nc',LON,LAT,TIME,Ref_time,u_values,v_values)

ct.create_srf_nc(w_path+'ERA5_monthly_srf_1980_2021.nc',LON,LAT,TIME,Ref_time,srf_values)
ct.create_lwrad_nc(w_path+'ERA5_monthly_lrf_1980_2021.nc',LON,LAT,TIME,Ref_time,lrf_values)
ct.create_lwrad_down_nc(w_path+'ERA5_monthly_lrf_down_1980_2021.nc',LON,LAT,TIME,Ref_time,lrf_down_values)

ct.create_sst_nc(w_path+'ERA5_monthly_sst_1980_2021.nc',LON,LAT,TIME,Ref_time,sst_values)
ct.create_pair_nc(w_path+'ERA5_monthly_pair_1980_2021.nc',LON,LAT,TIME,Ref_time,pair_values)
ct.create_tair_nc(w_path+'ERA5_monthly_tair_1980_2021.nc',LON,LAT,TIME,Ref_time,T2_values)
ct.create_qair_nc(w_path+'ERA5_monthly_qair_1980_2021.nc',LON,LAT,TIME,Ref_time,qair_values)
ct.create_dqdsst_nc(w_path+'ERA5_monthly_dqsst_1980_2021.nc',LON,LAT,TIME,Ref_time,dqdsst_values)
ct.create_rain_nc(w_path+'ERA5_monthly_rain_1980_2021.nc',LON,LAT,TIME,Ref_time,rain_values)
ct.create_cloud_nc(w_path+'ERA5_monthly_cloud_1980_2021.nc',LON,LAT,TIME,Ref_time,cloud_values)















