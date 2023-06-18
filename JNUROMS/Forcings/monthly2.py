# -*- coding: utf-8 -*-
"""
Created on Tue May  2 01:29:51 2023

@author: shjo9
"""

import sys
sys.path.append('D:/Onedrive/base142/Factory/ERA5_daily_shjo/')
from netCDF4 import Dataset, date2num, num2date
import xarray as xr
import numpy as np
import create_tools2 as ct

r_path='G:/ERA5_low/'
w_path='G:/ERA5_monthly_85/'

if 1:
    nc_=xr.open_dataset(r_path+'ERA5_monthly_msl_1979_2022.nc',decode_times=False)
    nc_s=nc_.loc[dict(latitude=slice(10,-85))]
    # Daily=nc_s.groupby('time.dayofyear').mean()
    # Daily= Daily.rename({'dayofyear':'time'})
    # Daily.to_netcdf('/data4/base158/Warehouse01/daily/EAR5_daily_1993.nc')


# A=nc.sst[1,:,:]-273.15

LON=nc_s.longitude.values
LAT=np.flip(nc_s.latitude.values,axis=0)
TIME=nc_s.time.values/24

import datetime as dt 


tmp_time_var='time'
t_rng=['1979-01','2022-12']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=nc_s[tmp_time_var]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),28)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),'days since 1900-1-1 00:00:00')-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
Frc_time_time=num2date(OGCM_TIMES[TIMES_co]/24-tmp_dif,My_time_ref)
Frc_time_num=OGCM_TIMES[TIMES_co]/24-tmp_dif



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



# =============================================================================
# 
# =============================================================================
ct.create_wind_nc(w_path+'ERA5_monthly_wind_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,u_values,v_values)

ct.create_srf_nc(w_path+'ERA5_monthly_srf_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,srf_values)
ct.create_lwrad_nc(w_path+'ERA5_monthly_lrf_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,lrf_values)
ct.create_lwrad_down_nc(w_path+'ERA5_monthly_lrf_down_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,lrf_down_values)

ct.create_sst_nc(w_path+'ERA5_monthly_sst_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,sst_values)
ct.create_pair_nc(w_path+'ERA5_monthly_pair_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,pair_values)
ct.create_tair_nc(w_path+'ERA5_monthly_tair_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,T2_values)
ct.create_qair_nc(w_path+'ERA5_monthly_qair_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,qair_values)
# ct.create_dqdsst_nc(w_path+'ERA5_monthly_dqsst_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,dqdsst_values)
ct.create_rain_nc(w_path+'ERA5_monthly_rain_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,rain_values)
ct.create_cloud_nc(w_path+'ERA5_monthly_cloud_1980_2021.nc',LON,LAT,Frc_time_num,My_time_ref,cloud_values)















