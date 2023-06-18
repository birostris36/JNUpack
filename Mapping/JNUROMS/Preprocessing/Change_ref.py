# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:41:48 2023

@author: shjo9
"""

from netCDF4 import Dataset, date2num, num2date
import datetime as dt
A=Dataset('G:/MODEL_DATA/TK_Inputs/roms_bdy_1970-1992.nc')['bry_time'][:]
B=Dataset('G:/MODEL_DATA/TK_Inputs/roms_bdy_1993-2019.nc')['bry_time'][:]


nc_pth='G:/MODEL_DATA/TK_Inputs/roms_bdy_1970-1992_ref70.nc'

nc01=Dataset(nc_pth,'r')


My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=nc01['bry_time'][:]
TIME_UNIT='days since 1950-1-1 00:00:00'
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)

# =============================================================================
tmp_dif=date2num(dt.datetime(1993,1,1),TIME_UNIT)-date2num(dt.datetime(1993,1,1),My_time_ref)
Bry_time_time=num2date(OGCM_TIMES-tmp_dif,My_time_ref)
Bry_time_num=OGCM_TIMES-tmp_dif

nc01['bry_time'][:]=Bry_time_num

nc01.close()













import xarray as xr


A=Dataset('G:/MODEL_DATA/TK_Inputs/roms_bdy_1970-1992.nc')['bry_time'][:]
B=Dataset('G:/MODEL_DATA/TK_Inputs/roms_bdy_1993-2019.nc')['bry_time'][:]




















