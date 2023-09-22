from scipy import io
from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import pandas as pd
from numpy import ma
import os

npth='D:/ERA5_low/monthly_full/ERA5_monthly.nc'
wnpth='J:/Reanalysis/myERA5_198001_202212_sst.nc'

Sample=Dataset(npth)
LAT,LON=np.flipud(Sample['latitude'][:]),Sample['longitude'][:]

# Process Times
tmp_time_var='time'
t_rng=['1980-01','2022-12']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=Sample[tmp_time_var]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)/24-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
my_time=num2date(OGCM_TIMES[TIMES_co]/24-tmp_dif,My_time_ref)
my_time_num=OGCM_TIMES[TIMES_co]/24-tmp_dif

SST=np.flip(Dataset(npth)['sst'][TIMES_co,0,:,:],axis=1)-273.15
SST.data[SST.data<-1000]=np.nan

def mySST(nc_save_name,LON,LAT,TIME,Ref_time,values1):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My ERA5 SST data '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('time', np.float64, ('time',))
    time.units=Ref_time
    time.field='time, scalar, series'
    
    DATA1 = ncfile.createVariable('sst',np.float64,('time','lat','lon'),compression='zlib') #
    DATA1.units = 'degree_C' 
    DATA1.long_name = 'ERA5 SST' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME 
     
    DATA1[:] = values1

    ncfile.close()

mySST(wnpth,LON,LAT,my_time_num,My_time_ref,SST)







