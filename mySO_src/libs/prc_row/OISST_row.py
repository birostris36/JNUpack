from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import os
from numpy import ma

npth='J:/Monyhly_oisst.nc'
wnpth='J:/Reanalysis/myOISST_198001_202212_sst.nc'

Sample=Dataset(npth)
LAT,LON=Sample['lat'][:],Sample['lon'][:]

# Process Times
tmp_time_var='time'
t_rng=['1980-01','2022-12']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=Sample['time'][:]
TIME_UNIT=Sample['time'].units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
my_time=num2date(OGCM_TIMES[TIMES_co]-tmp_dif,My_time_ref)
my_time_num=OGCM_TIMES[TIMES_co]-tmp_dif

OISST_ice = Dataset(npth)['ice'][TIMES_co,:,:].squeeze()
OISST_sst = Dataset(npth)['sst'][TIMES_co,:,:].squeeze()
ERROR = Dataset(npth)['err'][TIMES_co,:,:].squeeze()
anom = Dataset(npth)['anom'][TIMES_co,:,:].squeeze()

OISST_ice[OISST_ice<=-1000]=np.nan
mask=(OISST_ice!=OISST_ice).data
OISST_ice=ma.array(OISST_ice,mask=mask)

OISST_sst[OISST_sst<=-1000]=np.nan
mask=(OISST_sst!=OISST_sst).data
OISST_sst=ma.array(OISST_sst,mask=mask)

ERROR[ERROR<=-1000]=np.nan
mask=(ERROR!=ERROR).data
OISST_ice=ma.array(ERROR,mask=mask)

anom[anom<=-1000]=np.nan
mask=(anom!=anom).data
anom=ma.array(anom,mask=mask)


def mySST(nc_save_name,LON,LAT,TIME,Ref_time,values1,values2,values3,values4):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My OISST data '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('time', np.float64, ('time',))
    time.units=Ref_time
    time.field='time, scalar, series'
    
    DATA1 = ncfile.createVariable('sst',np.float64,('time','lat','lon'),compression='zlib') #
    DATA1.units = 'degree_C' 
    DATA1.long_name = 'OISST' 
    DATA1.coordinates = "time, lat, lon"
    
    DATA2 = ncfile.createVariable('ice',np.float64,('time','lat','lon'),compression='zlib') #
    DATA2.units = '%' 
    DATA2.long_name = 'OISST ice' 
    DATA2.coordinates = "time, lat, lon"
    
    DATA3 = ncfile.createVariable('anom',np.float64,('time','lat','lon'),compression='zlib') #
    DATA3.units = 'degree_C' 
    DATA3.long_name = 'OISST Daily sea surface temperature anomalies' 
    DATA3.coordinates = "time, lat, lon"
    
    DATA4 = ncfile.createVariable('err',np.float64,('time','lat','lon'),compression='zlib') #
    DATA4.units = 'degree_C' 
    DATA4.long_name = 'OISST Estimated error standard deviation of analysed_sst' 
    DATA4.coordinates = "time, lat, lon"
    
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME 
     
    DATA1[:] = values1
    DATA2[:] = values2
    DATA3[:] = values3
    DATA4[:] = values4

    ncfile.close()

mySST(wnpth,LON,LAT,my_time_num,My_time_ref,OISST_sst,OISST_ice,anom,ERROR)
















