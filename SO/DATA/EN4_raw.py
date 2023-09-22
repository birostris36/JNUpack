from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import os

EN4_pth='E:/_data/EN4_422_c14/'
wnpth='J:/Reanalysis/myEN4_198001_202012_ts.nc'

t_rng=['1980-01','2020-12']
My_time_ref='days since 1970-1-1 00:00:00'
d_list=list(np.sort([EN4_pth+i for i in os.listdir(EN4_pth) if i.endswith('.nc')]))

Sample=Dataset(d_list[0])
LON,LAT=Sample['lon'][:],Sample['lat'][:]
DEPTH,DEPTH_BNDS=Sample['depth'][:],Sample['depth_bnds'][:]

TIME=MFDataset(d_list)['time']
TIMES=TIME[:]
UNITS=TIME.units

# Process Times
tmp_time_var='time'
t_rng=['1980-01','2022-12']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=TIMES
TIME_UNIT=UNITS
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
my_time=num2date(OGCM_TIMES[TIMES_co]-tmp_dif,My_time_ref)
my_time_num=OGCM_TIMES[TIMES_co]-tmp_dif

### Read data ================================
print('!!! Read data... !!!')
temp=MFDataset(d_list)['temperature'][TIMES_co]-273.15
salt=MFDataset(d_list)['salinity'][TIMES_co]

### Write data ================================
print('!!! Write data... !!!')

def myOGCM(nc_save_name,LON,LAT,DEPTH,DEPTH_BNDS,TIME,Ref_time,values1,values2):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('nv', 2)
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My EN4 data '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    depth = ncfile.createVariable('depth', np.float32, ('depth',))
    depth.units = 'depth_m'
    
    depth_bnds = ncfile.createVariable('depth_bnds', np.float32, ('depth','nv'))
    depth_bnds.units = 'depth_m'
    
    time = ncfile.createVariable('time', np.float64, ('time',))
    time.units=Ref_time
    time.field='time, scalar, series'
    # time.cycle_length=cycle
    
    DATA1 = ncfile.createVariable('temp',np.float64,('time','depth','lat','lon'),compression='zlib') #
    DATA1.units = 'degree_C' 
    DATA1.long_name = 'EN4 potential_temperature' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    DATA2 = ncfile.createVariable('salt',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA2.units = 'g kg-1' 
    DATA2.long_name = 'EN4 salt' 
    DATA2.coordinates = "time, depth, lat, lon"

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    depth_bnds[:]=DEPTH_BNDS
    time[:] = TIME 
     
    DATA1[:] = values1
    DATA2[:] = values2

    ncfile.close()
    
myOGCM(wnpth,LON,LAT,DEPTH,DEPTH_BNDS,my_time_num,My_time_ref,temp,salt)























