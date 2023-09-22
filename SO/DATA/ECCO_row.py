from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import os
ECCO_ts_pth='/data4/tmp/TS/'
ECCO_vel_pth='/data4/tmp/VEL/'

wnpth='./myECCO_199201_201712_tsuvw.nc'

t_rng=['1980-01','2020-12']
My_time_ref='days since 1970-1-1 00:00:00'

### Read ECCO data ===============================
ts_list=list(np.sort([ECCO_ts_pth+i for i in os.listdir(ECCO_ts_pth) if i.endswith('nc')]))
vl_list=list(np.sort([ECCO_vel_pth+i for i in os.listdir(ECCO_vel_pth) if i.endswith('nc')]))
DATA=xr.open_mfdataset(ts_list,decode_times=False,parallel=True)
DATA2=xr.open_mfdataset(vl_list,decode_times=False)

TIME=DATA.time.values
TIME_units=DATA.time.attrs['units']
LON,LAT=DATA.longitude.values,DATA.latitude.values
DEPTH=DATA.Z.values
tmp_dif = date2num(dt.datetime(2000,1,1),My_time_ref) - date2num(dt.datetime(2000,1,1),'days since 1992-01-01T12:00:00')
TIME_re=TIME/24+tmp_dif
ECCO_THETA=DATA.THETA.values
ECCO_SALT=DATA.SALT.values

ECCO_EVEL=DATA2.EVEL.values
ECCO_NVEL=DATA2.NVEL.values
ECCO_WVEL=DATA2.WVEL.values

LAT_BNDS=DATA.latitude_bnds.values
LON_BNDS=DATA.longitude_bnds.values

posi_co,nega_co=np.where(LON>=0)[0],np.where(LON<0)[0]
LON_re=np.concatenate( [LON[posi_co],LON[nega_co]+360],axis=0 )

ECCO_THETA_re=np.concatenate( [ECCO_THETA[:,:,:,posi_co], ECCO_THETA[:,:,:,nega_co]] ,axis=3)
ECCO_SALT_re=np.concatenate( [ECCO_SALT[:,:,:,posi_co], ECCO_SALT[:,:,:,nega_co]] ,axis=3)

ECCO_EVEL_re=np.concatenate( [ECCO_EVEL[:,:,:,posi_co], ECCO_EVEL[:,:,:,nega_co]] ,axis=3)
ECCO_NVEL_re=np.concatenate( [ECCO_NVEL[:,:,:,posi_co], ECCO_NVEL[:,:,:,nega_co]] ,axis=3)
ECCO_WVEL_re=np.concatenate( [ECCO_WVEL[:,:,:,posi_co], ECCO_WVEL[:,:,:,nega_co]] ,axis=3)


def myOGCM(nc_save_name,LON,LAT,LON_bnds,LAT_bnds,DEPTH,TIME,Ref_time,values1,values2,values3,values4,values5):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('nv', 2)
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My ECCO data '
    
    lat_bnds = ncfile.createVariable('lat_bnds', np.float32, ('lat','nv'))
    lon_bnds = ncfile.createVariable('lon_bnds', np.float32, ('lon','nv'))
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    depth = ncfile.createVariable('depth', np.float32, ('depth',))
    depth.units = 'depth_m'
    time = ncfile.createVariable('time', np.float64, ('time',))
    time.units=Ref_time
    time.field='time, scalar, series'
    # time.cycle_length=cycle
    
    DATA1 = ncfile.createVariable('temp',np.float64,('time','depth','lat','lon'),compression='zlib') #
    DATA1.units = 'degree_C' 
    DATA1.long_name = 'ECCO temp' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    DATA2 = ncfile.createVariable('salt',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA2.units = '1e-3' 
    DATA2.long_name = 'ECCO salt' 
    DATA2.coordinates = "time, depth, lat, lon" 

    DATA3 = ncfile.createVariable('u',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA3.units = 'meter second-1' 
    DATA3.long_name = 'Eastward velocity' 
    DATA3.coordinates = "time, depth, lat, lon"
    
    DATA4 = ncfile.createVariable('v',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA4.units = 'meter second-1' 
    DATA4.long_name = 'Northward velocity' 
    DATA4.coordinates = "time, depth, lat, lon"
    
    DATA5 = ncfile.createVariable('w',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA5.units = 'meter second-1' 
    DATA5.long_name = 'Vertical velocity' #
    DATA5.coordinates = "time, depth, lat, lon"

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    time[:] = TIME 
    
    lat_bnds[:] = LAT_bnds
    lon_bnds[:] = LON_bnds
     
    DATA1[:] = values1
    DATA2[:] = values2
    DATA3[:] = values3
    DATA4[:] = values4
    DATA5[:] = values5

    ncfile.close()

myOGCM(wnpth,LON_re,LAT,LON_BNDS,LAT_BNDS,DEPTH,TIME,TIME_units,\
    ECCO_THETA_re,ECCO_SALT_re,\
        ECCO_EVEL_re,ECCO_NVEL_re,ECCO_WVEL_re)




