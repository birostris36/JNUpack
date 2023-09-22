from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import os

GECCO_pth='E:/_data/GECCO/'
wnpth='J:/Reanalysis/myGECCO_199201_201712_ztsuv.nc'

Sample=Dataset(GECCO_pth+'GECCO_SALT_194801_201812_ori.nc')

LAT,LON=Sample['lat'][:],Sample['lon'][:]
DEPTH=Sample['Depth'][:]

posi_co,nega_co=np.where(LON>=0)[0],np.where(LON<0)[0]
LON_re=np.concatenate( [LON[posi_co],LON[nega_co]+360],axis=0 )

# Process Times
tmp_time_var='time'
t_rng=['1980-01','2020-12']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=Sample[tmp_time_var]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
my_time=num2date(OGCM_TIMES[TIMES_co]-tmp_dif,My_time_ref)
my_time_num=OGCM_TIMES[TIMES_co]-tmp_dif

GECCO_ZETA = Dataset(GECCO_pth+'GECCO_ZETA_194801_201812_ori.nc')['zeta'][TIMES_co,:,:]

GECCO_SALT = Dataset(GECCO_pth+'GECCO_SALT_194801_201812_ori.nc')['salt'][TIMES_co,:,:,:]
GECCO_TEMP = Dataset(GECCO_pth+'GECCO_TEMP_194801_201812_ori.nc')['temp'][TIMES_co,:,:,:]

GECCO_UVEL = Dataset(GECCO_pth+'GECCO_UVEL_194801_201812_ori.nc')['u'][TIMES_co,:,:,:]
GECCO_VVEL = Dataset(GECCO_pth+'GECCO_VVEL_194801_201812_ori.nc')['v'][TIMES_co,:,:,:]

GECCO_ZETA_re=np.concatenate( [GECCO_ZETA[:,:,posi_co], GECCO_ZETA[:,:,nega_co]] ,axis=2)

GECCO_TEMP_re=np.concatenate( [GECCO_TEMP[:,:,:,posi_co], GECCO_TEMP[:,:,:,nega_co]] ,axis=3)
GECCO_SALT_re=np.concatenate( [GECCO_SALT[:,:,:,posi_co], GECCO_SALT[:,:,:,nega_co]] ,axis=3)
GECCO_UVEL_re=np.concatenate( [GECCO_UVEL[:,:,:,posi_co], GECCO_UVEL[:,:,:,nega_co]] ,axis=3)
GECCO_VVEL_re=np.concatenate( [GECCO_VVEL[:,:,:,posi_co], GECCO_VVEL[:,:,:,nega_co]] ,axis=3)


def myOGCM(nc_save_name,LON,LAT,DEPTH,TIME,Ref_time,values1,values2,values3,values4,values5):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My GECCO data '
    
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
    DATA1.long_name = 'GECCO temp' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    DATA2 = ncfile.createVariable('salt',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA2.units = 'g kg-1' 
    DATA2.long_name = 'GECCO salt' 
    DATA2.coordinates = "time, depth, lat, lon"

    DATA3 = ncfile.createVariable('u',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA3.units = 'meter second-1' 
    DATA3.long_name = 'Eastward velocity' 
    DATA3.coordinates = "time, depth, lat, lon"
    
    DATA4 = ncfile.createVariable('v',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA4.units = 'meter second-1' 
    DATA4.long_name = 'Northward velocity' 
    DATA4.coordinates = "time, depth, lat, lon"
    
    DATA5 = ncfile.createVariable('zeta',np.float64,('time','lat','lon'),compression='zlib') 
    DATA5.units = 'meter' 
    DATA5.long_name = 'model ssh' #
    DATA5.coordinates = "time, lat, lon"

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    time[:] = TIME 
     
    DATA1[:] = values1
    DATA2[:] = values2
    DATA3[:] = values3
    DATA4[:] = values4
    DATA5[:] = values5

    ncfile.close()


myOGCM(wnpth,LON_re,LAT,DEPTH,my_time_num,My_time_ref,\
    GECCO_TEMP_re,GECCO_SALT_re,\
        GECCO_UVEL_re,GECCO_VVEL_re,GECCO_ZETA_re)