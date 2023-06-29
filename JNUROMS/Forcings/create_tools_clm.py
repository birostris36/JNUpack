# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:56:58 2023

@author: shjo9
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:25:43 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

from netCDF4 import Dataset
import xarray as xr
import numpy as np

def get_dqdsst(sst,sat,rho_atm,U,qsea):
    '''
    %  sst     : sea surface temperature (Celsius)
    %  sat     : sea surface atmospheric temperature (Celsius)
    %  rho_atm : atmospheric density (kilogram meter-3) 
    %  U       : wind speed (meter s-1)
    %  qsea    : sea level specific humidity
    '''
    # Specific heat of atmosphere.
    Cp=1004.8
    # Sensible heat transfert coefficient (stable condition)
    Ch = 0.66e-3
    # Latent heat transfert coefficient (stable condition)
    Ce = 1.15e-3
    # Emissivity coefficient
    eps = 0.98
    # Stefan constant
    stef = 5.6697e-8;
    # SST (KELVIN)
    SST = sst + 273.15;
    # Latent heat of vaporisation (J.kg-1)
    L = 2.5008e6 - 2.3e3 * sat
    # Infrared contribution
    q1 = -4. * stef * (SST**3)
    # Sensible heat contribution
    q2 = -rho_atm * Cp * Ch * U
    # Latent heat contribution
    dqsdt = 2353.* np.log(10.) * qsea / (SST**2)
    q3 = -rho_atm * Ce * L * U * dqsdt
    dqdsst = q1 + q2 + q3 
    return dqdsst

def create_wind_nc(nc_save_name,LON,LAT,TIME,Ref_time,values1,values2):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('wind_time',len(TIME))
    
    ncfile.title='ERA5 daily winds'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('wind_time', np.float64, ('wind_time',))
    time.units=Ref_time
    time.field='qair_time, scalar, series'
    time.cycle_length=360
    
    Uwind = ncfile.createVariable('Uwind',np.float64,('wind_time','lat','lon')) # note: unlimited dimension is leftmost
    Uwind.units = 'meter second-1' 
    Uwind.long_name = 'surface u-wind component' # this is a CF standard name
    Uwind.time='wind_time'
    Uwind.coordinates = "lon lat"
    Uwind.field='u-wind, scalar, series'
    
    Vwind = ncfile.createVariable('Vwind',np.float64,('wind_time','lat','lon')) # note: unlimited dimension is leftmost
    Vwind.units = 'meter second-1' 
    Vwind.long_name = 'surface v-wind component' # this is a CF standard name
    Vwind.time='wind_time'
    Vwind.coordinates = "lon lat"
    Vwind.field='v-wind, scalar, series'


    # Data.field=Var.field

    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Uwind[:] = values1
    Vwind[:] = values2

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
    
def create_tair_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('tair_time',len(TIME))
    
    ncfile.title='ERA5 daily tair'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('tair_time', np.float64, ('tair_time',))
    time.units=Ref_time
    time.field='tair_time, scalar, series'
    time.cycle_length=360

    Tair = ncfile.createVariable('Tair',np.float64,('tair_time','lat','lon')) # note: unlimited dimension is leftmost
    Tair.units = 'Celsius' 
    Tair.long_name = 'surface air temperature at 2m, for CORE' # this is a CF standard name
    Tair.time='tair_time'
    Tair.coordinates = "lon lat"
    Tair.field='Tair, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Tair[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
        
def create_pair_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('pair_time',len(TIME))
    
    ncfile.title='ERA5 daily pair (1993)'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('pair_time', np.float64, ('pair_time',))
    time.units=Ref_time
    time.field='pair_time, scalar, series'
    time.cycle_length=360

    
    Pair = ncfile.createVariable('Pair',np.float64,('pair_time','lat','lon')) # note: unlimited dimension is leftmost
    Pair.units = 'millibar' 
    Pair.long_name = 'surface air pressure, for CORE' # this is a CF standard name
    Pair.time='pair_time'
    Pair.coordinates = "lon lat"
    Pair.field='Tair, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Pair[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    

def create_qair_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('qair_time',len(TIME))
    
    ncfile.title='ERA5 daily qair'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('qair_time', np.float64, ('qair_time',))
    time.units=Ref_time
    # time.cycle_length=cycle
    time.field='qair_time, scalar, series'
    time.cycle_length=360

    
    Qair = ncfile.createVariable('Qair',np.float64,('qair_time','lat','lon')) # note: unlimited dimension is leftmost
    Qair.units = 'percentage' 
    Qair.long_name = 'surface air specific humidity, for CORE' # this is a CF standard name
    Qair.time='qair_time'
    Qair.coordinates = "lon lat"
    Qair.field='Qair, scalar, series'
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Qair[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
    
def create_cloud_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('cloud_time',len(TIME))
    
    ncfile.title='ERA5 daily cloud'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('cloud_time', np.float64, ('cloud_time',))
    time.units=Ref_time
    # time.cycle_length=cycle
    time.field='cloud_time, scalar, series'
    time.cycle_length=360

    
    Cloud = ncfile.createVariable('Cloud',np.float64,('cloud_time','lat','lon')) # note: unlimited dimension is leftmost
    Cloud.units = 'nondimensional' 
    Cloud.long_name = '-' # this is a CF standard name
    Cloud.time='cloud_time'
    Cloud.coordinates = "lon lat"
    Cloud.field='Cloud, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Cloud[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    

def create_sst_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('sst_time',len(TIME))
    
    ncfile.title='ERA5 daily sst'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('sst_time', np.float64, ('sst_time',))
    time.units=Ref_time
    time.field='sst_time, scalar, series'
    time.cycle_length=360

    
    SST = ncfile.createVariable('sst',np.float64,('sst_time','lat','lon')) # note: unlimited dimension is leftmost
    SST.units =  'Celsius'
    SST.long_name = 'sea surface temperature' # this is a CF standard name
    SST.time='sst_time'
    SST.coordinates = "lon lat"
    SST.field='SST, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    SST[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
def create_dqdsst_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('sst_time',len(TIME))
    
    ncfile.title='ERA5 daily dqdsst '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('sst_time', np.float64, ('sst_time',))
    time.units=Ref_time
    time.field='sst_time, scalar, series'
    time.cycle_length=360

    
    dqdsst = ncfile.createVariable('dqdsst',np.float64,('sst_time','lat','lon')) # note: unlimited dimension is leftmost
    dqdsst.units = 'Watts meter-2 Celsius-1' 
    dqdsst.long_name = 'surface net heat flux sensitivity to SST' # this is a CF standard name
    dqdsst.time='sst_time'
    dqdsst.coordinates = "lon lat"
    dqdsst.field='dQdSST, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    dqdsst[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
def create_srf_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('srf_time',len(TIME))
    
    ncfile.title='ERA5 daily srf '
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('srf_time', np.float64, ('srf_time',))
    time.units=Ref_time
    time.field='srf_time, scalar, series'
    time.cycle_length=360

    
    srf = ncfile.createVariable('srf',np.float64,('srf_time','lat','lon')) # note: unlimited dimension is leftmost
    srf.units = 'Watt meter-2'
    srf.long_name = 'shortwave radiation, scalar, series' # this is a CF standard name
    srf.time='srf_time'
    srf.positive='downward flux, heating'
    srf.negative='upward flux, cooling'
    srf.coordinates = "lon lat"
    srf.field='shortwave radiation, scalar, series'


    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    srf[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()
    
    
def create_lwrad_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('lrf_time',len(TIME))
    
    ncfile.title='ERA5 daily lwrad'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('lrf_time', np.float64, ('lrf_time',))
    time.units=Ref_time
    time.field='lrf_time, scalar, series'
    time.cycle_length=360

    
    lrf = ncfile.createVariable('lwrad',np.float64,('lrf_time','lat','lon')) # note: unlimited dimension is leftmost
    lrf.units = 'Watt meter-2'
    lrf.long_name = 'Net longwave radiation' # this is a CF standard name
    lrf.time='lrf_time'
    lrf.positive='downward flux, heating'
    lrf.negative='upward flux, cooling'
    lrf.coordinates = "lon lat"
    lrf.field = 'longwave radiation, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    lrf[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()   
    
    
def create_lwrad_down_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('lrf_time',len(TIME))
    
    ncfile.title='ERA5 daily lwrad_down'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('lrf_time', np.float64, ('lrf_time',))
    time.units=Ref_time
    time.field='lrf_time, scalar, series'
    time.cycle_length=360

    lwrad_down = ncfile.createVariable('lwrad_down',np.float64,('lrf_time','lat','lon')) # note: unlimited dimension is leftmost
    lwrad_down.units = 'Watt meter-2'
    lwrad_down.long_name = 'Downward longwave radiation' # this is a CF standard name
    lwrad_down.time='lrf_time'
    lwrad_down.positive='downward flux, heating'
    lwrad_down.negative='upward flux, cooling'
    lwrad_down.coordinates = "lon lat"
    lwrad_down.field = 'longwave radiation, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    lwrad_down[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()   
    
    
    
   
def create_rain_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('rain_time',len(TIME))
    
    ncfile.title='ERA5 daily rain'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = 'degrees_east'
    time = ncfile.createVariable('rain_time', np.float64, ('rain_time',))
    time.units=Ref_time
    time.field='rain_time, scalar, series'
    time.cycle_length=360

    
    rain = ncfile.createVariable('rain',np.float64,('rain_time','lat','lon')) # note: unlimited dimension is leftmost
    rain.units = 'kg mm-2 s-1'
    rain.long_name = 'rain fall rate' # this is a CF standard name
    rain.time='rain_time'
    rain.coordinates = "lon lat"
    rain.field = 'rain, scalar, series'

    
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    rain[:] = values

    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()   
    
    
    
    
    
    