# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:20:57 2022

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

import numpy as np
import datetime as dt
import os
from netCDF4 import Dataset, num2date, date2num

r_path='D:/Working_hub/OneDrive/base142/Warehouse01/ERA_Frc/'

file_list = os.listdir(r_path)
ERA_list=[i for i in file_list if i.endswith('.nc')]

for i in ERA_list:
    nc=Dataset(r_path+i)
    # nc.variables.keys()

    var_name=list(nc.variables.keys())[-1]
    time_name=list(nc.variables.keys())[-2]
    nc_save_name=r_path+'CLM'+'/ERA5_'+var_name+'_CLM.nc'
        
    Make_clm(nc,var_name,nc_save_name,time_name)
    nc.close()



# A=Dataset('D:/Working_hub/OneDrive/Sources/Data/ERA_GLO_mon/adaptor.mars.internal.nc')


# nc_save_name='ERA5_dQdSST_CLM.nc'
def Make_clm(nc,var_name,nc_save_name,time_name):
    print('!!! '+var_name+' !!!')
    LAT=nc['lat']
    LON=nc['lon']
    Time_=nc[time_name]
    Time_unit='days since 1950-01-01'
    Time=num2date(Time_[:],Time_unit)
    
    Time_array=np.array([[i.year, i.month, i.day] for i in Time])
    Var=nc[var_name]
    T,L,M=Var.shape
    
    empty_cell=[]
    
    for i in range(1,13):
        month_co=np.where(Time_array[:,1]==i)[0]
        empty_cell.append(np.mean(Var[month_co],axis=0))
    Full_cell=np.array(empty_cell)
    
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF3_CLASSIC')


    ncfile.createDimension('lat', len(nc['lat'][:]))
    ncfile.createDimension('lon', len(nc['lon'][:]))
    ncfile.createDimension(time_name, 12)

    ncfile.title='ERA5 Climatology (1979~2021)'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',))
    lat.units = LAT.units
    # lat.long_name = LAT.long_name
    lon = ncfile.createVariable('lon', np.float32, ('lon',))
    lon.units = LON.units
    # lon.long_name = LON.long_name
    time = ncfile.createVariable(time_name, np.float64, (time_name,))
    time.long_name = time_name
    time.units='time since initialization'
    time.cycle_length=365.25
    
    Data = ncfile.createVariable(var_name,np.float64,(time_name,'lat','lon')) # note: unlimited dimension is leftmost
    Data.units = Var.units # degrees Kelvin
    Data.standard_name = Var.long_name # this is a CF standard name
    # Data.field=Var.field
    Data.time=Var.time
    Data.coordinates= Var.coordinates

    lat[:] = LAT[:]
    lon[:] = LON[:]
    time[:] = np.arange(15.2188,351.0313,30.4375)  
    Data[:] = Full_cell
    
    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()

# Make_clm(nc,var_name,nc_save_name,time_name)



