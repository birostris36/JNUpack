import xarray as xr
import os
import numpy as np
from netCDF4 import Dataset,num2date
from scipy.interpolate import NearestNDInterpolator
from copy import deepcopy
import sys

wnpth='/home/shjo/_data/myARMOR3D_SO_199301_202312_ts.nc'
pth='/home/shjo/_data/AMOGO/'
ARMOR3D_lst=list(np.sort([pth+i for i in os.listdir(pth) if i.endswith('.nc')]))

Sample=Dataset(ARMOR3D_lst[0])
LAT,LON,DEPTH=Sample['latitude'][:],Sample['longitude'][:],Sample['depth'][:]
LAT_co=np.where( (LAT>=-80)&(LAT<=-10) )[0]
LAT=LAT[LAT_co]
DEPTH_co=np.where( (DEPTH>=0)&(DEPTH<=2000) )[0]
DEPTH=DEPTH[DEPTH_co]
TIME_ref=Sample['time'].units.replace('hours','days')

T,D,A,O=len(ARMOR3D_lst),len(DEPTH),len(LAT),len(LON)

myTEMP=np.zeros([T,D,A,O])
mySALT=np.zeros_like(myTEMP)
myTIME=np.zeros([T])
### Read DATA ==============================================
print('!!! Read Data !!!')

for i,nn in zip(ARMOR3D_lst,range(T)):
    tmp=xr.open_dataset(i,decode_times=False).loc[dict(depth=slice(0,2000),latitude=slice(-80,-10))]
    myTEMP[nn]=tmp['to'].values
    mySALT[nn]=tmp['so'].values
    myTIME[nn]=tmp['time'].values
    if nn%30==0:
        print( '!!! '+str(num2date(myTIME[nn],TIME_ref))+' !!!'   )
        print('!!! '+i+' !!!')
        
### Writing data =============================================
print('!!! Writing data... !!!')
def myOGCM(nc_save_name,LON,LAT,DEPTH,TIME,Ref_time,values1,values2):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My Multi Observation Global Ocean 3D data '
    
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
    DATA1.long_name = 'Multi Observation Global Ocean temp' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    DATA2 = ncfile.createVariable('salt',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA2.units = 'g kg-1' 
    DATA2.long_name = 'Multi Observation Global Ocean (practical) salt' 
    DATA2.coordinates = "time, depth, lat, lon"
    
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    time[:] = TIME 
     
    DATA1[:] = values1
    DATA2[:] = values2

    ncfile.close()

myOGCM(wnpth,LON,LAT,DEPTH,myTIME,TIME_ref,myTEMP,mySALT)