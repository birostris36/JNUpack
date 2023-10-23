import os
import numpy as np
from netCDF4 import Dataset,MFDataset
import time
import xarray as xr

wnpth='/home/shjo/_data/ARGO/myARGO_200201_201712_ts.nc'
myARGO_pth='/home/shjo/_data/ARGO/field/'

f_lst=[myARGO_pth+i for i in os.listdir(myARGO_pth)]

myARGO=[]
for i in f_lst:
    [myARGO.append(i+'/'+j) for j in os.listdir(i)]
    
myARGO_temp=[i for i in myARGO if i.endswith('TEMP.nc')]
myARGO_salt=[i for i in myARGO if i.endswith('PSAL.nc')]

print('!!! Argo TEMP list : !!!')
print(myARGO_temp[0])
print(myARGO_temp[-1])
print('DATA len : ',len(myARGO_temp))
time.sleep(10)
print('!!! Argo TEMP list : !!!')
print(myARGO_salt[0])
print(myARGO_salt[-1])
print('DATA len : ',len(myARGO_salt))
time.sleep(10)

Sample=Dataset(myARGO_temp[0])

LON,LAT,DEPTH,TIME=Sample['longitude'][:],Sample['latitude'][:],Sample['depth'][:],Sample['time']

TIME_units=TIME.units
TIME_val  =TIME[:]

print('!!! Read Data !!!')
TEMP=xr.open_mfdataset(myARGO_temp).TEMP.loc[dict(depth=slice(0,2000))].values

print('!!! Reshape !!!')
posi_co,nega_co=np.where(LON>=0)[0],np.where(LON<0)[0]
LON_re=np.concatenate( [LON[posi_co],LON[nega_co]+360],axis=0 )

TEMP=np.concatenate( [TEMP[:,:,:,posi_co], TEMP[:,:,:,nega_co]] ,axis=3)

def myOGCM(nc_save_name,LON,LAT,DEPTH,TIME,Ref_time,values1):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My ARGO data '
    
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
    DATA1.long_name = 'ARGO temp' 
    DATA1.coordinates = "time, depth, lat, lon"

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    time[:] = TIME 
     
    DATA1[:] = values1
    # DATA2[:] = values2

    ncfile.close()
    
print('!!! Write Data !!!')

myOGCM(wnpth,LON_re,LAT,DEPTH,TIME_val,TIME_units,\
    TEMP)




