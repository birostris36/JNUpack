import xarray as xr
import os
import numpy as np
from netCDF4 import Dataset,num2date
from scipy.interpolate import NearestNDInterpolator
from copy import deepcopy
import sys
# sys.path.append('C:/Users/shjo/Bridge/JNUpack/mySO_src/libs/')
# from myTools import NearesND
    
wnpth='J:/myCORA_1rg_196001_202212_ts.nc'
pth='D:/CORA/Grid/'
cora_temp_lst=list(np.sort([pth+i for i in os.listdir(pth) if i.endswith('TEMP.nc')]))
cora_salt_lst=list(np.sort([pth+i for i in os.listdir(pth) if i.endswith('PSAL.nc')]))

Sample=xr.open_dataset(cora_temp_lst[0], decode_times=False)
LAT,LON=Sample.latitude.values,Sample.longitude.values
DEPTH=Sample.depth.values
posi_co,nega_co=np.where(LON>=0)[0],np.where(LON<0)[0]
LON_re=np.concatenate( [LON[posi_co],LON[nega_co]+360],axis=0 )

TIME_units=Sample.time.units

myLAT=np.arange(-77,80,1 )
myLON=np.arange(0.5,360,1)

def NearesND(lon1, lat1, data, lon2, lat2):

    if len(np.shape(lon1)) == 1 or len(np.shape(lat1)) == 1:
        lon1, lat1 = np.meshgrid(lon1, lat1)
    
    if len(np.shape(lon2)) == 1 or len(np.shape(lat2)) == 1:
        lon2, lat2 = np.meshgrid(lon2, lat2)

    dom_lon_flat = lon1.ravel()
    dom_lat_flat = lat1.ravel()
    linear_filled_flat = data.ravel()
    valid_indices = ~np.isnan(linear_filled_flat)
    NearInterp = NearestNDInterpolator((dom_lon_flat[valid_indices], 
                                     dom_lat_flat[valid_indices]),
                                     linear_filled_flat[valid_indices])
    neares_filled = NearInterp(lon2, lat2)
    return neares_filled

print('!!! Read Data !!!')
temp_re=np.zeros([len(cora_temp_lst),len(DEPTH),len(myLAT),len(myLON)])
salt_re=np.zeros_like(temp_re)
myTIME=np.zeros(len(cora_temp_lst))
for temp_,salt_,nn in zip(cora_temp_lst,cora_salt_lst,range(len(cora_temp_lst))):

    temp,salt=Dataset(temp_)['TEMP'][:],Dataset(salt_)['PSAL'][:]
    temp=np.concatenate( [temp[:,:,:,posi_co], temp[:,:,:,nega_co]] ,axis=3).squeeze()
    salt=np.concatenate( [salt[:,:,:,posi_co], salt[:,:,:,nega_co]] ,axis=3).squeeze()
    
    myTIME[nn]=Dataset(temp_)['time'][:]
    
    if nn%30==0:
        print('!!! '+str(num2date(myTIME[nn],TIME_units)+' !!!' ))
    
    for nd in range(len(DEPTH)):
        temp_re[nn,nd] = NearesND(LON_re, LAT, temp[nd], myLON, myLAT)
        salt_re[nn,nd] = NearesND(LON_re, LAT, salt[nd], myLON, myLAT)

### Writing data =============================================
print('!!! Writing data... !!!')
def myOGCM(nc_save_name,LON,LAT,DEPTH,TIME,Ref_time,values1,values2):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My CORE regrided (x1) data '
    
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
    DATA1.long_name = 'CORA temp' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    DATA2 = ncfile.createVariable('salt',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA2.units = 'g kg-1' 
    DATA2.long_name = 'CORA salt' 
    DATA2.coordinates = "time, depth, lat, lon"
    
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    time[:] = TIME 
     
    DATA1[:] = values1
    DATA2[:] = values2

    ncfile.close()

myOGCM(wnpth,myLON,myLAT,DEPTH,myTIME,TIME_units,temp_re,salt_re)


