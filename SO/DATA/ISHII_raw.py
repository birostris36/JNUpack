from scipy import io
from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import pandas as pd
from numpy import ma
import os

pth='D:/JMA/CHA_MAT/'
wnpth='J:/Obs/myISHII_198001_201912_ts.nc'

d_lst=[pth+i for i in os.listdir(pth) if i.endswith('mat')]

Sample=io.loadmat(d_lst[0])
LAT,LON=Sample['lat'],Sample['lon']
DEPTH=np.flipud(Sample['depth'])
### Read mat files =================================================
TEMP=np.zeros([len(d_lst)*12,len(DEPTH),len(LAT),len(LON)])
SALT=np.zeros_like(TEMP)
for i,fac in zip(d_lst,range(len(d_lst))):
    tmp=io.loadmat(i)
    TEMP[fac*12:(fac+1)*12,:,:,:]=tmp['temp']
    SALT[fac*12:(fac+1)*12,:,:,:]=tmp['salt']

### Create time num =================================================
TIMES=pd.date_range('1980-01','2020-01',freq='M').strftime('%Y%m')
Ref_time='days since 1970-1-1 00:00:00'
my_time_num=date2num([dt.datetime(int(str(i)[:4]),int(str(i)[4:6]),15)\
    for i in TIMES],Ref_time)

### Data Masking ====================================================
mask=(TEMP!=TEMP).data
TEMP=ma.array(TEMP,mask=mask)
mask=(SALT!=SALT).data
SALT=ma.array(SALT,mask=mask)
TEMP,SALT=np.flip(TEMP,axis=1),np.flip(SALT,axis=1)


### Data Writing =====================================================
def myOGCM(nc_save_name,LON,LAT,DEPTH,TIME,Ref_time,values1,values2):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('depth',len(DEPTH))
    ncfile.createDimension('time',len(TIME))
    
    ncfile.title='My Ishii (JMA) data '
    
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
    DATA1.long_name = 'Ishii (JMA) temp' 
    DATA1.coordinates = "time, depth, lat, lon"
    
    DATA2 = ncfile.createVariable('salt',np.float64,('time','depth','lat','lon'),compression='zlib') 
    DATA2.units = 'g kg-1 (?)' 
    DATA2.long_name = 'Ishii (JMA) salt' 
    DATA2.coordinates = "time, depth, lat, lon"
    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    depth[:]= DEPTH
    time[:] = TIME 
     
    DATA1[:] = values1
    DATA2[:] = values2

    ncfile.close()
    
myOGCM(wnpth,LON,LAT,DEPTH,my_time_num,Ref_time,TEMP,SALT)
    









