from netCDF4 import MFDataset,Dataset,num2date,date2num
import numpy as np
import datetime as dt
import xarray as xr
import os
ECCO_npth='/data4/tmp/myECCO_199201_201712_ts.nc'
ECCO_vel_pth='/data4/tmp/VEL/'

wnpth='./myECCO_199201_201712_tsuvw.nc'

t_rng=['1980-01','2022-12']
My_time_ref='days since 1970-1-1 00:00:00'

### Read ECCO data ===============================
vl_list=list(np.sort([ECCO_vel_pth+i for i in os.listdir(ECCO_vel_pth) if i.endswith('nc')]))
DATA2=xr.open_mfdataset(vl_list,decode_times=False)

LON=DATA2.longitude.values

ECCO_EVEL=DATA2.EVEL.values
ECCO_NVEL=DATA2.NVEL.values
ECCO_WVEL=DATA2.WVEL.values

posi_co,nega_co=np.where(LON>=0)[0],np.where(LON<0)[0]
LON_re=np.concatenate( [LON[posi_co],LON[nega_co]+360],axis=0 )

ECCO_EVEL_re=np.concatenate( [ECCO_EVEL[:,:,:,posi_co], ECCO_EVEL[:,:,:,nega_co]] ,axis=3)
ECCO_NVEL_re=np.concatenate( [ECCO_NVEL[:,:,:,posi_co], ECCO_NVEL[:,:,:,nega_co]] ,axis=3)
ECCO_WVEL_re=np.concatenate( [ECCO_WVEL[:,:,:,posi_co], ECCO_WVEL[:,:,:,nega_co]] ,axis=3)

nc=Dataset(ECCO_npth,'a')

DATA3 = nc.createVariable('u',np.float64,('time','depth','lat','lon'),compression='zlib') 
DATA3.units = 'meter second-1' 
DATA3.long_name = 'Eastward velocity' 
DATA3.coordinates = "time, depth, lat, lon"

DATA4 = nc.createVariable('v',np.float64,('time','depth','lat','lon'),compression='zlib') 
DATA4.units = 'meter second-1' 
DATA4.long_name = 'Northward velocity' 
DATA4.coordinates = "time, depth, lat, lon"

DATA5 = nc.createVariable('w',np.float64,('time','depth','lat','lon'),compression='zlib') 
DATA5.units = 'meter second-1' 
DATA5.long_name = 'Vertical velocity' #
DATA5.coordinates = "time, depth, lat, lon"

DATA3[:] = ECCO_EVEL_re
DATA4[:] = ECCO_NVEL_re
DATA5[:] = ECCO_WVEL_re

nc.close()





