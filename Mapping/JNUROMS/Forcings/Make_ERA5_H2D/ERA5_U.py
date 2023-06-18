#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:48:37 2023

@author: birostris
@Email: birostris36@gmail.com 

@Category:  
@Reference: None
@Descriptions:
"""

import sys
sys.path.append('/data2/base132/tmp_MantaROMS/Make_ERA5_H2D/')
import matplotlib.pyplot as plt
from netCDF4 import Dataset, date2num, num2date
import pandas as pd
import numpy as np
from tqdm import tqdm
import datetime as dt
import time

read_path='/data2/ERA5/'
vname='u10'
save_path='/data2/ERA5_model/'
save_name='ERA5_daily_uwind_1979_2022.nc'


nc=Dataset(read_path+'ERA5_uwind_S_1979_2022.nc')

Sample_data=nc[vname][0]
at,on=Sample_data.shape
LON,LAT=nc['longitude'][:],nc['latitude'][:]

TIMES=nc['time']

# =============================================================================
t_rng=['1979-01','2022-12']
My_time_ref='days since 1970-1-1 00:00:00'
Times_units=TIMES.units.replace('hours','days')
TIMES_values=TIMES[:]

TIMES_s2d= np.sort(np.array(list(set([int(i) for i in TIMES_values/24]))))
TIMES_real01=num2date(TIMES_s2d,Times_units)
TIMES_real02=num2date(TIMES_values,TIMES.units)

Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
TIMES_co=np.where( (TIMES_real01>=Tst)&(TIMES_real01<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),Times_units)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
My_time=num2date(TIMES_s2d[TIMES_co]-tmp_dif,My_time_ref)
My_time_num=TIMES_s2d[TIMES_co]-tmp_dif
# =============================================================================

TIMES_str1=[ i.strftime('%Y-%m-%d').split('-') for i in TIMES_real02 ]
TIMES_str2=[ i.strftime('%Y-%m-%d') for i in TIMES_real02 ]

Times_set=np.sort(list(set(TIMES_str2)))
Times_guide=pd.DataFrame([i.split('-') for i in Times_set],columns=['Year','Month','Day'])
Times=pd.DataFrame(TIMES_str1,columns=['Year','Month','Day'])

# =============================================================================
if len(Times)%4!=0:
    raise ValueError

LAT_new=np.flipud(LAT)
# Make Pacific centered
posi_co,nega_co=np.where(LON>=0)[0],np.where(LON<0)[0]
LON_new=np.concatenate([LON[posi_co],360+LON[nega_co]])

stT=time.time()
N=0; MY_Var=np.zeros([int(len(TIMES_str1)/4),at,on])
for YY in [f'{i:02d}' for i in range(1979,2023)]:
    print('!!! Year : ',YY)
    for MM in [f'{i:02d}' for i in range(1,13)]:
        for DD in Times_guide['Day'].loc[(Times_guide['Year']==YY)&(Times_guide['Month']==MM)]:
            tmp_co=np.where( (Times['Year']==YY) & (Times['Month']==MM) &\
                            (Times['Day']==DD) )[0]
            tmp_var=np.flip(np.mean(nc[vname][tmp_co,:,:],axis=0),axis=0)
            posi_var,nega_var = tmp_var[:,posi_co],tmp_var[:,nega_co]
            MY_Var[N]=np.concatenate([posi_var,nega_var],axis=1) 
            N+=1
print('\n!!! END Total Elapsed time : '+str((time.time()-stT)/60)[:4]+'min !!!')

# MY_Var[MY_Var<-1000]=0

stT=time.time()
create_uwind_nc(save_path+save_name,LON_new,LAT_new,My_time_num,My_time_ref,MY_Var)
print('\n!!! END Total Elapsed time : '+str((time.time()-stT)/60)[:4]+'min !!!')

def create_uwind_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('wind_time',len(TIME))
    
    ncfile.title='ERA5 daily uwind'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',),zlib=False)
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',),zlib=False)
    lon.units = 'degrees_east'
    time = ncfile.createVariable('wind_time', np.float64, ('wind_time',),zlib=False)
    time.units=Ref_time
    time.field='u-wind_time, scalar, series'

    Uwind = ncfile.createVariable('Uwind',np.float64,('wind_time','lat','lon'),\
                                zlib=False, fill_value=1e+20) #
    Uwind.units =  'meter secind-1'
    Uwind.long_name = 'surface u-wind component' # this is a CF standard name
    Uwind.time='wind_time'
    Uwind.coordinates = "lon lat"
    Uwind.field='u-wind, scalar, series'

    # Data.field=Var.field
    lat[:] = LAT
    lon[:] = LON
    time[:] = TIME  
    Uwind[:] = values
    # print("-- Wrote data, temp.shape is now ", temp.shape)
    # print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())
    ncfile.close()











































