#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 17:34:08 2023

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
import xarray as xr

read_path='/data2/ERA5/'

# =============================================================================
# A=Dataset('/data2/ERA5_model/ERA5_daily_dQdSST_1979_2022.nc')
# A=xr.open_dataset('/data2/ERA5_model/ERA5_daily_dQdSST_1979_2022.nc')

# MY_Var=np.load('/data2/MY_dQdSST.npy')

# =============================================================================

vname01='sst'
vname02='t2m'
vname03='d2m'

vname04='u10'
vname05='v10'

save_path='/data2/ERA5_model/'
save_name='ERA5_daily_dQdSST_1979_2022.nc'

SST=Dataset(read_path+'ERA5_SST_S_1979_2022.nc')
T2=Dataset(read_path+'ERA5_2temp_S_1979_2022.nc')
D2=Dataset(read_path+'ERA5_2dew_S_1979_2022.nc')
Uwind=Dataset(read_path+'ERA5_uwind_S_1979_2022.nc')
Vwind=Dataset(read_path+'ERA5_vwind_S_1979_2022.nc')
# Qair=Dataset(read_path02+'.nc')



Sample_data=SST[vname01][0]
at,on=Sample_data.shape
LON,LAT=SST['longitude'][:],SST['latitude'][:]

TIMES=SST['time']

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
                
            tmp_sst=np.flip(np.mean(SST[vname01][tmp_co,:,:],axis=0),axis=0)-273.15

            tmp_u=np.flip(np.mean(Uwind[vname04][tmp_co,:,:],axis=0),axis=0)
            tmp_v=np.flip(np.mean(Vwind[vname05][tmp_co,:,:],axis=0),axis=0)
                
            tmp_T2=np.flip(np.mean(T2[vname02][tmp_co,:,:],axis=0),axis=0)-273.15
            tmp_D2=np.flip(np.mean(D2[vname03][tmp_co,:,:],axis=0),axis=0)-273.15
            
            tmp_qair = 100 * ( np.exp( (17.625 * tmp_D2 ) /  (243.04 + tmp_D2 ) ) / 
                     np.exp( ( 17.625 * tmp_T2 ) / (243.04 + tmp_T2 ) ) )
            
            tmp_dqdsst = get_dqdsst(tmp_sst.data,tmp_T2.data,1.2,(tmp_u**2+tmp_v**2)**(1/2),\
                                    tmp_qair.data/1000)
    
            posi_var,nega_var = tmp_dqdsst[:,posi_co],tmp_dqdsst[:,nega_co]
            MY_Var[N]=np.concatenate([posi_var,nega_var],axis=1) 
            N+=1
print('\n!!! END Total Elapsed time : '+str((time.time()-stT)/60)[:4]+'min !!!')

MY_Var[MY_Var<-10**5]=-1e+20

# A=np.concatenate([posi_var,nega_var],axis=1) 
# A[A<-10*5]=-1e+20
# plt.pcolor(A); plt.colorbar()

stT=time.time()
create_dqdsst_nc(save_path+save_name,LON_new,LAT_new,My_time_num,My_time_ref,MY_Var)
print('\n!!! END Total Elapsed time : '+str((time.time()-stT)/60)[:4]+'min !!!')

def create_dqdsst_nc(nc_save_name,LON,LAT,TIME,Ref_time,values):
    
    ncfile = Dataset(nc_save_name,mode='w',format='NETCDF4')

    ncfile.createDimension('lat', len(LAT))
    ncfile.createDimension('lon', len(LON))
    ncfile.createDimension('sst_time',len(TIME))
    
    ncfile.title='ERA5 daily dQdSST'
    
    lat = ncfile.createVariable('lat', np.float32, ('lat',),zlib=False)
    lat.units = 'degrees_north'
    lon = ncfile.createVariable('lon', np.float32, ('lon',),zlib=False)
    lon.units = 'degrees_east'
    time = ncfile.createVariable('sst_time', np.float64, ('sst_time',),zlib=False)
    time.units=Ref_time
    time.field='sst_time, scalar, series'

    dqdsst = ncfile.createVariable('dQdSST',np.float64,('sst_time','lat','lon'),\
                                zlib=False, fill_value=1e+20) # note: unlimited dimension is leftmost
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











































