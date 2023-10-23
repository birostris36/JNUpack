import xarray as xr
import os
import numpy as np
from netCDF4 import Dataset,num2date
import sys
import gsw.density as gsw_d
import gsw.conversions as gsw_c
from gsw._wrapped_ufuncs import cp_t_exact
# sys.path.append('C:/Users/shjo/Bridge/JNUpack/mySO_src/libs/')
# from myTools import NearesND,calc_ohc

npth='/home/shjo/_data/myARMOR3D_SO_199301_202312_ts.nc'

### ===========================================
Sample=Dataset(npth,'r')
TIME=Sample['time'][:]
LAT,LON=Sample['lat'][:],Sample['lon'][:]
DEPTH=Sample['depth'][:]



T,D,AT,ON=len(TIME),len(DEPTH),len(LAT),len(LON)

def calc_ohc(DEPTH,LAT,LON,SA,TEMP,Abs_salt=True):
    '''
    All values --> xarray.dataarray 
    '''
    if Abs_salt:
        pass
    else:
        #Practical Salinity --> Absolute salinity
        SA = gsw_c.SA_from_SP(SA,DEPTH,LON,LAT) # [SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat)
    
    # Potential temp --> Conservative temp
    CT=gsw_c.CT_from_pt(SA,TEMP) #CT = gsw_CT_from_pt(SA,pt)
    rho = gsw_d.rho(SA,CT,DEPTH)

    # Potential temperature --> In-situ temperature
    t=gsw_c.t_from_CT(SA,CT,DEPTH)

    # Heat capacity
    CP = cp_t_exact(SA,t,DEPTH)
    
    CP,CT,rho=CP.values,CT.values,rho.values
    CP=(CP[1:,:,:]+CP[:-1,:,:])/2
    CT=(CT[1:,:,:]+CT[:-1,:,:])/2
    rho=(rho[1:,:,:]+rho[:-1,:,:])/2

    DZ=np.diff(DEPTH)
    
    D,A,O=CT.shape

    dz=np.tile(np.tile(DZ.reshape([D,1]), A )\
                    .reshape([D,A,-1]), O)

    OHC_=CP*CT*rho*dz
    
    # Integrates from ref depth (2000m) 
    OHC=np.sum(OHC_,axis=0)

    return OHC

try: 
    Sample['OHC700']
    Sample.close()
except:
    Sample.close()
    print('!!! zeros !!!')
    OHC700=np.zeros([T,AT,ON])
    OHC2000=np.zeros_like(OHC700)
    OHC7002000=np.zeros_like(OHC700)
    
    NC=Dataset(npth,'a')
    myOHC700 = NC.createVariable('OHC700',np.float64,('time','lat','lon'),compression='zlib') #
    myOHC7002000 = NC.createVariable('OHC7002000',np.float64,('time','lat','lon'),compression='zlib') #
    myOHC2000 = NC.createVariable('OHC2000',np.float64,('time','lat','lon'),compression='zlib') #

    myOHC700.units = 'Jm^-2' 
    myOHC7002000.units = 'Jm^-2' 
    myOHC2000.units = 'Jm^-2' 
    
    myOHC700.long_name = 'ARMOR3D OHC 700'
    myOHC7002000.long_name = 'ARMOR3D OHC 7002000'
    myOHC2000.long_name = 'ARMOR3D OHC 2000'
       
    myOHC700.coordinates = "time, lat, lon"
    myOHC7002000.coordinates = "time, lat, lon"
    myOHC2000.coordinates = "time, lat, lon"

    myOHC700[:]=OHC700
    myOHC7002000[:]=OHC7002000
    myOHC2000[:]=OHC2000

    NC.close()

### Calc OHC ========================================
print('!!! Open nc ... !!!')
NC=Dataset(npth,'a')
XNC=xr.open_dataset(npth)
XDEPTH_700=XNC.depth.loc[dict(depth=slice(0,700))]
XDEPTH_7002000=XNC.depth.loc[dict(depth=slice(700,2000))]
XDEPTH_2000=XNC.depth.loc[dict(depth=slice(0,2000))]
XLAT,XLON=XNC.lat,XNC.lon

XTMEP700=XNC.temp.loc[dict(depth=slice(0,700))]
XTMEP7002=XNC.temp.loc[dict(depth=slice(700,2000))]
XTMEP2000=XNC.temp.loc[dict(depth=slice(0,2000))]

XSALT700=XNC.salt.loc[dict(depth=slice(0,700))]
XSALT7002=XNC.salt.loc[dict(depth=slice(700,2000))]
XSALT2000=XNC.salt.loc[dict(depth=slice(0,2000))]


for tt in range(T):

    OHC700    =calc_ohc(XDEPTH_700    ,XLAT,XLON,XSALT700[tt] ,XTMEP700[tt] ,Abs_salt=False)
    OHC7002000=calc_ohc(XDEPTH_7002000,XLAT,XLON,XSALT7002[tt],XTMEP7002[tt],Abs_salt=False)
    OHC2000   =calc_ohc(XDEPTH_2000   ,XLAT,XLON,XSALT2000[tt],XTMEP2000[tt],Abs_salt=False)

    NC['OHC700'][tt,:,:]=OHC700
    NC['OHC7002000'][tt,:,:]=OHC7002000
    NC['OHC2000'][tt,:,:]=OHC2000

    if tt%30==0:
        print('!!! '+f'{tt/T*100:0.1f}'+' !!!')
NC.close()











    
