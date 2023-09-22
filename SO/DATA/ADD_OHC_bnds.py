import matplotlib.pyplot as plt
import gsw.density as gsw_d
import gsw.conversions as gsw_c
from gsw._wrapped_ufuncs import cp_t_exact
import xarray as xr
from netCDF4 import Dataset
import os
from numpy import ma
import numpy as np

npth='J:/ALL/myISHII_198001_201912_tsh.nc'

DATA=xr.open_dataset(npth).loc[dict(depth=slice(0,750))]

DEPTH,LAT,LON=DATA.depth,DATA.lat,DATA.lon 
SA,TEMP=DATA.salt,DATA.temp
DEPTH_BND=DATA['depth_bnds']

def calc_ohc_bnds(DEPTH,DEPTH_BND,LAT,LON,SA,TEMP,Abs_salt=True):
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
    
    dz=DEPTH_BND[:,1]-DEPTH_BND[:,0]
    OHC_=CP*CT*rho*dz
    
    # Integrates from ref depth (2000m) 
    OHC=OHC_.sum(dim='depth',skipna=False).values
    return OHC

OHC=calc_ohc_bnds(DEPTH,LAT,LON,SA,TEMP,Abs_salt=True)
mask=(OHC!=OHC).data
OHC=ma.array(OHC,mask=mask)

DATA.close()
del DATA

NC=Dataset(npth,'a')
myOHC = NC.createVariable('OHC_700',np.float64,('time','lat','lon'),compression='zlib') #
myOHC.units = 'Jm^-2' 
myOHC.long_name = 'ISHII OHC 700' 
myOHC.coordinates = "time, lat, lon"
myOHC[:]=OHC
NC.close()

