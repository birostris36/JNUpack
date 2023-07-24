# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 11:46:09 2023

@author: shjo9
"""

import xarray as xr


pth1='D:/HEAT/DATA/'
pth2='E:/_data/MyOHC/'


GECCO_OHC_700=xr.open_dataset(pth2+'GECCO_OHC_SO_c14_700m_1980_2018.nc')\
    .loc[dict(lat=slice(-70,-50),lon=slice(180,270),time=slice('1980-01','2018-12') )].mean(dim=['lat','lon']).OHC
EN4_OHC_700=xr.open_dataset(pth1+'EN4_OHC_GLOBAL_c14_700m_1980_2023.nc')\
    .loc[dict(lat=slice(-70,-50),lon=slice(180,270),time=slice('1980-01','2018-12') )].mean(dim=['lat','lon']).OHC

GECCO_OHC_7002000=xr.open_dataset(pth2+'GECCO_OHC_SO_c14_700m2000m_1980_2018.nc')\
    .loc[dict(lat=slice(-70,-50),lon=slice(180,270),time=slice('1980-01','2018-12') )].mean(dim=['lat','lon']).OHC
EN4_OHC_7002000=xr.open_dataset(pth1+'EN4_OHC_GLOBAL_c14_700m2000m_1980_2023.nc')\
    .loc[dict(lat=slice(-70,-50),lon=slice(180,270),time=slice('1980-01','2018-12') )].mean(dim=['lat','lon']).OHC

GECCO_OHC_700_Y=GECCO_OHC_700.rolling(time=24,center=12).mean()
EN4_OHC_700_Y=EN4_OHC_700.rolling(time=24,center=12).mean()
GECCO_OHC_7002000_Y=GECCO_OHC_7002000.rolling(time=12,center=12).mean()
EN4_OHC_7002000_Y=EN4_OHC_7002000.rolling(time=12,center=12).mean()



EN4_OHC_700_Y.plot()
GECCO_OHC_700_Y.plot()





EN4_OHC_7002000_Y.plot()

GECCO_OHC_7002000_Y.plot()










