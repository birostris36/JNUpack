import sys
sys.path.append('C:/Users/shjo/Bridge/JNUpack/SO/libs/')
import matplotlib as mpl
# mpl.use('agg')
import os
import numpy as np
import haversine
from netCDF4 import Dataset
import pickle
from copy import deepcopy
import xarray as xr
from myTrend import myfitting2d_sttcs,myfitting1d_sttcs
from myPlot import  figmaster,myClrbr, dta_colr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata 
from haversine import haversine
import warnings
warnings.filterwarnings('ignore')

pthMO='J:/tmp_proc/Obs/'

myMDOB=[pthMO+i for i in os.listdir(pthMO) if i.endswith('.nc')]

score = lambda x: x-360 if x > 180 else x

for i in myMDOB:
    print(i)
    tmp=Dataset(i,'a')
    lon,lat=tmp['lon'][:],tmp['lat'][:]
    ### ==================================================
    lon_dist,lat_dist=[],[]
    for i in lat:
        lon_dist.append(haversine( (i,lon[0]),(i,lon[1]) ))
    for i in lon:
        lat_dist.append(haversine( (lat[0],score(i)),(lon[1],score(i)) ))
    lon_dist=np.array(lon_dist)
    lat_dist=np.array(lat_dist)
    ### ==================================================
    lon_dist_=deepcopy(lon_dist).reshape(len(lat),1)
    i=0
    while i < len(lon)-1:
        i+=1
        lon_dist_=np.concatenate([lon_dist_,lon_dist.reshape(len(lat),1)],axis=1)
        
    lat_dist_=deepcopy(lat_dist).reshape(1,len(lon))
    i=0
    while i < len(lat)-1:
        i+=1
        lat_dist_=np.concatenate([lat_dist_,lat_dist.reshape(1,len(lon))],axis=0)
    ### ==================================================
    Area=lon_dist_*lon_dist_
    
    myArea = tmp.createVariable('Area',np.float64,('lat','lon'),compression='zlib') #
    myArea.units = 'km^2' 
    myArea.long_name = 'Area' 
    myArea.coordinates = "lat, lon"
    myArea[:]=Area
    tmp.close()
    
    