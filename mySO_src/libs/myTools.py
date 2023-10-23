import time
import os
from shutil import copyfile
import gsw.density as gsw_d
import gsw.conversions as gsw_c
from gsw._wrapped_ufuncs import cp_t_exact
import sys
from haversine import haversine
from scipy.interpolate import interp2d, griddata, LinearNDInterpolator, NearestNDInterpolator
import numpy as np

def myInfo(loc,wpth):
    if not wpth.endswith('/'):
        raise
    with open(loc,'r') as mC:
        myCode=mC.read()
    with open(wpth+"/myLog.txt", 'w') as f:
        info='Ran time:\t'+time.strftime('%c')+'\nScript loc:\t'+os.getcwd().replace('\\','/')\
            +'\nScript name:\t'+loc.replace('\\','/')+'\nwpth:\t\t'+ wpth+'\n\n\n\n\n'+\
                'Code :\n\n\n'+myCode
        f.write(info)
    
    # copyfile(loc,wpth+loc.split("\\")[-1].replace('py','txt'))


def latlon_dist(LAT,LON):
    '''
    LAT,LON-->1d numpy.array
    '''
    score = lambda x: x-360 if x > 180 else x

    at,on=LAT.shape,LON.shape
        
    lat_dist=np.zeros([at,on])
    ni,nj=-1,-1
    while ni<at-2:
        ni+=1
        while nj < on-1:
            nj+=1
            lat_dist[ni,nj]=haversine( ( LAT[ni],score(LON[nj]) ), ( LAT[ni+1],score(LON[nj]) ) )
        nj=-1
    lat_dist[-1,:]=lat_dist[-2,:]

    lon_dist=np.zeros([at,on])
    ni,nj=-1,-1
    while ni<at-1:
        ni+=1
        while nj < on-2:
            nj+=1
            lon_dist[ni,nj]=haversine( ( LAT[ni],score(LON[nj]) ), ( LAT[ni],score(LON[nj+1]) ) )
        nj=-1
    lon_dist[:,-1]=lon_dist[:,-2]
        
    return lat_dist, lon_dist
    
    
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
    CP=(CP[:,1:,:,:]+CP[:,:-1,:,:])/2
    CT=(CT[:,1:,:,:]+CT[:,:-1,:,:])/2
    rho=(rho[:,1:,:,:]+rho[:,:-1,:,:])/2

    DZ=np.diff(DEPTH)
    
    T,D,A,O=CT.shape

    dz=np.tile(np.tile(np.tile(DZ.reshape([D,1]), A )\
                    .reshape([D,A,-1]), O), (T,1,1,1) )
    
    OHC_=CP*CT*rho*dz
    
    # Integrates from ref depth (2000m) 
    OHC=np.sum(OHC_,axis=1)

    return OHC
