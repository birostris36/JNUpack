# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 16:42:59 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""


PKG_path = 'D:/Working_hub/OneDrive/base142/Factory/MantaROMS/src_d/'
import sys 
sys.path.append(PKG_path)
import JNUROMS as jr
from JNU_create import create_bry
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
from scipy.interpolate import griddata
import datetime as dt
import multiprocessing as mp
from tqdm import tqdm
import time
import datetime as dt


My_Bry='D:/Working_hub/OneDrive/base142/Factory/MantaROMS/Test_nc/test_ini.nc'
My_Grd='D:/Working_hub/OneDrive/base142/Warehouse01/Grd_v2_10d_24S.nc'
OGCM_PATH='D:/Working_hub/DATA_dep/SODA/'

Bry_title='test_bry'
# My Variables
MyVar={'Layer_N':36,'Vtransform':1,'Vstretching':1,\
       'Theta_s':5,'Theta_b':.4,'Tcline':50,'hmin':50}

# OGCM Variables
OGCMVar={'lon_rho':'xt_ocean','lat_rho':'yt_ocean','depth':'st_ocean','time':'time',\
         'lon_u':'xu_ocean','lat_u':'yu_ocean','lon_v':'xu_ocean','lat_v':'yu_ocean',
         'temp':'temp','salt':'salt','u':'u','v':'v','zeta':'ssh'}
    
OGCMS=[OGCM_PATH+'/'+i for i in os.listdir(OGCM_PATH) if i.endswith('.nc')]

# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
ncG.close()

atG,onG=lonG.shape
cosa,sina=np.cos(angle),np.sin(angle)

# Get OGCM Grid info
ncO=Dataset(OGCMS[0])
lonO,latO=ncO[OGCMVar['lon_rho']][:],ncO[OGCMVar['lat_rho']][:];
depthO=ncO[OGCMVar['depth']][:]
ncO.close()

# Get OGCM lon lat coordinates for slicing





# =============================================================================
# Process Times
tmp_time_var='time'
t_rng=['1980-01','2017-12']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=MFDataset(OGCM_PATH+'*.nc')[tmp_time_var]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
Bry_time_time=num2date(OGCM_TIMES[TIMES_co]-tmp_dif,My_time_ref)
Bry_time_num=OGCM_TIMES[TIMES_co]-tmp_dif

# Make an Inifile
# create_bry(My_Bry,mask,topo,MyVar,Bry_time_num,Bry_title,ncFormat='NETCDF3_CLASSIC')

latO_co=np.where( (latO>=latG[-2,0]) & (latO<=latG[-1,0]) )[0]
lonO_co=np.where( (lonO>=lonG[0,0]) & (lonO<=lonG[0,-1]) )[0]

lonO_s_m,latO_s_m=np.meshgrid(lonO[lonO_co],latO[latO_co])
Tt=len(TIMES_co)
# Get OGCM data for initial
OGCM_Data={}#,OGCM_Mask={}

# trans Coordinate for SODA (1 file --> 12month)
TIMES_co_set=list(set([int(i/12) for i in TIMES_co]))

def main(child_co):
    
    global parent_co, OGCMS, Counter
    

    for i in ['zeta','u','v','temp','salt']:
        print('!!! Data processing : '+i+' !!!')
        if i=='zeta':
            data=np.zeros([Tt,onG])
            for t in tqdm(range(Tt)):      
                tmp_data=np.squeeze(MFDataset(OGCMS[parent_co])[OGCMVar[i]][t,latO_co,lonO_co])
                tmp_mask=np.invert(tmp_data.mask)
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                              tmp_data[tmp_mask].flatten(),(lonG[-1:,:].flatten(),latG[-1:,:].flatten()),'nearest')
                data_=data_.reshape(lonG[-1:,:].shape)
                data[t]=data_
        else:
            data,n=np.zeros([Tt,len(depthO),onG]),0
            for t in tqdm(range(Tt)):  
                tmp_data=np.squeeze(Dataset(OGCMS[parent_co])[OGCMVar[i]][t,:,latO_co,lonO_co])
                tmp_mask=np.invert(tmp_data.mask)
                n=0;
                for j,k in zip(tmp_data,tmp_mask):
                    data_=griddata((lonO_s_m[k].flatten(),latO_s_m[k].flatten()),\
                                  j[k].flatten(),(lonG[-1:,:].flatten(),latG[-1:,:].flatten()),'nearest')
                    data[t,n]=data_.reshape(lonG[-1:,:].shape)
                    n+=1  
        OGCM_Data[i]=data
                
    # Process vector elements
    def rho2u_3d(var):
        N,Mp,Lp=var.shape
        L=Lp-1
        var_u=0.5*(var[:,:,:L]+var[:,:,1:Lp])
        return var_u
    def rho2v_3d(var):
        N,Mp,Lp=var.shape
        M=Mp-1
        var_v=0.5*(var[:,:M,:]+var[:,1:Lp,:])
        return var_v
    
    u=rho2u_3d(OGCM_Data['u']*cosa+OGCM_Data['v']*sina)
    v=rho2v_3d(OGCM_Data['v']*cosa-OGCM_Data['u']*sina)
    
    # Calculates barotropic uv
    dz=np.diff(depthO)
    u_tmp=u[:-1,:,:]+u[1:,:,:]
    v_tmp=v[:-1,:,:]+v[1:,:,:]
    dz4uv=np.zeros([len(depthO)-1,atG,onG])
    for ii in range(atG):
        for jj in range(onG):
            dz4uv[:,ii,jj]=dz
    # Barotropic velocities1
    ubar=np.sum(u_tmp*dz4uv[:,:,1:],axis=0)/depthO[-1]
    vbar=np.sum(v_tmp*dz4uv[:,1:,:],axis=0)/depthO[-1]
    
    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000
    
    Rzeta=OGCM_Data['zeta']
    zr=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 1,topo,Rzeta,1);
        
    zu=rho2u_3d(zr);
    zv=rho2v_3d(zr);
    zw=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 5,topo,Rzeta,1);
        
    dzr=zw[1:,:,:]-zw[:-1,:,:]
    dzu=rho2u_3d(dzr);
    dzv=rho2v_3d(dzr);
    
    # Add a level on top and bottom with no-gradient
    temp,salt=OGCM_Data['temp'],OGCM_Data['salt']
    
    u1=np.vstack((np.expand_dims(u[0,:,:],axis=0)\
                  ,u,np.expand_dims(u[-1,:,:],axis=0)))
    v1=np.vstack((np.expand_dims(v[0,:,:],axis=0)\
                  ,v,np.expand_dims(v[-1,:,:],axis=0)))
    temp=np.vstack((np.expand_dims(temp[0,:,:],axis=0)\
                  ,temp,np.expand_dims(temp[-1,:,:],axis=0)))
    salt=np.vstack((np.expand_dims(salt[0,:,:],axis=0)\
                  ,salt,np.expand_dims(salt[-1,:,:],axis=0)))
    
    u=jr.ztosigma(np.flip(u1,axis=0),zu,np.flipud(Z));
    v=jr.ztosigma(np.flip(v1,axis=0),zv,np.flipud(Z));
    temp=jr.ztosigma(np.flip(temp,axis=0),zr,np.flipud(Z));
    salt=jr.ztosigma(np.flip(salt,axis=0),zr,np.flipud(Z));
    
    # Barotropic velocities2
    ubar_=np.sum(u*dzu,axis=0)/np.sum(dzu,axis=0)
    vbar_=np.sum(v*dzv,axis=0)/np.sum(dzv,axis=0)
    
    ncI=Dataset(My_Bry,mode='a')
    ncI['zeta'][:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncI['temp'][:]=temp
    ncI['salt'][:]=salt
    ncI['u'][:]=u
    ncI['v'][:]=v
    ncI['ubar'][:]=ubar_
    ncI['vbar'][:]=vbar_
    ncI.close()


# =============================================================================
# Main
# =============================================================================
stT=time.time()
for parent_co in TIMES_co_set:
    
    with mp.Pool(12) as MyPool:
        results=MyPool.map(main,range(12))

















