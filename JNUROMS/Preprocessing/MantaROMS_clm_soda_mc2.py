#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:38:43 2023

@author: birostris
@Email: birostris36@gmail.com 

@Category:  
@Reference: None
@Descriptions:
"""




PKG_path = '/data2/base132/Factory/tmp_MantaROMS/src_d/'
import sys 
sys.path.append(PKG_path)
import JNUROMS as jr
from JNU_create import create_ini,create_clm_SODA
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
from scipy.interpolate import griddata
import datetime as dt
from tqdm import tqdm
import matplotlib.pyplot as plt
import multiprocessing as mp
import time

My_Grd='/data2/base132/Warehouse01/Grd_SO_05d.nc'
# My_Grd='/data2/base132/Warehouse01/Grd_JYGO.nc'

My_Clm='/data2/base132/Warehouse01/Clm_SODA_1980_2017.nc'
OGCM_PATH='/data2/SODA/'

Clm_title='test_CLM'
# My Variables
MyVar={'Layer_N':50,'Vtransform':2,'Vstretching':5,\
       'Theta_s':7,'Theta_b':0.1,'Tcline':300,'hmin':50}

# OGCM Variables
OGCMVar={'lon_rho':'xt_ocean','lat_rho':'yt_ocean','depth':'st_ocean','time':'time',\
         'lon_u':'xu_ocean','lat_u':'yu_ocean','lon_v':'xu_ocean','lat_v':'yu_ocean',
         'temp':'temp','salt':'salt','u':'u','v':'v','zeta':'ssh'}
    
OGCMS=[OGCM_PATH+'/'+i for i in os.listdir(OGCM_PATH) if i.endswith('.nc')]

# My CPUs
print('MY CPU NUM : ',mp.cpu_count())
print(mp.current_process().name,mp.current_process().pid)

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
lonGmax=np.max(np.abs(np.diff(lonG[0,:])))
latGmax=np.max(np.abs(np.diff(latG[:,0])))

lonOmax=np.max(np.abs(np.diff(lonO[:])))
latOmax=np.max(np.abs(np.diff(latO[:])))

lonEval=np.max([lonO[-1],lonG[0,-1]+lonOmax])
lonSval=np.min([lonO[0],lonG[0,0]+lonOmax])

latEval=np.max([latO[-1],latG[-1,0]+latOmax])
latSval=np.min([latO[0],latG[0,0]+latOmax])

lonO_co=np.where(lonO[(lonO>=lonSval-lonOmax)&(lonO<=lonEval)])[0]

latO_co=np.where(latO[(latO>=latG[0,0]-latOmax)&(latO<=latG[-1,0]+0.5)])[0]
latO_s=latO[latO_co]
lonO_s=lonO[lonO_co]

# =============================================================================
if latG[0,0]<latO_s[0]:
    OGCM_lat_diff=np.diff(latO_s[0:2].data)
    tmp_l=np.arange(latO_s[0]-OGCM_lat_diff,latG[0,0]-OGCM_lat_diff,-OGCM_lat_diff)
    tmp_array=np.full((len(tmp_l), len(lonO_s[:])), False)
    tmp_array_val=np.zeros_like(tmp_array)
    latO_s_=np.concatenate([np.flipud(tmp_l),latO_s[:]],axis=0)
   
# =============================================================================
lonO_s_m,latO_s_m=np.meshgrid(lonO_s,latO_s_)
# =============================================================================
# Process Times
tmp_time_var='time'
t_rng=['1980-01','2017-12']
My_time_ref='days since 1950-1-1 00:00:00'
OGCM_TIMES=MFDataset(OGCM_PATH+'*.nc')[tmp_time_var]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),30)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
OGCM_TIMES_re=num2date(OGCM_TIMES[:]-tmp_dif,My_time_ref)
Clm_TIMES=OGCM_TIMES[TIMES_co]-tmp_dif

# Make an Inifile
st=time.time()
# create_clm_SODA(My_Clm,mask,topo,MyVar,Clm_TIMES,My_time_ref,Clm_title,ncFormat='NETCDF4')
print('\n!!! Elapsed time : '+str((time.time()-st)/60)[:4]+'min !!!')

# trans_Coordinate for SODA (1 file --> 12 month)
TIMES_co_set=list(set([int(i/12) for i in TIMES_co]))


def main(child_co):
    
    global parent_co,OGCMS,Counter 
    
    st=time.time()
    print('\n!!! Processing month (child_co+1) : '+str(child_co+1)+' !!!')
    OGCM_Data={}#,OGCM_Mask={}
    for i in ['zeta','u','v','temp','salt']:
        # print('!!! Data processing : '+i+' !!!')
        if i == 'zeta':
            tmp_data=np.squeeze(Dataset(OGCMS[parent_co])[OGCMVar[i]][child_co,latO_co,lonO_co])
            tmp_mask=np.invert(tmp_data.mask)
        else:
            tmp_data=np.squeeze(Dataset(OGCMS[parent_co])[OGCMVar[i]][child_co,:,latO_co,lonO_co])
            tmp_mask=np.invert(tmp_data.mask)
        ## Another way to create mask
        # mv=ncO[OGCMVar[i]].missing_value
        if len(tmp_data.shape)==2:
            tmp_data=np.concatenate([tmp_array_val,tmp_data],axis=0)
            tmp_mask=np.concatenate([tmp_array,tmp_mask])
            data=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                          tmp_data[tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
            data_ex=data.reshape(latO_s_m.shape)
            data=griddata((lonO_s_m.flatten(),latO_s_m.flatten()),\
                          data_ex.flatten(),(lonG.flatten(),latG.flatten()),'cubic')
            data=data.reshape(lonG.shape)
           
        elif len(tmp_data.shape)==3:
            data,n=np.zeros([len(depthO),atG,onG]),0
            for j,k in zip(tmp_data,tmp_mask):
                tmp_data_j=np.concatenate([tmp_array_val,j],axis=0)
                tmp_mask_k=np.concatenate([tmp_array,k])
                data_=griddata((lonO_s_m[tmp_mask_k].flatten(),latO_s_m[tmp_mask_k].flatten()),\
                              tmp_data_j[tmp_mask_k].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_ex=data_.reshape(latO_s_m.shape)
                data_=griddata((lonO_s_m.flatten(),latO_s_m.flatten()),\
                              data_ex.flatten(),(lonG.flatten(),latG.flatten()),'cubic')
                data[n]=data_.reshape(latG.shape)
                n+=1
        OGCM_Data[i]=data
        
    u=jr.rho2u_3d(OGCM_Data['u']*cosa+OGCM_Data['v']*sina)
    v=jr.rho2v_3d(OGCM_Data['v']*cosa-OGCM_Data['u']*sina)
    
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
                 1,topo,Rzeta);
    zu=jr.rho2u_3d(zr);
    zv=jr.rho2v_3d(zr);
    zw=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 5,topo,Rzeta);
    dzr=zw[1:,:,:]-zw[:-1,:,:]
    dzu=jr.rho2u_3d(dzr);
    dzv=jr.rho2v_3d(dzr);
    
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
    
    time.sleep(child_co*6)
    ncI=Dataset(My_Clm,mode='a')

    ncI['zeta'][parent_co*12+child_co,:,:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncI['temp'][parent_co*12+child_co,:,:,:]=temp
    ncI['salt'][parent_co*12+child_co,:,:,:]=salt
    ncI['u'][parent_co*12+child_co,:,:,:]=u
    ncI['v'][parent_co*12+child_co,:,:,:]=v
    ncI['ubar'][parent_co*12+child_co,:,:]=ubar
    ncI['vbar'][parent_co*12+child_co,:,:]=vbar
      
    ncI.close()
    # =============================================================================
    ET=time.time()-st
    print('\n!!! Elapsed time : '+str(ET/60)[:4]+'min !!!')
    # =============================================================================

TIMES_co_set=TIMES_co_set[25:]
Counter=25
Error_counter=0


stT=time.time()
for parent_co in TIMES_co_set:
    try:
        with mp.Pool(12) as MyPool:
            results=MyPool.map(main,range(12))
        print('\n!!! Counter : '+str(Counter)+' !!!')
        Counter+=1;
    except :
        print('\n!!! Error_count : '+str(Error_counter)+' !!!')
        Error_counter+=1
        with mp.Pool(12) as MyPool:
            results=MyPool.map(main,range(12))
        print('\n!!! Counter : '+str(Counter)+' !!!')
        Counter+=1;
    if Error_counter==20:
        break
print('\n!!! END Total Elapsed time : '+str((time.time()-stT)/60)[:4]+'min !!!')  


 













