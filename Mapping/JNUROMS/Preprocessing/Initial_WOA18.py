# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:44:08 2023

@author: shjo9
"""


PKG_path = 'D:/OneDrive/JNUpack/JNUROMS/'
import sys 
sys.path.append(PKG_path)
import Tools.JNUROMS as jr
from Tools.JNU_create import create_ini_WOA
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
from scipy.interpolate import griddata
import datetime as dt
from tqdm import tqdm
import matplotlib.pyplot as plt

My_Ini='G:/MODEL_DATA/Ini/TEST_INI_WOA.nc'
My_Grd='G:/MODEL_DATA/Grd/Grd_SO_05d_sponge.nc'
OGCM_PATH='G:/WOA18/'

TEMP_PATH,SALT_PATH=OGCM_PATH+'temp_M025/',OGCM_PATH+'salt_M025/'
TEMP_ann_PATH,SALT_ann_PATH=OGCM_PATH+'temp_M025/',OGCM_PATH+'salt_M025/'

temp_M_WOAs=TEMP_PATH+'woa18_decav_t01_04.nc'
salt_M_WOAs=SALT_PATH+'woa18_decav_s01_04.nc'

temp_A_WOAs=TEMP_PATH+'woa18_decav_t00_04.nc'
salt_A_WOAs=SALT_PATH+'woa18_decav_s00_04.nc'


Ini_title=''
# My Variables
MyVar={'Layer_N':50,'Vtransform':2,'Vstretching':5,\
       'Theta_s':7,'Theta_b':.4,'Tcline':300,'hmin':50}

# OGCM Variables
OGCMVar={'lon_rho':'lon','lat_rho':'lat','depth':'depth','time':'time',\
         'lon_u':'lon','lat_u':'lat','lon_v':'lon','lat_v':'lat'}

# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
ncG.close()

atG,onG=lonG.shape
cosa,sina=np.cos(angle),np.sin(angle)

# Get OGCM Grid info
ncO=Dataset(TEMP_PATH+'woa18_decav_t01_04.nc')
lonO,latO=ncO[OGCMVar['lon_rho']][:],ncO[OGCMVar['lat_rho']][:];
depthO=ncO[OGCMVar['depth']][:]
ncO.close()

# Make an Inifile
create_ini_WOA(My_Ini,mask,topo,MyVar,Ini_title,ncFormat='NETCDF4')

# Slice lat northern bound
latO_co=np.where( (latO>=-85.5) & (latO<=-23) )[0]

# Change data to Pacific centered 
p_co,n_co=np.where(lonO>=0)[0],np.where(lonO<0)[0]

lonO_re=np.concatenate([lonO[p_co],lonO[n_co]+360],axis=0)

lonO_re_m,latO_m=np.meshgrid(lonO_re,latO[latO_co])

# =============================================================================
# 
# =============================================================================
depth_ann=Dataset(salt_A_WOAs)['depth'][:]
depth_ann_co=np.where(depth_ann>1500)[0]
depth_ann_re=depth_ann[depth_ann_co]
salt_A=np.squeeze(Dataset(salt_A_WOAs)['s_an'][0,depth_ann_co,latO_co,:])
temp_A=np.squeeze(Dataset(temp_A_WOAs)['t_an'][0,depth_ann_co,latO_co,:])
depthO_re=np.concatenate([depthO, depth_ann_re],axis=0)


for i_temp, j_salt in zip([temp_M_WOAs],[salt_M_WOAs]):
    
    tmp_temp=np.squeeze(Dataset(i_temp)['t_an'][0,:,latO_co,:])
    tmp_salt=np.squeeze(Dataset(j_salt)['s_an'][0,:,latO_co,:])
    
    tmp_temp_re=np.concatenate( [tmp_temp,temp_A ], axis=0)
    tmp_salt_re=np.concatenate( [tmp_salt,salt_A ], axis=0)
    
    tmp_temp3d=np.concatenate( (tmp_temp_re[:,:,p_co], tmp_temp_re[:,:,n_co]),axis=2)
    tmp_salt3d=np.concatenate( (tmp_salt_re[:,:,p_co], tmp_salt_re[:,:,n_co]),axis=2)
    
    tmp_t_mask3d=np.concatenate( (np.invert(tmp_temp.mask)[:,:,p_co], np.invert(tmp_temp.mask)[:,:,n_co]),axis=2)
    tmp_s_mask3d=np.concatenate( (np.invert(tmp_salt.mask)[:,:,p_co], np.invert(tmp_salt.mask)[:,:,n_co]),axis=2)

    temp3d,n=np.zeros([len(depthO_re),atG,onG]),0
    salt3d=np.zeros_like(temp3d)

    for i_temp2d, j_salt2d, k_mask2d, l_mask2d in tqdm(zip(tmp_temp3d,tmp_salt3d,tmp_t_mask3d,tmp_s_mask3d)): 
        # 3D interpolation for temperature
        temp_=griddata((lonO_re_m[k_mask2d].flatten(),latO_m[k_mask2d].flatten()),\
                      i_temp2d[k_mask2d].flatten(),(lonO_re_m.flatten(),latO_m.flatten()),'nearest')
        temp_ex=temp_.reshape(latO_m.shape)
        temp_=griddata((lonO_re_m.flatten(),latO_m.flatten()),\
                      temp_ex.flatten(),(lonG.flatten(),latG.flatten()),'cubic')
        temp3d[n]=temp_.reshape(latG.shape)
        
        # 3D interpolation for temperature
        salt_=griddata((lonO_re_m[l_mask2d].flatten(),latO_m[l_mask2d].flatten()),\
                      j_salt2d[l_mask2d].flatten(),(lonO_re_m.flatten(),latO_m.flatten()),'nearest')
        salt_ex=salt_.reshape(latO_m.shape)
        salt_=griddata((lonO_re_m.flatten(),latO_m.flatten()),\
                      salt_ex.flatten(),(lonG.flatten(),latG.flatten()),'cubic')
        salt3d[n]=salt_.reshape(latG.shape)
        # Counter n
        n+=1
# =============================================================================            
# Process ROMS Vertical grid
Z=np.zeros(len(depthO_re)+2)
Z[0]=100;Z[1:-1]=-depthO_re;Z[-1]=-100000

Rzeta=np.zeros([atG,onG])
zr=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
         MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
             1,topo,Rzeta);
zw=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
         MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
             5,topo,Rzeta);
dzr=zw[1:,:,:]-zw[:-1,:,:]

# Add a level on top and bottom with no-gradient
temp=np.vstack((np.expand_dims(temp3d[0,:,:],axis=0)\
              ,temp3d,np.expand_dims(temp3d[-1,:,:],axis=0)))
salt=np.vstack((np.expand_dims(salt3d[0,:,:],axis=0)\
              ,salt3d,np.expand_dims(salt3d[-1,:,:],axis=0)))

print('!!! Transformming z --> sigma... !!!')
temp=jr.ztosigma(np.flip(temp,axis=0),zr,np.flipud(Z));
salt=jr.ztosigma(np.flip(salt,axis=0),zr,np.flipud(Z));

# Barotropic velocities2

ncI=Dataset(My_Ini,mode='a')
ncI['zeta'][:]=Rzeta
# ncI['SSH'][:]=Rzeta
ncI['temp'][:]=temp
ncI['salt'][:]=salt
ncI.close()




















