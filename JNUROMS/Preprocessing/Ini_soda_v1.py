# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 20:28:14 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

PKG_path = 'D:/OneDrive/JNUpack/JNUROMS/'
import sys 
sys.path.append(PKG_path)
import Tools.JNUROMS as jr
from Tools.JNU_create import create_ini
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
from scipy.interpolate import griddata
import datetime as dt
from tqdm import tqdm
import matplotlib.pyplot as plt

My_Ini='G:/MODEL_DATA/Ini/Ini_Q0_Rtopo30S_S_soda_64250701_198001.nc'
My_Grd='G:/MODEL_DATA/GRD/Grd_Q0_Rtopo30S_Smooth.nc'
OGCM_PATH='G:/SODA/'

Ini_title='Ini_Q0_Rtopo30S_S_soda_64250701_198001'
# My Variables
MyVar={'Layer_N':64,'Vtransform':2,'Vstretching':5,\
       'Theta_s':7,'Theta_b':.1,'Tcline':450,'hmin':50}
conserv=1
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
t_rng=['1980-01','1980-01']
My_time_ref='days since 1970-1-1 00:00:00'
OGCM_TIMES=MFDataset(OGCM_PATH+'*.nc')[tmp_time_var]
TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),28)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
tmp_y,tmp_m=int(t_rng[0].split('-')[0]),int(t_rng[0].split('-')[-1])
tmp_dif=date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-date2num(dt.datetime(tmp_y,tmp_m,1),My_time_ref)
Ini_time_num=((date2num(dt.datetime(tmp_y,tmp_m,1),TIME_UNIT)-tmp_dif)+16)*86400
print('!!! Ini_time + 16d !!!')
# Make an Inifile
create_ini(My_Ini,mask,topo,MyVar,Ini_time_num,Ini_title,ncFormat='NETCDF4')

# Get OGCM data for initial
OGCM_Data={}#,OGCM_Mask={}
for i in ['u','v','temp','salt','zeta','ubar','vbar']:
    print('!!! Data processing : '+i+' !!!')
    if i == 'zeta':
        tmp_data=np.squeeze(MFDataset(OGCMS)[OGCMVar[i]][TIMES_co,latO_co,lonO_co])
        tmp_mask=np.invert(tmp_data.mask)
    elif i=='ubar':
        tmp_u=np.squeeze(MFDataset(OGCMS)[OGCMVar['u']][TIMES_co,:,latO_co,lonO_co])
        tmp_mask_=np.invert(tmp_u.mask)[:,:,:]
        
        tmp_u[tmp_u.mask]=0
        
        du=np.zeros([tmp_u.shape[1],tmp_u.shape[2]])
        zu=np.zeros_like(du)
        dz=np.gradient(-depthO)
        for n in range(len(depthO)):
            du=du+dz[n]*tmp_u[n,:,:].data
            zu=zu+dz[n]*tmp_mask_[n,:,:]
        DATA=du/zu
        # DATA[DATA==0]=np.nan
        tmp_mask=tmp_mask_[0,:,:]
        tmp_data=DATA
        
    elif i=='vbar':
        tmp_v=np.squeeze(MFDataset(OGCMS)[OGCMVar['v']][TIMES_co,:,latO_co,lonO_co])
        tmp_mask_=np.invert(tmp_v.mask)[:,:,:]
        
        tmp_v[tmp_v.mask]=0
        
        dv=np.zeros([tmp_v.shape[1],tmp_v.shape[2]])
        zv=np.zeros_like(dv)
        dz=np.gradient(-depthO)
        for n in range(len(depthO)):
            dv=dv+dz[n]*tmp_v[n,:,:].data
            zv=zv+dz[n]*tmp_mask_[n,:,:]
        DATA=dv/zv
        # DATA[DATA==0]=np.nan
        tmp_mask=tmp_mask_[0,:,:]
        tmp_data=DATA

    else:
        tmp_data=np.squeeze(MFDataset(OGCMS)[OGCMVar[i]][TIMES_co,:,latO_co,lonO_co])
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
        for j,k in tqdm(zip(tmp_data,tmp_mask)):
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
            
# Process vector elements
def rho2u_3d(var):
    N,Mp,Lp=var.shape
    L=Lp-1
    var_u=0.5*(var[:,:,:L]+var[:,:,1:Lp])
    return var_u
def rho2v_3d(var):
    N,Mp,Lp=var.shape
    M=Mp-1
    var_v=0.5*(var[:,:M,:]+var[:,1:Mp,:])
    return var_v
def rho2u_2d(var):
    Mp,Lp=var.shape
    L=Lp-1
    var_u=0.5*(var[:,:L]+var[:,1:Lp])
    return var_u
def rho2v_2d(var):
    Mp,Lp=var.shape
    M=Mp-1
    var_v=0.5*(var[:M,:]+var[1:Mp,:])
    return var_v

u=rho2u_3d(OGCM_Data['u']*cosa+OGCM_Data['v']*sina)
v=rho2v_3d(OGCM_Data['v']*cosa-OGCM_Data['u']*sina)

ubar=rho2u_2d(OGCM_Data['ubar']*cosa+OGCM_Data['vbar']*sina)
vbar=rho2v_2d(OGCM_Data['vbar']*cosa-OGCM_Data['ubar']*sina)

'''
# Calculates barotropic uv
dz=np.diff(depthO)
u_tmp=u[:-1,:,:]+u[1:,:,:]
v_tmp=v[:-1,:,:]+v[1:,:,:]
dz4uv=np.zeros([len(depthO)-1,atG,onG])
for ii in range(atG):
    for jj in range(onG):
        dz4uv[:,ii,jj]=dz
        
# Barotropic velocities2
# ubar=np.sum(u_tmp*dz4uv[:,:,1:],axis=0)/depthO[-1]
# vbar=np.sum(v_tmp*dz4uv[:,1:,:],axis=0)/depthO[-1]
'''

# Process ROMS Vertical grid
Z=np.zeros(len(depthO)+2)
Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000

Rzeta=OGCM_Data['zeta']
zr=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
         MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
             1,topo,Rzeta);
zu=rho2u_3d(zr);
zv=rho2v_3d(zr);
zw=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
         MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
             5,topo,Rzeta);
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

print('!!! Transformming z --> sigma... !!!')
u=jr.ztosigma(np.flip(u1,axis=0),zu,np.flipud(Z));
v=jr.ztosigma(np.flip(v1,axis=0),zv,np.flipud(Z));
temp=jr.ztosigma(np.flip(temp,axis=0),zr,np.flipud(Z));
salt=jr.ztosigma(np.flip(salt,axis=0),zr,np.flipud(Z));

if conserv==1:
    u = u - np.sum(u*dzu,axis=0)/np.sum(dzu,axis=0) + ubar
    v = v - np.sum(v*dzv,axis=0)/np.sum(dzv,axis=0) + vbar

# Barotropic velocities2
ubar_=np.sum(u*dzu,axis=0)/np.sum(dzu,axis=0)
vbar_=np.sum(v*dzv,axis=0)/np.sum(dzv,axis=0)

ncI=Dataset(My_Ini,mode='a')
ncI['zeta'][:]=Rzeta
# ncI['SSH'][:]=Rzeta
ncI['temp'][:]=temp
ncI['salt'][:]=salt
ncI['u'][:]=u
ncI['v'][:]=v
ncI['ubar'][:]=ubar_
ncI['vbar'][:]=vbar_
ncI.close()




