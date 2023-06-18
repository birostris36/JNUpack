# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 21:30:42 2023

@author: shjo9
"""
PKG_path = 'D:/OneDrive/JNUpack/JNUROMS/'
import sys 
sys.path.append(PKG_path)
import numpy as np
from netCDF4 import Dataset,MFDataset
import os
import xarray as xr
import dask
import matplotlib as mpl
import cmocean
import Tools.JNUROMS as jr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from scipy.interpolate import interp2d, griddata
from tqdm import tqdm

My_Clm='G:/MODEL_DATA/TEST_CLM_CLM.nc'
My_Grd='G:/MODEL_DATA/GRD/Grd_Q0_Rtopo30S_Smooth.nc'
OGCM_PATH='G:/SODA/'

Ini_title='TEST_CLM_CLM'
# My Variables
MyVar={'Layer_N':64,'Vtransform':2,'Vstretching':5,\
       'Theta_s':7,'Theta_b':.1,'Tcline':450,'hmin':50}
conserv=1
# OGCM Variables
OGCMVar={'lon_rho':'xt_ocean','lat_rho':'yt_ocean','depth':'st_ocean','time':'time',\
         'lon_u':'xu_ocean','lat_u':'yu_ocean','lon_v':'xu_ocean','lat_v':'yu_ocean',
         'temp':'temp','salt':'salt','u':'u','v':'v','zeta':'ssh'}
    
# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
ncG.close()

atG,onG=lonG.shape
cosa,sina=np.cos(angle),np.sin(angle)

# Get OGCM Grid info
ncO=Dataset('G:/MODEL_DATA/'+'SODA_CLM_MEAN.nc')
lonO,latO=ncO[OGCMVar['lon_rho']][:],ncO[OGCMVar['lat_rho']][:];
depthO=ncO[OGCMVar['depth']][:]
thO=ncO[OGCMVar['depth']].shape[0]

# ncO.close()

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
lonO_s_m,latO_s_m=np.meshgrid(lonO[lonO_co],latO[latO_co])


# =============================================================================
if latG[0,0]<latO_s[0]:
    OGCM_lat_diff=np.diff(latO_s[0:2].data)
    tmp_l=np.arange(latO_s[0]-OGCM_lat_diff,latG[0,0]-OGCM_lat_diff,-OGCM_lat_diff)
    tmp_array=np.full((len(tmp_l), len(lonO_s[:])), False)
    tmp_array_val=np.zeros_like(tmp_array)
    latO_s_=np.concatenate([np.flipud(tmp_l),latO_s[:]],axis=0)
    
    
    
OGCM_Data={}#,OGCM_Mask={}
def main(child_co):

    for i in ['u','v','temp','salt','zeta','ubar','vbar']:
        print('!!! Data processing : '+i+' !!!')
        
        if (i=='zeta') or (i=='ubar') or (i=='vbar'):
            data=np.zeros([lonG.shape[0],lonG.shape[-1]])
            if i=='zeta':
                DATA=np.squeeze(ncO[OGCMVar[i]][child_co,latO_co,lonO_co])
                tmp_mask_=np.invert(DATA.mask)
                
            elif i=='ubar':
                tmp_u=np.squeeze(ncO[OGCMVar['u']][child_co,:,latO_co,lonO_co])
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
                tmp_mask_=tmp_mask_[0,:,:]
                
            elif i=='vbar':
                tmp_v=np.squeeze(ncO[OGCMVar['v']][child_co,:,latO_co,lonO_co])
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
                tmp_mask_=tmp_mask_[0,:,:]
                
            tmp_mask=tmp_mask_
            data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                          DATA[tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
            data_=data_.reshape(latO_s_m.shape)
    
            # Interp 4 Grid
            data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:,:].flatten(),latG[:,:].flatten()) ,'cubic' )
            data=data_re_.reshape(lonG[:,:].shape)
            OGCM_Data[i]=data
            
        else:
            data=np.zeros([len(depthO),lonG.shape[0],lonG.shape[-1]])
            DATA=np.squeeze(ncO[OGCMVar[i]][child_co,:,latO_co,lonO_co])
            tmp_mask_=np.invert(DATA.mask)
        
            for d in range(len(depthO)):
                # Interp mask
                tmp_mask=tmp_mask_[d]
                data_=griddata((lonO_s_m[tmp_mask].flatten(),latO_s_m[tmp_mask].flatten()),\
                              DATA[d][tmp_mask].flatten(),(lonO_s_m.flatten(),latO_s_m.flatten()),'nearest')
                data_=data_.reshape(latO_s_m.shape)
        
                # Interp 4 Grid
                data_re_=griddata( (lonO_s_m.flatten(),latO_s_m.flatten()), data_.flatten(), (lonG[:,:],latG[:,:]) ,'cubic' )
                data[d]=data_re_.reshape(lonG.shape) #.reshape(-1)
        OGCM_Data[i]=data

    def rho2u_3d(var):
        N,Mp,Lp=var.shape
        L=Lp-1
        var_u=0.5*(var[:,:,:L]+var[:,:,1:Lp])
        return var_u
    def rho2v_3d(var):
        T,Mp,Lp=var.shape
        M=Mp-1
        var_v=0.5*(var[:,:M,:]+var[:,1:Mp,:])
        return var_v
    def rho2u_2d(var):
        N,Lp=var.shape
        L=Lp-1
        var_u=0.5*(var[:,:L]+var[:,1:Lp])
        return var_u
    def rho2v_2d(var):
        Mp,Lp=var.shape
        M=Mp-1
        var_u=0.5*(var[:M,:]+var[1::Mp,:])
        return var_u
    
    cosa,sina=np.cos(angle)[:],np.sin(angle)[:] 
    
    # cosa=np.tile( np.tile(cosa,(thO,1,1)), (12,1,1,1) )
    # sina=np.tile( np.tile(sina,(thO,1,1)), (12,1,1,1) )
    cosa= np.tile(cosa,(thO,1,1))
    sina= np.tile(sina,(thO,1,1))
    
    #Process 2D vectors
    ubar= rho2u_2d(OGCM_Data['ubar'][:,:]*cosa[0,:,:]+OGCM_Data['vbar'][:,:]*sina[0,:,:]).squeeze()
    vbar= rho2v_2d(OGCM_Data['vbar'][:,:]*cosa[0,:,:]+OGCM_Data['ubar'][:,:]*sina[0,:,:]).squeeze()
    
    #Process 3D vectors
    u=rho2u_3d(OGCM_Data['u'][:,:,:]*cosa[:,:,:]+OGCM_Data['v'][:,:,:]*sina[:,:,:])
    v=rho2v_3d(OGCM_Data['v'][:,:,:]*cosa[:,:,:]-OGCM_Data['u'][:,:,:]*sina[:,:,:])
    
    # Process ROMS Vertical grid
    Z=np.zeros(len(depthO)+2)
    Z[0]=100;Z[1:-1]=-depthO;Z[-1]=-100000
    
    Rzeta=OGCM_Data['zeta'][:,:] # -1 for northern BRY
    
    zr=np.zeros([MyVar['Layer_N'],OGCM_Data['zeta'].shape[0],OGCM_Data['zeta'].shape[-1]])
    zw=np.zeros([MyVar['Layer_N']+1,OGCM_Data['zeta'].shape[0],OGCM_Data['zeta'].shape[-1]])
    
    zr[:,:,:]=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 1,topo,Rzeta);        
    zw[:,:,:]=jr.zlevs(MyVar['Vtransform'],MyVar['Vstretching'],MyVar['Theta_s'],\
             MyVar['Theta_b'],MyVar['Tcline'],MyVar['Layer_N'],\
                 5,topo,Rzeta);   
        
    zu=rho2u_3d(zr)[:,:,:]
    zv=rho2v_3d(zr)
    dzr=zw[1:,:,:]-zw[:-1,:,:] # [t,depth,lat,lon]
    dzu=rho2u_3d(dzr)
    dzv=rho2v_3d(dzr)
    
    # Add a level on top and bottom with no-gradient
    temp,salt=OGCM_Data['temp'][:,:,:],OGCM_Data['salt'][:,:,:]
    
    u1=np.hstack((np.expand_dims(u[0,:,:],axis=1)\
                  ,u,np.expand_dims(u[-1,:,:],axis=1)))
    v1=np.hstack((np.expand_dims(v[0,:,:],axis=1)\
                  ,v,np.expand_dims(v[-1,:,:],axis=1)))
    temp=np.hstack((np.expand_dims(temp[0,:,:],axis=1)\
                  ,temp,np.expand_dims(temp[-1,:,:],axis=1)))
    salt=np.hstack((np.expand_dims(salt[0,:,:],axis=1)\
                  ,salt,np.expand_dims(salt[-1,:,:],axis=1)))
    
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

    ncC=Dataset(My_Clm,mode='a')
    ncC['zeta'][child_co,:,:]=Rzeta
    # ncI['SSH'][:]=Rzeta
    ncC['temp'][child_co,:,:,:]=temp
    ncC['salt'][child_co,:,:,:]=salt
    ncC['u'][child_co,:,:,:]=u
    ncC['v'][child_co,:,:,:]=v
    ncC['ubar'][child_co,:,:]=ubar_
    ncC['vbar'][child_co,:,:]=vbar_
    ncC.close()
    

# TIMES_co_set=TIMES_co_set[25:]
Counter=0
Error_counter=0

# TIMES_co_set=[0]

stT=time.time()

while True: 

    try:
        with mp.Pool(12) as MyPool:
            results=MyPool.map(main,range(12))
        print('\n!!! Counter : '+str(Counter)+' !!!')
        Counter+=1;
        break

    except :
        Error_counter+=1
        print('\n!!! Error_count : '+str(Error_counter)+' !!!')
        
    if Error_counter==35:
        raise
print('\n!!! END Total Elapsed time : '+str((time.time()-stT)/60)[:4]+'min !!!')  





