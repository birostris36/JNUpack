# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:30:36 2023

@author: shjo9
"""

import xarray as xr
from netCDF4 import Dataset
import os 
import numpy as np
import matplotlib.pyplot as plt

Ross_lat_rng,Ross_lon_rng=[-73,-60], [180,215]
Wedd_lat_rng,Wedd_lon_rng=[-73,-60], [300,360]

# Avg_pth='G:/Models/Cases_CLM/'
# save_pth='G:/Models/Cases_CLM/tmp/'

# =============================================================================
#Avg_pth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_XSICE/Outputs_srfd/Avg/'
#Log_npth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_XSICE/LOGS/Log_V205701_M_XSICE_sr     fd.txt'
#save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/Control/'

#Avg_pth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_XSICE/Outputs_srfd_xclm/Avg/'
#Log_npth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_XSICE/LOGS/Log_V205701_M_XSICE_sr     fd_xclm.txt'
#save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/Control_xclm/'

#Avg_pth='/home/jejunu/JNU_shjo/_data/base201/ROMSLAB/room_ks01/Outputs/Avg/'
#Log_npth='/home/jejunu/JNU_shjo/_data/base201/ROMSLAB/room_ks01/LOGS/C_Log_V205701_XSICE.txt'
#save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/Kate_svn/'

#Avg_pth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_LSXICE_2/Outputs/Avg/'
#Log_npth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_LSXICE_2/LOGS/C_Log_V205701_XSICE     _02.txt'
#save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/ROMS1053_LSXICE/'

Avg_pth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_LSXICE/Outputs/Avg/'
Log_npth='/home/jejunu/JNU_shjo/_data/base201/Cases_A/V205701_LSXICE/LOGS/Log_V205701_LSXICE.txt     '
save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/V205701_LSXICE/'

#Avg_pth='/home/jejunu/JNU_shjo/_data/base201/Cases_CLM/TEST/Outputs/Avg/'
#Log_npth='/home/jejunu/JNU_shjo/_data/base201/Cases_CLM/TEST/LOGS/Log_TEST_01.txt'
#save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/Control_CLM/'

#Avg_pth='/home/jejunu/JNU_shjo/_data/base201/Cases_CLM/Q0701_CLM/Outputs_CLM/Avg01/'
#Log_npth='/home/jejunu/JNU_shjo/_data/base201/Cases_CLM/Q0701_CLM/LOGS/LOG_Q0701_CLM_02.txt'
#save_pth='/home/jejunu/JNU_shjo/_data/base201/Warehouse02/Q0701_CLM/'

# =============================================================================


fig_bool=1

AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]

Sample_Data=Dataset(AVGS[0])
lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
Ross_lat_co=np.where((lat_rho[:,0]>=Ross_lat_rng[0])&(lat_rho[:,0]<=Ross_lat_rng[-1]))[0]
Ross_lon_co=np.where((lon_rho[0,:]>=Ross_lon_rng[0])&(lon_rho[0,:]<=Ross_lon_rng[-1]))[0]

Wedd_lat_co=np.where((lat_rho[:,0]>=Wedd_lat_rng[0])&(lat_rho[:,0]<=Wedd_lat_rng[-1]))[0]
Wedd_lon_co=np.where((lon_rho[0,:]>=Wedd_lon_rng[0])&(lon_rho[0,:]<=Wedd_lon_rng[-1]))[0]


ROSS=xr.open_mfdataset(Avg_pth+'*.nc').zeta.loc[dict(xi_rho=Ross_lon_co,eta_rho=Ross_lat_co)].mean(dim=['eta_rho','xi_rho']).values
WEDD=xr.open_mfdataset(Avg_pth+'*.nc').zeta.loc[dict(xi_rho=Wedd_lon_co,eta_rho=Wedd_lat_co)].mean(dim=['eta_rho','xi_rho']).values


Model_Times1 = np.arange(len(ROSS))
Model_Times2=[str(i/12)[:2] for i in np.arange(len(ROSS))]

Label_size = 25
fig, axs = plt.subplots(2,1,figsize=(11,4.5),constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': [1, 1.]},dpi=200)
f1 = axs[0].plot(Model_Times1,ROSS, label='ROSS zeta',color='k',linewidth=2,zorder=0)
# axs[0].plot(np.array(Model_Times1)[t_co],f1_z2(tmp_t),color='r',linewidth=3,linestyle='dashed')
axs[0].tick_params(axis='y', labelsize=Label_size)
# axs[0].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
xtick_location = Model_Times1[5::12]
xtick_labels = Model_Times2[5::12]
axs[0].set_xticks(ticks=xtick_location)
axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
# axs[0].grid(axis='x', alpha=.3,linestyle='-.',color='k')
# axs[0].set_ylim(Zeta2d.data.mean()-Zeta2d.data.std()*2.5,\
#                 Zeta2d.data.mean()+Zeta2d.data.std()*7)# axs[0].set_yticks(ticks=np.arange(18,23,1))
axs[0].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs[0].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
axs[0].legend(fontsize=18,loc='upper left')
# axs[0].set_ylim(-0.4,0.3)
#! Fig2 
f1 = axs[1].plot(Model_Times1,WEDD, label='WEDD zeta',color='k',linewidth=2,zorder=0)
# axs[1].plot(np.array(Model_Times1)[t_co],f1_m2(tmp_t),color='r',linewidth=3,linestyle='dashed')
axs[1].tick_params(axis='y', labelsize=Label_size)
xtick_location = Model_Times1[5::12*10]
xtick_labels =Model_Times2[5::12*10]
axs[1].set_xticks(ticks=xtick_location)
axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[1].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
# axs[1].set_ylim(Momentum2d.data.mean()-Momentum2d.data.std()*2.5,\
#                 Momentum2d.data.mean()+Momentum2d.data.std()*3.)
# axs[1].set_yticks(ticks=np.arange(18,23,1))
axs[1].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs[1].legend(fontsize=18,loc='upper left')
plt.tight_layout()
if fig_bool:
    plt.savefig(save_pth+'ppt/'+'ROSS_WEDD_zeta_time',
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(save_pth+'ROSS_WEDD_zeta_time')
plt.show()