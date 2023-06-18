# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:15:39 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

import xarray as xr
from netCDF4 import Dataset,MFDataset
import os
import numpy as np
import matplotlib.pyplot as plt

data_path='G:/Models/TK0525ED_CLM/'
save_dir1='/data4/base158/Warehouse01/MyPro01/CaseTK_11CHA/'

fig_bool=0

ncs=[data_path+i for i in os.listdir(data_path)]

SampleNC=Dataset(ncs[0])

Variables=SampleNC.variables.keys()

PN,PM=SampleNC['pn'][:]*10**3,SampleNC['pm'][:]*10**3 # km
CELL_size=1/PN*1/PM # m**2


# ZETA
ZETA=MFDataset(ncs)['zeta']
ZETA_values=ZETA[:]
ZETA_area=np.zeros_like(ZETA_values)
n=0
for i in ZETA_values:
    ZETA_area[n]=i*CELL_size
    n+=1

Total_size=np.sum(CELL_size)

ZETA_timeseries=np.sum(ZETA_area,axis=2).sum(axis=1)/Total_size

plt.plot(ZETA_timeseries[-500:])

# AICE
AICE=MFDataset(ncs)['aice']
AICE_values=AICE[:]
AICE_area=np.zeros_like(AICE_values)
n=0
for i in AICE_values:
    AICE_area[n]=i*CELL_size
    n+=1

AICE_timeseries=np.sum(AICE_area,axis=2).sum(axis=1)/Total_size


t=np.arange(len(AICE_timeseries))

AICE_fp1=np.polyfit(t[200:],AICE_timeseries[200:],1)
AICE_trend=np.polyval(AICE_fp1,t[200:])

plt.plot(AICE_timeseries[200:])
plt.plot(AICE_trend)


# SHFLUX 
SHFLUX=MFDataset(ncs)['shflux']
SHFLUX_values=SHFLUX[:]
SHFLUX_area=np.zeros_like(SHFLUX_values)
n=0
for i in SHFLUX_values:
    SHFLUX_area[n]=i*CELL_size
    n+=1

SHFLUX_timeseries=np.sum(SHFLUX_area,axis=2).sum(axis=1)/Total_size


t=range(len(SHFLUX_timeseries))

SHFLUX_fp1=np.polyfit(t,SHFLUX_timeseries,1)
SHFLUX_trend=np.polyval(SHFLUX_fp1,t)

plt.plot(SHFLUX_timeseries)
plt.plot(SHFLUX_trend)

# SST
SST=MFDataset(ncs)['temp'][:,-1,:,:]
SST_values=SST[:]
SST_area=np.zeros_like(SST_values)
n=0
for i in SST_values:
    SST_area[n]=i*CELL_size
    n+=1

SST_timeseries=np.sum(SST_area,axis=2).sum(axis=1)/Total_size


t=range(len(SST_timeseries))

SST_fp1=np.polyfit(t,SST_timeseries,1)
SST_trend=np.polyval(SST_fp1,t)

plt.plot(SST_timeseries)
plt.plot(SST_trend)


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True

plt.rcParams["font.family"] = 'Helvetica'
# mpl.rcParams['axes.unicode_minus'] = False


Model_Times1 = t
Model_Times2=[str(i/12)[:2] for i in t]


Label_size = 25
fig, axs = plt.subplots(4,1,figsize=(11,9.),constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': [1, 1.,1,1]},dpi=200)
f1 = axs[0].plot(Model_Times1,ZETA_timeseries, label='ZETA',color='k',linewidth=2,zorder=0)
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
axs[0].legend(fontsize=18)
axs[0].set_ylim(-0.4,0.3)
#! Fig2 
f1 = axs[1].plot(Model_Times1,AICE_timeseries, label='AICE',color='k',linewidth=2,zorder=0)
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
axs[1].legend(fontsize=18,loc=4)

f1 = axs[2].plot(Model_Times1,SHFLUX_timeseries, label='SHFLX',color='k',linewidth=2,zorder=0)
# axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
axs[2].tick_params(axis='y', labelsize=Label_size)
xtick_location = Model_Times1[0::12*10]
xtick_labels =Model_Times2[0::12*10]
axs[2].set_xticks(ticks=xtick_location)
axs[2].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
# axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
#                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
# axs[1].set_yticks(ticks=np.arange(18,23,1))
axs[2].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs[2].legend(fontsize=18,loc=4)
axs[2].set_ylim(-100,100)


f1 = axs[3].plot(Model_Times1,SST_timeseries, label='SST',color='k',linewidth=2,zorder=0)
# axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
axs[3].tick_params(axis='y', labelsize=Label_size)
xtick_location = Model_Times1[0::12*10]
xtick_labels =Model_Times2[0::12*10]
axs[3].set_xticks(ticks=xtick_location)
axs[3].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
# axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
#                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
# axs[1].set_yticks(ticks=np.arange(18,23,1))
axs[3].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs[3].legend(fontsize=18,loc=4)
axs[3].set_ylim(7,12.6)

plt.tight_layout()
if fig_bool:
    plt.savefig(save_dir1+'ppt/'+'Model_stability2',
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(save_dir1+'Model_stability2')
plt.show()



