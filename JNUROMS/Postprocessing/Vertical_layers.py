# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 18:05:27 2023

@author: shjo9
"""


PKG_path = 'D:/JNUpack/JNUROMS'
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
from matplotlib.colors import ListedColormap,LinearSegmentedColormap


Grd_name='D:/OneDrive/base142/Warehouse01/Grd_SO_05d.nc'
Data_name='D:/OneDrive/base142/Factory/MantaROMS/Test_nc/Ini_soda_05d_jhlee_198002.nc'

fig_bool=0

w_path_sig='D:/Working_hub/OneDrive/base142/Warehouse01/MyPro01/TK0525EM/'

ncG,ncD=Dataset(Grd_name),Dataset(Data_name)

LON,LAT=ncG['lon_rho'][:],ncG['lat_rho'][:]

Coord1=np.where((LAT[:,0]>=-60.1)&(LAT[:,0]<=-59.9))[0]

zr1=jr.zlevs(ncD['Vtransform'][:], 4,10, 4,\
            ncD['hc'][:], ncD['sc_r'][:].shape[0],1, ncG['h'][Coord1,:], ncD['zeta'][0,Coord1,:])

Coord2=np.where( ncG['h'][:]>5999)

zr2=jr.zlevs(ncD['Vtransform'][:], ncD['Vstretching'][:],ncD['theta_s'][:], ncD['theta_b'][:],\
            ncD['hc'][:], ncD['sc_r'][:].shape[0],1, ncG['h'][Coord2[0][:],Coord2[1][:]], ncD['zeta'][0,Coord2[0][:],Coord2[1][:]])

# zr[zr>1000]=np.nan

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams["font.family"] = 'Helvetica'
cmap_list=['#333333','#7F00BF','#9200E4','#AD07FF','#C23EFF','#DA87FF','#BC0D0F','#B30000','#C91F25',\
   '#D9343E','#E84A56','#F65E67','#FF6E6E','#FF8051','#FF8B1B','#FF9C00','#FFAA09','#FFBC1E',\
       '#FFD039','#FFE256','#FFF26F','#008000','#00A400','#00D500','#1EF31E','#69FC69','#0077B3',\
           '#008DDE','#07ABFF','#3EC1FF','#000000','#000000']
MyCmap=ListedColormap(cmap_list).reversed()

Label_size = 14

# Figure 1
fig, axs = plt.subplots(figsize=(4.5,5),constrained_layout = True,dpi=200)
axs.plot(np.diff(-zr2[:,0,0]), label='pc01',color='C0',linewidth=1.5,zorder=0)
plt.scatter(range(49),np.diff(-zr2[:,0,0]),c=np.arange(49), label='pc01',\
            cmap=MyCmap)# axs.set_yticks(ticks=np.arange(18,23,1))
axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
# axs.set_ylim(-2.5,3.5)
axs.set_xlim(-1,50)
# xtick_location = KVTe.index.values.tolist()[5::12*2]
# xtick_labels = XL.values.tolist()[5::12*2]
# axs.set_xticks(ticks=xtick_loca1, rotation=0, fontsize=Label_size, alpha=0)
axs.tick_params(axis='both', labelsize=Label_size)
axs.spines["top"].set_alpha(1.)
axs.spines["bottom"].set_alpha(1.)
axs.spines["right"].set_alpha(1.0)
axs.spines["left"].set_alpha(1.)
plt.tight_layout()
if fig_bool:
    plt.savefig(w_path_sig+'ppt/Model_vertical_layers1',
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path_sig+'Model_vertical_layers1')
plt.show()

# Figure 2
fig, axs = plt.subplots(figsize=(4.5,5),constrained_layout = True,dpi=200)
axs.plot(zr2[:,0,0], label='pc01',color='C0',linewidth=1.5,zorder=0)
plt.scatter(range(50),zr2[:,0,0],c=np.arange(50), label='pc01',\
            cmap=MyCmap)# axs.set_yticks(ticks=np.arange(18,23,1))
axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
# axs.set_ylim(-2.5,3.5)
axs.set_xlim(-1,50)
# xtick_location = KVTe.index.values.tolist()[5::12*2]
# xtick_labels = XL.values.tolist()[5::12*2]
# axs.set_xticks(ticks=xtick_loca1, rotation=0, fontsize=Label_size, alpha=0)
axs.tick_params(axis='both', labelsize=Label_size)
axs.spines["top"].set_alpha(1.)
axs.spines["bottom"].set_alpha(1.)
axs.spines["right"].set_alpha(1.0)
axs.spines["left"].set_alpha(1.)
plt.tight_layout()
if fig_bool:
    plt.savefig(w_path_sig+'ppt/Model_vertical_layers2',
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path_sig+'Model_vertical_layers2')
plt.show()


Label_size=9
# Figure 3
fig, axs = plt.subplots(2,1,figsize=(6,4),constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': [1, 1]},dpi=200)
for i in zr1:
    axs[0].plot(LON[0,:],i, label='pc01',color='k',linewidth=.5,zorder=0)
axs[0].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
axs[0].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size, color='k')
axs[0].set_ylim(-ncD['hc'][:],0)
axs[0].set_xlim(0,LON[0,-1])

for i in zr1:
    axs[1].plot(LON[0,:],i, label='pc01',color='k',linewidth=.5,zorder=0)
axs[1].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
axs[1].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size,  color='k')
axs[1].set_ylim(-6000,-ncD['hc'][:])
axs[1].set_xlim(0,LON[0,-1])

plt.tight_layout(pad=.1)
if fig_bool:
    plt.savefig(w_path_sig+'ppt/Model_vertical_layers3',
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path_sig+'Model_vertical_layers3')
plt.show()





    
    
    
    