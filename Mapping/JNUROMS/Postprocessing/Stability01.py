# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:15:09 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

PKG_path = 'D:/OneDrive/base142/Factory/Model_results/Modules/Model_pkg'
import sys 
sys.path.append(PKG_path)
# from Manta_ROMS import Log_Manager, Sample_Manager
import numpy as np
# from pptx import Presentation # ????? 
# from pptx.util import Inches,Cm, Pt # ??, ??? ??? ??
# from pptx.enum.text import PP_ALIGN
from netCDF4 import Dataset,MFDataset
import os
import xarray as xr
import dask
import matplotlib as mpl
import cmocean
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from scipy.interpolate import interp2d, griddata

data_path='G:/Models/'
Log_name='test_log.txt'

tmp_name=''

fig_bool=0
save_dir1='/data4/base158/Warehouse01/CaseSODA/TK0525EM/'

with open(data_path+Log_name) as f:
    Model_Log=f.readlines()
    Logs1000 = Model_Log[:1000]


def get_model_stability(Model_Log):
        i=0;j=0; stid=10**5; 
        while i<=len(Model_Log):
            if i>=len(Model_Log):
                break
            elif Model_Log[i].find('rx0') !=-1:
                rx0=Model_Log[i].split(':')[-1].strip().replace(' ','')
                rx1=Model_Log[i+1].split(':')[-1].strip().replace(' ','')
                rx=rx0+'\n'+rx1
            if Model_Log[i].find('STEP')!=-1:
                header_name = [ii.strip() for ii in Model_Log[i].split(' ') if len(ii)]
                [header_name.append(ii.strip()) for ii in Model_Log[i+1].split('  ') if len(ii)]
                PD = pd.DataFrame({},columns=header_name)
                stid=i
            tmp1=Model_Log[i].lstrip()
            if i>stid and len(tmp1) and tmp1[0].isnumeric() and not i>len(Model_Log)-300:
                A= [ii.strip() for ii in tmp1.split(' ') if len(ii)]
                [A.append(ii.strip()) for ii in Model_Log[i+1].lstrip().split(' ') if len(ii)]
                PD.loc[j] = A
                i+=2; j+=1
            i+=1
        for i in PD.columns:
            try:
                PD[i] = PD[i].astype(float)
            except :
                pass
        return PD,rx

PD,rx=get_model_stability(Model_Log)


PD['NET_VOLUME'].plot()
PD.columns


RX=rx.split('\n')
RX1=RX[0][:16]
RX2=RX[-1][:-7]

Title_name='Topo ('+RX1+' / '+RX2+')'



plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True

plt.rcParams["font.family"] = 'Helvetica'
mpl.rcParams['axes.unicode_minus'] = False


Model_Times1 = PD['YYYY-MM-DD']
Model_Times2 = [i[2:4] for i in Model_Times1]


Label_size = 25
fig, axs = plt.subplots(4,1,figsize=(11,8.5),constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': [1,1, 1.,1]},dpi=200)
f1 = axs[0].plot(Model_Times1,PD['KINETIC_ENRG'], label='KINETIC_ENRG',color='k',linewidth=2,zorder=0)
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
#! Fig2 
f1 = axs[1].plot(Model_Times1,PD['POTEN_ENRG'], label='POTEN_ENRG',color='k',linewidth=2,zorder=0)
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

f1 = axs[2].plot(Model_Times1,PD['TOTAL_ENRG'], label='TOTAL_ENRG',color='k',linewidth=2,zorder=0)
# axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
axs[2].tick_params(axis='y', labelsize=Label_size)
xtick_location = Model_Times1[5::12*40]
xtick_labels =Model_Times2[5::12*40]
axs[2].set_xticks(ticks=xtick_location)
axs[2].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
# axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
#                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
# axs[1].set_yticks(ticks=np.arange(18,23,1))
axs[2].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs[2].legend(fontsize=18,loc=4)

f1 = axs[3].plot(Model_Times1,PD['NET_VOLUME'], label='NET_VOLUME',color='k',linewidth=2,zorder=0)
# axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
axs[3].tick_params(axis='y', labelsize=Label_size)
xtick_location = Model_Times1[5::12*60]
xtick_labels =Model_Times2[5::12*60]
axs[3].set_xticks(ticks=xtick_location)
axs[3].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
# axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
#                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
# axs[1].set_yticks(ticks=np.arange(18,23,1))
axs[3].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs[3].legend(fontsize=18,loc=4)



plt.tight_layout()
if fig_bool:
    plt.savefig(save_dir1+'ppt/'+'Model_momentum_logs'+tmp_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(save_dir1+'Model_momentum_logs'+tmp_name)
plt.show()













