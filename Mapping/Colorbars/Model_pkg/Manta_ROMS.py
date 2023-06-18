# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:13:49 2022

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean
import matplotlib as mpl
from copy import deepcopy   
import dask
import pandas as pd
from netCDF4 import Dataset,MFDataset,num2date,date2num
import datetime as dt

# Read Sample
    
class Log_Manager():
    def __init__(self,Log_dir,Log_name):
        with open(Log_dir+Log_name) as f:
            Model_Log=f.readlines()
        self.Model_Log = Model_Log
        
    def get_parameters(self,NN):
        Logs1000 = self.Model_Log[:1000]
        for i,j in zip(Logs1000,range(len(Logs1000))):
            if i.find('Physical Parameters')!=-1:
                stid=j
            elif i.find('Hout(idFsur)')!=-1:
                edid=j
                
            elif i.find('Aout(idFsur)')!=-1:
                edid=j
                 
        test2 = Logs1000[stid:edid]           
        test3 = [i.strip().split('    ')[0] for i in test2]
        test4 = [i.split(' '*2) for i in test3 if len(i.split(' '*2))==2]
        test5=pd.DataFrame(test4,columns=['Values','ID'])
        test5.index =test5['ID']
        test5.drop(columns='ID',inplace=True)
        Parameters01=np.concatenate([test5.index.values.reshape([len(test5),1]),test5.values],1)
        NN_=NN-len(Parameters01)*2%NN
        tmp_A=np.char.mod('%d', np.zeros([int(NN_/2),2])); tmp_A[tmp_A=='0']='-'
        Parameters01_re=np.concatenate([Parameters01,tmp_A],0)
        Parameters=Parameters01_re.reshape([int(len(Parameters01_re.reshape(-1))/NN),NN])
        return Parameters
        
def get_model_stability(self):
        i=0;j=0; stid=10**5; 
        Model_Log = self.Model_Log
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
            if i>stid and len(tmp1) and tmp1[0].isnumeric():
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
    

class Sample_Manager():
    def __init__(self,r_path,input_name='ocean_avg_0001.nc'):
        Sample = Dataset(r_path+input_name,'r')
        nc_attrs = Sample.ncattrs()
        Model_info = []
        for i in nc_attrs:
            Model_info.append(Sample.getncattr(i).replace(',','\n\t\t'))
        Model_info = pd.DataFrame(Model_info,index=nc_attrs)
        self.Model_info = Model_info

    def get_cpp_Sample(self,Mode_N=7):
        # CPP
        CPP_ = self.Model_info.loc['CPP_options']
        CPP_ = CPP_[0].split('\n\t\t')[1:]
        tv = int(len(CPP_)/Mode_N)
        tmp = np.zeros([tv+1,Mode_N])
        empty_list = [list(i) for i in tmp]
        i=j=0
        for k in range(len(CPP_)):
            empty_list[i][j]=CPP_[k].strip()
            j+=1
            if j==Mode_N:
                j=0;i+=1
        CPP = pd.DataFrame(empty_list).replace(0.0,'-')
        return CPP

    def get_LBC_Sample(self):
        # NLM_LBC
        NLM_LBC_ = self.Model_info.loc['NLM_LBC'][0]
        NLM_LBC_ = NLM_LBC_.split('\n')
        LBC_Col_ = NLM_LBC_[1].split(':')[1].rstrip()
        LBC_Col = [i for i in LBC_Col_.split(' ') if len(i)]
    
        LBC_index=[]; LBC=[]  
        for i in NLM_LBC_[1:]:        
            tmp=i.split(':')
            LBC_index.append(tmp[0])
            if i==NLM_LBC_[1:][0]:
                LBC.append(tmp[1].strip().split(' '*2))
            else:
                LBC.append(tmp[1].strip().split(' '*4))
        
        NLM_LBC=pd.DataFrame(LBC,columns=LBC_Col,index=LBC_index)
        NLM_LBC['index'] = LBC_index
        NLM_LBC=NLM_LBC.reindex(columns=['index','WEST','EAST','SOUTH','NORTH'])
        NLM_LBC.replace('EDGE','',inplace=True)
        # NLM_LBC.concat(LBC_Col,axis=1)
        # NLM_LBC.transpose()
        # NLM_LBC_Col = NLM_LBC.columns.values
        return NLM_LBC
    
    def get_inputs_files(self):
        Files=self.Model_info.loc[[i for i in self.Model_info.index if i.find('file')!=-1]]
        Files=np.concatenate([Files.index.values.reshape(len(Files),1),Files.values],1)
        return Files


    



    
    
class Fig_manager():
    def __init__(self):
        # =============================================================================
        # Figure configuration
        # =============================================================================
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



# PD.columns

# plt.plot(PD['Day']/365-PD['Day'][0]/365,PD['TOTAL_ENRG'])

# # =============================================================================
# # Figure 1. Check model stability   
# # =============================================================================
# Label_size = 25
# fig, axs = plt.subplots(3,1,figsize=(11,6.5),constrained_layout = True,
#                         sharex=True,gridspec_kw={'height_ratios': [1, 1.,1]},dpi=200)
# f1 = axs[0].plot(Model_Times1,Zeta2d.data, label='Zeta',color='k',linewidth=2,zorder=0)
# axs[0].plot(np.array(Model_Times1)[t_co],f1_z2(tmp_t),color='r',linewidth=3,linestyle='dashed')
# axs[0].tick_params(axis='y', labelsize=Label_size)
# axs[0].set_xlim('1993-01','2015-04')
# xtick_location = Model_Times1[0::12*2]
# xtick_labels = Model_Times2[0::12*2]
# axs[0].set_xticks(ticks=xtick_location)
# axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
# # axs[0].grid(axis='x', alpha=.3,linestyle='-.',color='k')
# # axs[0].set_ylim(Zeta2d.data.mean()-Zeta2d.data.std()*2.5,\
# #                 Zeta2d.data.mean()+Zeta2d.data.std()*7)# axs[0].set_yticks(ticks=np.arange(18,23,1))
# axs[0].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
# axs[0].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
# #! Fig2 
# f1 = axs[1].plot(Model_Times1,Momentum2d.data, label='Momentum2d',color='k',linewidth=2,zorder=0)
# axs[1].plot(np.array(Model_Times1)[t_co],f1_m2(tmp_t),color='r',linewidth=3,linestyle='dashed')
# axs[1].tick_params(axis='y', labelsize=Label_size)
# xtick_location = Model_Times1[5::12*2]
# xtick_labels =Model_Times2[5::12*2]
# axs[1].set_xticks(ticks=xtick_location)
# axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[1].set_xlim('1993-01','2015-04')
# # axs[1].set_ylim(Momentum2d.data.mean()-Momentum2d.data.std()*2.5,\
# #                 Momentum2d.data.mean()+Momentum2d.data.std()*3.)
# # axs[1].set_yticks(ticks=np.arange(18,23,1))
# axs[1].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
# f1 = axs[2].plot(Model_Times1,Momentum3d.data, label='Momentum3d',color='k',linewidth=2,zorder=0)
# axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
# axs[2].tick_params(axis='y', labelsize=Label_size)
# xtick_location = Model_Times1[5::12*2]
# xtick_labels =Model_Times2[5::12*2]
# axs[2].set_xticks(ticks=xtick_location)
# axs[2].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
# axs[2].set_xlim('1993-01','2015-04')
# # axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
# #                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
# # axs[1].set_yticks(ticks=np.arange(18,23,1))
# axs[2].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
# plt.tight_layout()
# if fig_bool:
#     # plt.savefig(w_path_sig+'ppt/03Corr_indices',
#     #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
#     plt.savefig(w_dir+'Model_momentum')
# plt.show()














