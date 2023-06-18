# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:08:42 2023

@author: shjo9
"""

# Module import

import matplotlib
matplotlib.use('Agg') #Generates figures in Backend
import sys
sys.path.append('D:/OneDrive/JNUpack/JNUROMS/')
from Tools.dianose_module import *
from netCDF4 import Dataset, num2date, date2num
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import cmocean 
import numpy as np
import matplotlib.pyplot as plt


# ZETA ===================================================================
zeta_NN=15
zeta_lim=[-2,1.5]
zeta_levels=np.arange(zeta_lim[0],zeta_lim[-1]+0.2/2,0.1)
zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, zeta_NN,endpoint=True)))
# =============================================================================

# TMEP ===================================================================
temp_NN=15
temp_lim=[-2.5,28]
temp_levels=np.arange(temp_lim[0],temp_lim[-1]+2/2,2.)
temp_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')(np.linspace(0, 1, temp_NN,endpoint=True)))
# =============================================================================

# SALT ===================================================================
salt_NN=15
salt_lim=[32,37]
salt_levels=np.arange(salt_lim[0],salt_lim[-1]+.5/2,.5)
salt_CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, salt_NN,endpoint=True)))
# =============================================================================

# AICE ===================================================================
aice_NN=15
aice_lim=[0,1]
aice_levels=np.arange(aice_lim[0],aice_lim[-1]+.2/2,.1)
aice_CMAP = ListedColormap(cmocean.cm.ice(np.linspace(0, 1, aice_NN,endpoint=True)))
# =============================================================================


''' Parameters, CPP & Inputs '''

''' Topo, grid, z_layers '''
# Diag_TopoGrd()

Stability01.Log_pth='G:/TEST/LogTK25EM_DQD.txt'


''' Stability check 1,2 '''
Stability01()
Stability02() # Temp trand, aice trend

''' ZETA '''
data_drift('zeta',[-80,-24],zeta_levels,zeta_CMAP,zeta_lim,\
           mean='ann',st='1980-01',ed='1990-12')

data_drift('zeta',[-80,-24],zeta_levels,zeta_CMAP,zeta_lim,\
           mean='season',st='1980-01',ed='1990-12')

data_drift('temp',[-80,-24],temp_levels,temp_CMAP,temp_lim,\
           mean='monthly_clm',st='1980-01',ed='1990-12')

data_drift('salt',[-80,-24],salt_levels,salt_CMAP,salt_lim,\
           mean='monthly_clm',st='1980-01',ed='1990-12')

data_drift('aice',[-80,-24],aice_levels,aice_CMAP,aice_lim,\
           mean='monthly_clm',st='2000-01',ed='2001-12')

      
    
    
# zeta_AVISO(diff=True,{'1980-01','1990-12'})


''' Surface temperature drift '''
# SST_drift({'1980-01','1990-12'})
# SST_drift({'1980-01','1990-12'},'DJF')
# SST_drift({'1980-01','1990-12'},'JJA')

# SST_SODA(diff=True,{'1980-01','1990-12'})


''' Surface salinity drift '''
# SSS_drift({'1980-01','1990-12'})
# SSS_drift({'1980-01','1990-12'},'DJF')
# SSS_drift({'1980-01','1990-12'},'JJA')

''' Ice concetration '''
# SSS_drift({'1980-01','1990-12'})
# SSS_drift({'1980-01','1990-12'},'DJF')
# SSS_drift({'1980-01','1990-12'},'JJA')

''' Zonal mean of Subsurface salinity drift cross section '''
# zonal_data_drift('v_northward',zeta_CMAP,[-.01,.01],\
#                  mean='monthly',st='1980-02',ed='1990-12')

''' Zonal mean of Subsurface temperature cross section '''
# Z_SubT_drift({'1980-01','1990-12'})

'''Auger Temp vertical section '''
# Auger_temp_section('temp',temp_levels,temp_CMAP,\
#                    mean='monthly_clm',st='1980-01',ed='1990-12')

# =============================================================================
# A=Dataset(AVGS[0])['temp']

# B=np.mean(A,axis=2)

# DATA=xr.open_mfdataset(AVGS)['zeta','temp']

# tmp=DATA.mean(dim='xi_rho')[:5].values



# =============================================================================



''' Sv & Current Drake passage (66W) '''

''' Sv & Current South of Africa (24.5E) '''

''' Sv & Current South of Australia (129E) '''

''' Meridional overturing streamfunction (zonal mean)'''


''' WOA '''
# 













