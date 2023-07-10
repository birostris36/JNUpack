# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:08:42 2023

@author: shjo9
"""

# # Module import

# import matplotlib
# matplotlib.use('Agg') #Generates figures in Backend
import sys
sys.path.append('D:/JNUpack/JNUROMS/')
from Tools.dianose_module4CLM import *
from netCDF4 import Dataset, num2date, date2num
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import cmocean 
import numpy as np
import matplotlib.pyplot as plt

# ZETA ===================================================================
zeta_NN=16
zeta_lim=[-2,1.5]
zeta_levels=np.arange(zeta_lim[0],zeta_lim[-1]+0.2/2,0.2)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 0.9, len(zeta_levels)+1,endpoint=True)))

zeta_diff_lim=[-.8,.8]
zeta_diff_levels=np.arange(zeta_diff_lim[0],zeta_diff_lim[-1]+0.1/2,0.1)
zeta_diff_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, len(zeta_diff_levels)+1,endpoint=True)))


# =============================================================================

# TMEP ===================================================================
temp_NN=15
temp_lim=[-2.5,28]
temp_levels=np.arange(temp_lim[0],temp_lim[-1]+2/2,2.)
temp_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')(np.linspace(0, 1, temp_NN,endpoint=True)))

temp_diff_lim=[-3,3]
temp_diff_levels=np.arange(temp_diff_lim[0],temp_diff_lim[-1]+.5/2,.5)

temp_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')(np.linspace(0, 1, temp_NN,endpoint=True)))

# =============================================================================

# SALT ===================================================================
salt_NN=15
salt_lim=[33,36.5]
salt_levels=np.arange(salt_lim[0],salt_lim[-1]+.2/2,.2)
salt_CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, len(salt_levels),endpoint=True)))

salt_diff_lim=[-.5,.5]
salt_diff_levels=np.arange(salt_diff_lim[0],salt_diff_lim[-1]+.05/2,.05)
salt_diff_CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, len(salt_diff_levels),endpoint=True)))

# =============================================================================

# AICE ===================================================================
aice_NN=15
aice_lim=[0,1]
aice_levels=np.arange(aice_lim[0],aice_lim[-1]+.2/2,.1)
aice_CMAP = ListedColormap(cmocean.cm.ice(np.linspace(0, 1, aice_NN,endpoint=True)))
# =============================================================================

u_NN=16
u_lim=[-.8,.8]
u_levels=np.arange(u_lim[0],u_lim[-1]+0.1/2,0.1)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
uv_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, len(u_levels)+1,endpoint=True)))

u_diff_lim=[-.5,.5]
u_diff_levels=np.arange(u_diff_lim[0],u_diff_lim[-1]+0.1/2,0.1)
uv_diff_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, len(u_diff_levels)+1,endpoint=True)))


CMAP_LT=CMAP=plt.get_cmap('seismic',15)

# data_drift('u_eastward',[-80,-24],u_levels,uv_CMAP,u_lim,\
#             mean='ann',st='2016-01',ed='2016-12')
# Surface_data_Soda_diff('u_eastward',[-80,-24],uv_CMAP,uv_diff_CMAP,u_levels,\
#                         u_diff_levels,u_lim,u_diff_lim,mean='monthly',st='2016-03',ed='2016-12')

# ''' My ppt '''
# check_inputs()
lon_rng= 295
lat_rng=[-65,-55]
# ''' Stability check 1,2 '''
# Stability01()
# Stability03(lon_rng,lat_rng) # Temp trand, aice trend

# ''' Horizontal mean of Surface layer data '''
# data_drift('zeta',[-80,-24],zeta_levels,zeta_CMAP,zeta_lim,\
#             mean='ann',st='0001-01',ed='0001-12')
# data_drift('zeta',[-80,-24],zeta_levels,zeta_CMAP,zeta_lim,\
#             mean='monthly_clm',st='2000-01',ed='2016-12')
    
# data_drift('u_eastward',[-80,-24],u_levels,uv_CMAP,u_lim,\
#             mean='ann',st='1980-02',ed='2016-12')
# data_drift('u_eastward',[-80,-24],u_levels,uv_CMAP,u_lim,\
#             mean='monthly_clm',st='2000-01',ed='2016-12')

# data_drift('temp',[-80,-24],temp_levels,temp_CMAP,temp_lim,\
#             mean='ann',st='1980-02',ed='2016-12')
# data_drift('temp',[-80,-24],temp_levels,temp_CMAP,temp_lim,\
#             mean='monthly_clm',st='2000-01',ed='2016-12')

# data_drift('salt',[-80,-24],salt_levels,salt_CMAP,salt_lim,\
#             mean='ann',st='1980-02',ed='2016-12')
# data_drift('salt',[-80,-24],salt_levels,salt_CMAP,salt_lim,\
#             mean='monthly_clm',st='2000-01',ed='2016-12')
    
# data_drift('aice',[-80,-24],aice_levels,aice_CMAP,aice_lim,\
#             mean='ann',st='1980-01',ed='2016-12')
# data_drift('aice',[-80,-24],aice_levels,aice_CMAP,aice_lim,\
#             mean='monthly_clm',st='2000-01',ed='2016-12')

# Surface_data_Soda_diff('zeta',[-80,-24],zeta_CMAP,zeta_diff_CMAP,zeta_levels,\
#                         zeta_diff_levels,zeta_lim,zeta_diff_lim,mean='ann',st='1980-02',ed='2016-12')
# Surface_data_Soda_diff('zeta',[-80,-24],zeta_CMAP,zeta_diff_CMAP,zeta_levels,\
#                         zeta_diff_levels,zeta_lim,zeta_diff_lim,mean='monthly_clm'\
#                             ,st='2000-01',ed='2016-12')
    
# Surface_data_Soda_diff('temp',[-80,-24],temp_CMAP,temp_CMAP,temp_levels,\
#                         temp_diff_levels,temp_lim,temp_diff_lim,mean='ann',st='1980-02',ed='2016-12')
# Surface_data_Soda_diff('temp',[-80,-24],temp_CMAP,temp_CMAP,temp_levels,\
#                         temp_diff_levels,temp_lim,temp_diff_lim,mean='monthly_clm'\
#                             ,st='2000-01',ed='2016-12')
    
# Surface_data_Soda_diff('salt',[-80,-24],salt_CMAP,salt_diff_CMAP,salt_levels,\
#                         salt_diff_levels,salt_lim,salt_diff_lim,mean='ann',st='1980-02',ed='2016-12')
# Surface_data_Soda_diff('salt',[-80,-24],salt_CMAP,salt_diff_CMAP,salt_levels,\
#                         salt_diff_levels,salt_lim,salt_diff_lim,mean='monthly_clm'\
#                             ,st='2000-01',ed='2016-12')
    
# Surface_data_Soda_diff('u_eastward',[-80,-24],uv_CMAP,uv_diff_CMAP,u_levels,\
#                         u_diff_levels,u_lim,u_diff_lim,mean='ann',st='1980-02',ed='2016-12')
# Surface_data_Soda_diff('u_eastward',[-80,-24],uv_CMAP,uv_diff_CMAP,u_levels,\
#                         u_diff_levels,u_lim,u_diff_lim,mean='monthly_clm'\
#                             ,st='2000-01',ed='2016-12')

# ''' Zonal mean of Subsurface salinity drift cross section '''
# zonal_data_drift('u_eastward',uv_CMAP,[-.1,.1],\
#                   mean='ann',st='1980-02',ed='2016-12')
# zonal_data_drift('u_eastward',uv_CMAP,[-.1,.1],mean='monthly_clm',\
#                   st='2000-01',ed='2016-02')   
    
# zonal_data_drift('v_northward',uv_CMAP,[-.01,.01],\
#                   mean='ann',st='1980-02',ed='2016-12')
# zonal_data_drift('v_northward',uv_CMAP,[-.01,.01],mean='monthly_clm',\
#                   st='2000-01',ed='2016-02')   
# zonal_data_drift('temp',temp_CMAP,temp_lim,mean='season',\
#                  st='0001-01',ed='0001-12')
# zonal_data_drift('salt',salt_CMAP,salt_lim,mean='season',\
#                  st='2000-01',ed='2016-12')   
# zonal_data_diff_Soda('temp',temp_CMAP,temp_CMAP,temp_lim,[-2,2],\
#                      mean='season',st='2000-01',ed='2016-12')
# zonal_data_diff_Soda('salt',salt_CMAP,salt_CMAP,salt_lim,[-.3,.3],\
#                      mean='season',st='2000-01',ed='2016-12')   
    
# zonal_data_drift('temp',temp_CMAP,temp_lim,\
#                   mean='ann',st='1980-02',ed='2016-12')
# zonal_data_drift('temp',temp_CMAP,temp_lim,mean='monthly_clm',\
#                  st='2000-01',ed='2016-12')

# zonal_data_drift('salt',salt_CMAP,salt_lim,\
#                   mean='ann',st='1980-02',ed='2016-12')
# zonal_data_drift('salt',salt_CMAP,salt_lim,mean='monthly_clm',\
#                  st='2000-01',ed='2016-12')   
    
# ''' Zonal mean of Subsurface cross section '''
# zonal_data_diff_Soda('temp',temp_CMAP,temp_CMAP,temp_lim,[-2,2],\
#                      mean='monthly_clm',st='2000-01',ed='2016-12')
# zonal_data_diff_Soda('temp',temp_CMAP,temp_CMAP,temp_lim,[-2,2],\
#                      mean='ann',st='1980-02',ed='2016-12')
    
# zonal_data_diff_Soda('salt',salt_CMAP,salt_CMAP,salt_lim,[-.3,.3],\
#                      mean='monthly_clm',st='2000-01',ed='2016-12')
# zonal_data_diff_Soda('salt',salt_CMAP,salt_CMAP,salt_lim,[-.3,.3],\
#                      mean='ann',st='1980-02',ed='2016-12')
    
    
# zonal_data_diff_Soda('u_eastward',uv_CMAP,uv_CMAP,[-.08,.08],[-.05,.05],\
#                      mean='monthly_clm',st='2000-01',ed='2016-12')
# zonal_data_diff_Soda('u_eastward',uv_CMAP,uv_diff_CMAP,[-.08,.08],[-.05,.05],\
#                      mean='ann',st='1980-02',ed='2016-12')

# zonal_data_diff_Soda('v_northward',uv_CMAP,uv_CMAP,[-.01,.01],[-.01,.01],\
#                      mean='monthly_clm',st='2000-01',ed='2016-12')
# zonal_data_diff_Soda('v_northward',uv_CMAP,uv_diff_CMAP,[-.01,.01],[-.01,.01],\
#                      mean='ann',st='1980-02',ed='2016-12')

    
# zeta_AVISO(diff=True,{'1980-01','1990-12'})

# Z_SubT_drift({'1980-01','1990-12'})

# '''Auger Temp vertical section '''
Auger_temp_section('temp',temp_levels,temp_CMAP,has_year_zero=True,\
                    mean='-',st='0001-01',ed='0005-12')
# Auger_temp_section('temp',temp_levels,temp_CMAP,\
#                     mean='ann',st='1980-02',ed='2016-12')

# ''' Linear trend '''
# Surface_data_trend('zeta',[-80,-24],CMAP_LT,st='0001-01',ed='0120-12')

# zonal_data_trend('salt',[-80,-23],CMAP_LT,st='0005-01',ed='0005-12')


# =============================================================================



''' Sv & Current Drake passage (66W) '''

''' Sv & Current South of Africa (24.5E) '''

''' Sv & Current South of Australia (129E) '''

''' Meridional overturing streamfunction (zonal mean)'''


''' WOA '''
# 













