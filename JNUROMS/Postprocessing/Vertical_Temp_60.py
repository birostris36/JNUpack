# -*- coding: utf-8 -*-
"""
Created on Tue May  2 05:25:29 2023

@author: shjo9
"""


PKG_path = 'D:/OneDrive/JNUpack/JNUROMS/Tools/'
import sys 
sys.path.append(PKG_path)
import JNUROMS as jr
from JNU_create import create_ini
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
from scipy.interpolate import griddata
import datetime as dt
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import colors,colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import cmocean
import pandas as pd



My_Grd='D:/OneDrive/base142/Warehouse01/Grd_SO_05d.nc'
AVG_PATH='G:/Models/TK0525ED_CLM'
    
w_path_sig='D:/OneDrive/base142/Model_results/TK0525ED_CLM/vertical_temp_60/'
save_name='Vertical_temp_60'
fig_bool=1



AVGS=np.sort([AVG_PATH+'/'+i for i in os.listdir(AVG_PATH) if i.endswith('.nc')])

# Get My Grid info
ncG=Dataset(My_Grd)
lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
ncG.close()

atG,onG=lonG.shape
cosa,sina=np.cos(angle),np.sin(angle)

# Get OGCM Grid info
Sample_ncO=Dataset(AVGS[0])
lonO,latO=Sample_ncO['lon_rho'][:],Sample_ncO['lat_rho'][:]
depthO=Sample_ncO['s_rho'][:]
Tcline=Sample_ncO['hc'][:]
Sample_ncO.close()

# =============================================================================
# Process Times
x_rng=[-70,-24]
VAR1,VAR2=[],[]
for i in tqdm(range(len(AVGS))):
    ncA,ncG=Dataset(AVGS[i]),Dataset(My_Grd)
    for n in [0,1,2,3,4,5,6,7,8,9,10,11]:
        X,Z,tmp1=jr.get_section(ncG,ncA,'temp',[60],x_rng,tindx=n)
        # _,_,tmp2=jr.get_section(ncG,ncA,'rho',[60],[-70,-34])
        VAR2.append(tmp1)
    tmp_JFD=np.mean(np.array(VAR2),0)
    VAR2=[]
    ncA.close(); ncG.close()
    VAR1.append(tmp_JFD);
VAR1=np.array(VAR1)

My_date=pd.date_range('1980-01','2018-01',freq='Y').strftime('%Y')


plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams["font.family"] = 'Arial'
cmap_list=['#333333','#7F00BF','#9200E4','#AD07FF','#C23EFF','#DA87FF','#BC0D0F','#B30000','#C91F25',\
   '#D9343E','#E84A56','#F65E67','#FF6E6E','#FF8051','#FF8B1B','#FF9C00','#FFAA09','#FFBC1E',\
       '#FFD039','#FFE256','#FFF26F','#008000','#00A400','#00D500','#1EF31E','#69FC69','#0077B3',\
           '#008DDE','#07ABFF','#3EC1FF','#87D9FF','#EEEEEE']
MyCmap=ListedColormap(cmap_list).reversed()

NN=31
data_lim=[-2,24]

My_levels=np.linspace(data_lim[0],data_lim[-1],NN)
MyCmap = ListedColormap(MyCmap(np.linspace(0, .95, len(My_levels)+1,endpoint=True)))

# MyCmap = ListedColormap(cmocean.cm.thermal(np.linspace(0, 1, len(My_levels)+1,endpoint=True)))

Label_size=13
t=0
for i,j in zip(VAR1,My_date):
    t+=1
    Title_name='Years: '+j+f' (+{t-1:02d}~{t:02d})'

    # Figures
    fig, axs = plt.subplots(2,1,figsize=(6,4),constrained_layout = True,
                            sharex=True,gridspec_kw={'height_ratios': [1, 1.3]},dpi=200)
    axs[0].set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
    im0=axs[0].contour(X,Z,i,colors='k',levels=[1,3])
    im0.collections[1].set_linestyle('dashed')
    im1=axs[0].pcolor(X,Z,i,cmap=MyCmap,vmin=-2,vmax=24)
    axs[0].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[0].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size, color='k',right=True)
    axs[0].set_ylim(-Tcline,0)
    axs[0].set_xlim(-70,-24)
    im3=axs[1].contour(X,Z,i,vmin=-2,vmax=6,colors='k',levels=[1,3])
    im3.collections[1].set_linestyle('dashed')
    # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
    im2=axs[1].pcolor(X,Z,i,cmap=MyCmap,vmin=-2,vmax=24)
    axs[1].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[1].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size,  color='k',right=True)
    axs[1].set_ylim(-5000,-Tcline)
    axs[1].set_xlim(-70,-24)
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("bottom", size="7%", pad=.35)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    if fig_bool:
        plt.savefig(w_path_sig+'ppt/'+save_name+'_'+j.replace('-','_'),
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path_sig+save_name+'_'+j.replace('-','_'),bbox_inches='tight')
    plt.show()
















