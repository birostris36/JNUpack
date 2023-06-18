# -*- coding: utf-8 -*-
"""
Created on Tue May  2 05:28:30 2023

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
    
w_path_sig='D:/OneDrive/base142/Model_results/TK0525ED_CLM/vertical_temp_140Auger/'
save_name='Auger_TK0525ED_CLM'
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
# =============================================================================
x_rng=[-66.5,-42.5]
VAR1,VAR2=[],[]
for i in tqdm(range(len(AVGS))):
    ncA,ncG=Dataset(AVGS[i]),Dataset(My_Grd)
    for n in [0,10,11]:
        X,Z,tmp1=jr.get_section(ncG,ncA,'temp',[140],x_rng,tindx=n)
        # _,_,tmp2=jr.get_section(ncG,ncA,'rho',[60],[-70,-34])
        VAR2.append(tmp1)
    tmp_JFD=np.mean(np.array(VAR2),0)
    VAR2=[]
    ncA.close(); ncG.close()
    VAR1.append(tmp_JFD);
VAR1=np.array(VAR1)


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
plt.rcParams['contour.negative_linestyle'] = 'solid'

cmap_list=['#7F00BF','#9200E4','#AD07FF','#C23EFF','#DA87FF','#BC0D0F','#B30000','#C91F25',\
   '#D9343E','#E84A56','#F65E67','#FF6E6E','#FF8051','#FF8B1B','#FF9C00','#FFAA09','#FFBC1E',\
       '#FFD039','#FFE256','#FFF26F','#008000','#00A400','#00D500','#1EF31E','#69FC69','#0077B3',\
           '#008DDE','#07ABFF','#3EC1FF','#87D9FF','#EEEEEE']
# MyCmap=ListedColormap(cmap_list).reversed()
# norm = colors.BoundaryNorm( np.arange(-3,28,.1), 900)
Label_size=18

data_lim=[-2.1,15]

NN=15

My_date=pd.date_range('1980-01','2018-01',freq='Y').strftime('%Y')

My_levels=np.linspace(data_lim[0],data_lim[-1],NN)
MyCmap = ListedColormap(cmocean.cm.thermal(np.linspace(0, 1, len(My_levels)+1,endpoint=True)))

# MyCmap=cmocean.cm.thermal

xtick_location=np.arange(-65,-43.5,2.5)
xtick_labels=[str(i)+'E' for i in xtick_location]
t=0
for i,j in zip(VAR1,My_date):
    t+=1
    Title_name='Years: '+j+f' (+{t-1:02d}~{t:02d})'
    # Figures
    fig, axs = plt.subplots(1,1,figsize=(10.5,4),constrained_layout = True,
                            sharex=True,gridspec_kw={'height_ratios': [1]},dpi=200)
    axs.set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
    im0=axs.contour(X,Z,i,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
    axs.clabel(im0, inline=1, fontsize=14)
    # im0.collections[1].set_linestyle('dashed')
    im1=axs.pcolor(X,Z,i,cmap=MyCmap,vmin=data_lim[0],vmax=data_lim[-1])
    axs.tick_params(axis='x', direction='in', length=5,width=1.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs.tick_params(axis='y', direction='in', length=5,width=1.5, pad=8, labelsize=Label_size, color='k',right=True)
    
    axs.set_xticks(ticks=xtick_location)
    axs.set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    
    axs.set_ylim(-800,0)
    axs.set_xlim(x_rng)
    
    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="3%", pad=.15)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs,label='',cax=cax,extend='both',aspect=.0001)
    if fig_bool:
        plt.savefig(w_path_sig+'ppt/'+save_name+'_'+j,
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path_sig+save_name+'_'+j,bbox_inches='tight')
    plt.show()






