# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 18:28:17 2023

@author: shjo9
"""



PKG_path = 'D:/OneDrive/JNUpack/JNUROMS/'
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
from matplotlib import colors,colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap,LinearSegmentedColormap


My_Grd='D:/OneDrive/base142/Warehouse01/Grd_SO_05d_sponge.nc'
AVG_PATH='G:/Models/TK0525ED_CLM'
    
w_path_sig='/data4/base158/Warehouse01/MyPro01/TK0525EM/Vertical_temp/'
save_name='Vertical_temp_lon60'
fig_bool=0

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
t_rng=['1980-02','2014-12']
OGCM_TIMES=MFDataset(list(AVGS))['ocean_time']
My_time_ref=OGCM_TIMES.units
# OGCM_TIMES=MFDataset(OGCM_PATH+'*.nc')['ocean_time']

TIME_UNIT=OGCM_TIMES.units
OGCM_times=num2date(OGCM_TIMES[:],TIME_UNIT)
Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),28)
TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
# =============================================================================
VAR1,VAR2=[],[]
for i in tqdm(range(len(TIMES_co))):
    ncA,ncG=Dataset(AVGS[i]),Dataset(My_Grd)
    X,Z,tmp1=jr.get_section(ncG,ncA,'temp',[0,360],[-60])
    # _,_,tmp2=jr.get_section(ncG,ncA,'rho',[60],[-70,-34])
    ncA.close(); ncG.close()
    VAR1.append(tmp1); #VAR2.append(tmp2)
VAR1,VAR2=np.array(VAR1),np.array(VAR2)


ncA['s_w']

plt.rcParams["font.weight"] = "regular"
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
           '#008DDE','#07ABFF','#3EC1FF','#87D9FF','#EEEEEE']
MyCmap=ListedColormap(cmap_list).reversed()
# norm = colors.BoundaryNorm( np.arange(-1,26.5,.1), 32)
Label_size=13

# for i,j in zip(VAR1,OGCM_times):
#     Title_name=j.strftime('%Y-%m')

#     # Figures
#     fig, axs = plt.subplots(2,1,figsize=(6,4),constrained_layout = True,
#                             sharex=True,gridspec_kw={'height_ratios': [1, 1.3]},dpi=200)
#     axs[0].set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
#     im0=axs[0].contour(X,Z,i,colors='k',levels=[1,3])
#     im0.collections[1].set_linestyle('dashed')
#     im1=axs[0].pcolor(X,Z,i,cmap=MyCmap,edgecolors='k', linewidths=1)
#     axs[0].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
#     axs[0].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size, color='k',right=True)
#     axs[0].set_ylim(-Tcline,0)
#     axs[0].set_xlim(-70,-35)
#     im3=axs[1].contour(X,Z,i,vmin=-2,vmax=6,colors='k',levels=[1,3])
#     im3.collections[1].set_linestyle('dashed')
#     # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
#     im2=axs[1].pcolor(X,Z,i,cmap=MyCmap,vmin=-2,vmax=24,edgecolors='k', linewidths=1)
#     axs[1].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
#     axs[1].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size,  color='k',right=True)
#     axs[1].set_ylim(-5000,-Tcline)
#     axs[1].set_xlim(-70,-35)
#     divider = make_axes_locatable(axs[1])
#     cax = divider.append_axes("bottom", size="7%", pad=.35)
#     cax.tick_params(labelsize=Label_size)
#     cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
#     h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
#     # plt.tight_layout(pad=0.1)
#     if fig_bool:
#         plt.savefig(w_path_sig+'ppt/'+save_name+'_'+Title_name.replace('-','_'),
#                     facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
#         plt.savefig(w_path_sig+save_name+'_'+Title_name.replace('-','_'))
#     plt.show()

for i,j in zip(VAR1,OGCM_times):
    Title_name=j.strftime('%Y-%m')

    # Figures
    fig, axs = plt.subplots(2,1,figsize=(6,4),constrained_layout = True,
                            sharex=True,gridspec_kw={'height_ratios': [1, 1.3]},dpi=200)
    axs[0].set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
    im0=axs[0].contour(X,Z,i,colors='k',levels=[1,3])
    im0.collections[1].set_linestyle('dashed')
    im1=axs[0].pcolor(X,Z,i,cmap=MyCmap,edgecolors='k', linewidths=1)
    axs[0].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[0].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size, color='k',right=True)
    axs[0].set_ylim(-Tcline,0)
    axs[0].set_xlim(0,360)
    im3=axs[1].contour(X,Z,i,vmin=-2,vmax=6,colors='k',levels=[1,3])
    im3.collections[1].set_linestyle('dashed')
    # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
    im2=axs[1].pcolor(X,Z,i,cmap=MyCmap,vmin=-2,vmax=24,edgecolors='k', linewidths=1)
    axs[1].tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs[1].tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size,  color='k',right=True)
    axs[1].set_ylim(-5000,-Tcline)
    axs[1].set_xlim(0,360)
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("bottom", size="7%", pad=.35)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    # plt.tight_layout(pad=0.1)
    if fig_bool:
        plt.savefig(w_path_sig+'ppt/'+save_name+'_'+Title_name.replace('-','_'),
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path_sig+save_name+'_'+Title_name.replace('-','_'))
    plt.show()




    Title_name=j.strftime('%Y-%m')

    # Figures
    fig, axs = plt.subplots(1,1,figsize=(6,4),constrained_layout = True,
                            sharex=True,gridspec_kw={'height_ratios': [1]},dpi=200)
    axs.set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
    im0=axs.contour(X,Z,i,colors='k',levels=[1,3])
    im0.collections[1].set_linestyle('dashed')
    im1=axs.pcolor(X,Z,i,cmap=MyCmap,edgecolors='k', linewidths=1)
    axs.tick_params(axis='x', direction='in', length=3, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs.tick_params(axis='y', direction='in', length=3, pad=8, labelsize=Label_size, color='k',right=True)
    axs.set_ylim(-5000,0)
    axs.set_xlim(0,360)

    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("bottom", size="7%", pad=.35)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    # plt.tight_layout(pad=0.1)
    if fig_bool:
        plt.savefig(w_path_sig+'ppt/'+save_name+'_'+Title_name.replace('-','_'),
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path_sig+save_name+'_'+Title_name.replace('-','_'))
    plt.show()










