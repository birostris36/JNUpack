# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:08:05 2023

@author: shjo9
"""
import sys 
sys.path.append('D:/OneDrive/JNUpack/')
sys.path.append('D:/OneDrive/JNUpack/JNUROMS')
import Tools.JNUROMS as jr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import numpy as np

Avg_npth='G:/Models/V205701_MW/Avg/'

ncG=Dataset('G:/MODEL_DATA/Grd/Grd_SO_05d_sponge.nc')
ncD=Dataset('G:/Models/V205701_MW/Avg/Avg_00001.nc')


X,Z,VAR=jr.get_section(ncG,ncD,'salt',[190],[-78,-24],tindx=4)

lat=X[0,:]

# lat_m,z_m=np.meshgrid(lat,Z[:,0])

# =============================================================================
# TMEP ===================================================================
temp_NN=15
temp_lim=[-2.5,28]
temp_levels=np.arange(temp_lim[0],temp_lim[-1]+2/2,2.)
temp_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')(np.linspace(0, 1, temp_NN,endpoint=True)))

temp_diff_lim=[-3,3]
temp_diff_levels=np.arange(temp_diff_lim[0],temp_diff_lim[-1]+.5/2,.5)

temp_CMAP = ListedColormap(plt.get_cmap('RdYlBu_r')(np.linspace(0, 1, temp_NN,endpoint=True)))

# SALT ===================================================================
salt_NN=15
salt_lim=[33,36.5]
salt_levels=np.arange(salt_lim[0],salt_lim[-1]+.2/2,.2)
salt_CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, len(salt_levels),endpoint=True)))

salt_diff_lim=[-.5,.5]
salt_diff_levels=np.arange(salt_diff_lim[0],salt_diff_lim[-1]+.05/2,.05)
salt_diff_CMAP = ListedColormap(plt.get_cmap('Spectral_r')(np.linspace(0, 1, len(salt_diff_levels),endpoint=True)))

# =============================================================================


data_lim=temp_lim
cmap=temp_CMAP
Label_size=12
xtick_location = np.linspace(lat[0], lat[-1],6)
xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]

fig, axs = plt.subplots(2,1,figsize=(6,4),
                        sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
# fig.subplots_adjust(wspace=0, hspace=0)
axs[0].set_title('temp',loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
im0=axs[0].contour(X,Z,VAR,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
axs[0].clabel(im0, inline=1, fontsize=14)
im0.collections[1].set_linestyle('dashed')
im1=axs[0].pcolor(X,Z,VAR,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
axs[0].set_ylim(-ncD['Tcline'][:],0)
axs[0].set_xlim(-80,-23.5)
im3=axs[1].contour(X,Z,VAR,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
axs[1].clabel(im0, inline=1, fontsize=14)
axs[0].set_xticks(ticks=xtick_location)
axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)

# im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
im2=axs[1].pcolor(X,Z,VAR,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
axs[1].set_ylim(-5000,-ncD['Tcline'][:])
axs[1].set_xlim(-80,-23.5)
axs[1].set_xticks(ticks=xtick_location)
axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("bottom", size="7%", pad=.35)
cax.tick_params(labelsize=Label_size)
cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
if 0:
    plt.savefig(save_pth+'Zonal_mean_'+kargs['mean']+'_'+data_nm+'/ppt/'+s_name_S,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(save_pth+'Zonal_mean_'+kargs['mean']+'_'+data_nm+'/'+s_name_S,bbox_inches='tight')
plt.show()
