# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 23:55:37 2023

@author: shjo9
"""

import sys
sys.path.append('D:/JNUpack/Mapping/Colorbars/bipolar/')
from bipolar import hotcold, bipolar
import os
import numpy as np
import xarray as xr
from eofs.xarray import Eof
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

My_str='PAC'

# from eofs.standard import Eof 
NN=10
pth='D:/HEAT/DATA/'
ncname='V_EN4_OHC_SO_c14_'+My_str+'_700m_1980_2023.nc'
w_path='D:/HEAT/EOF_V/'
Dir_pth='EOF_ohc_700m_1Y_'+My_str


try:
    os.mkdir(w_path+Dir_pth)
    os.mkdir(w_path+Dir_pth+'/ppt')
except:
    pass

CNN=16
fac=1
Mylim=[-0.1*fac,0.1*fac]
My_levels=np.linspace(Mylim[0], Mylim[-1], CNN+1,endpoint=True)

Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
 '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']      
Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)

CMAP = ListedColormap(Mycmap(
    np.linspace(0, 1, len(My_levels)-1,endpoint=True)) )

OHC = xr.open_dataset(pth+ncname).\
    loc[dict(time=slice('1980-01','2023-12'),lat=slice(-65,-10))].OHC

OHC_2Y=OHC.rolling(time=12,center=True).mean()[6:-5]

### EOF =======================================================================
solver=Eof(OHC_2Y)
eofs = -solver.eofs(neofs=NN, eofscaling=0)
pcs = -solver.pcs(npcs=NN,pcscaling=0)
var_=solver.varianceFraction(NN)*100
var=var_/np.sum(var_)*100

depth,lat=eofs.depth,eofs.lat
lat_m,depth_m=np.meshgrid(lat,depth)

### define =======================================================================
def plot_v_eofs(lat,depth,Mydata,t_name,CMAP,Mylim,My_levels,w_path,save_name,fig_bool=False):
    Label_size=18
    fig, axs = plt.subplots(1,1,figsize=(10.5,6),dpi=200)
    axs.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular','fontstyle':'italic'})
    # im0=axs.contour(lat_m,-depth_m,Mydata,colors='k',linestyle='-')
    # axs.clabel(im0, inline=1, fontsize=14)
    im1=axs.contourf(lat_m,-depth_m,Mydata,cmap=CMAP,vmin=Mylim[0],vmax=Mylim[-1],levels=My_levels)
    # im1=axs.pcolor(lat_m,-depth_m,Mydata,cmap=CMAP,vmin=Mylim[0],vmax=Mylim[-1])

    # im0.collections[1].set_linestyle('dashed')
    axs.tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    axs.tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
    # axs.set_ylim(-Tcline,0)
    # axs.set_xlim(-66.5,-42.5)
    # ax=gca()
    # ax.set_clim(-66.5,-42.5)

    divider = make_axes_locatable(axs)
    cax = divider.append_axes("bottom", size="7%", pad=.45)
    cax.tick_params(labelsize=Label_size)
    cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    h = fig.colorbar(im1, ax=axs,label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    h.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    if fig_bool:
        plt.savefig(w_path+Dir_pth+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+Dir_pth+'/'+save_name,bbox_inches='tight')
    plt.show()
    
def plot_pcs(time,time2,pc,t_name,w_path,save_name,fig_bool=True):
    Label_size = 18
    fig, axs = plt.subplots(1,1,figsize=(10,5),constrained_layout = True,
                        dpi=200)
    f1 = axs.plot(time,pc, label='KINETIC_ENRG',color='k',linewidth=2,zorder=0)
    axs.set_title(t_name,loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'italic'})
    axs.tick_params(axis='both', labelsize=Label_size)
    axs.grid(axis='x',linestyle='-.')
    xtick_location = time[5::12*4]
    xtick_labels = time2[5::12*4]
    axs.set_xticks(ticks=xtick_location)
    axs.set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+Dir_pth+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+Dir_pth+'/'+save_name,bbox_inches='tight')
    plt.show()
plt.rcParams["font.family"] = 'Arial'

TIME= [str(i)[0:7] for i in pcs.time.values]
TIME2=[str(i)[2:4] for i in pcs.time.values]

### Plot eofs =======================================================================
for i,j,n,m in zip(eofs[0:10].values*fac,np.arange(1,11),var.values,var_.values):
    save_name=Dir_pth+'_'+f'{j:02d}'+'mode'
    t_name=My_str+f' {j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'
    i[i>Mylim[-1]]=Mylim[-1]
    i[i<Mylim[0]]=Mylim[0]
    plot_v_eofs(lat_m,depth_m,i,t_name,CMAP,Mylim,My_levels,w_path,save_name,fig_bool=True)

### Plot pcs =======================================================================
for i,j,n,m in zip(pcs.values.transpose(),np.arange(1,11),var.values,var_.values):
    save_name='PC_'+f'{j:02d}'+'mode'
    t_name=My_str+f' {j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'
    plot_pcs(TIME,TIME2,i,t_name,w_path,save_name,fig_bool=True)


'''
NN=3
i1=-eofs[NN].values*fac
i1[i1>Mylim[-1]]=Mylim[-1]
i1[i1<Mylim[0]]=Mylim[0]
i2=-pcs.values.transpose()[NN]
j=np.arange(1,11)[NN]
n=var.values[NN]
m=var_.values[NN]
t_name=My_str+f' {j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'

plot_v_eofs(lat_m,depth_m,i1,t_name,CMAP,Mylim,My_levels,'w_path','save_name',fig_bool=False)
plot_pcs(TIME,TIME2,i2,t_name,'w_path','save_name',fig_bool=False)
'''

NN=1
i2=pcs.values.transpose()[NN]
j=np.arange(1,11)[NN]
n=var.values[NN]
m=var_.values[NN]
t_name=My_str+f' {j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'

plot_pcs(TIME,TIME2,i2,t_name,'w_path','save_name',fig_bool=False)


import pandas as pd
import matplotlib.pyplot as plt

ORI=pd.read_csv('D:/HEAT/Signals/MEI_ori.csv',header=None)
MEI_index=ORI.values.flatten()

MEI_index[MEI_index<-100]=np.nan
MEI_index=MEI_index[:-8]



mei_2Y=pd.DataFrame({'mei':MEI_index}).rolling(12,center=True).mean()
normal_mei=mei_2Y/np.max(mei_2Y)


normal_i=i2/np.max(i2)


np.corrcoef(normal_mei[6:-5].values.reshape(-1),normal_i.reshape(-1))

Label_size = 18
fig, axs = plt.subplots(1,1,figsize=(10,5),constrained_layout = True,
                    dpi=200)
f1 = axs.plot(TIME,-normal_i, label='pc 2 mode',color='k',linewidth=2,zorder=0)
f2 = axs.plot(TIME,normal_mei[6:-5], label='MEI index',color='r',linewidth=2,zorder=0)

axs.set_title(t_name,loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'italic'})
axs.tick_params(axis='both', labelsize=Label_size)
axs.grid(axis='x',linestyle='-.')
xtick_location = TIME[5::12*4]
xtick_labels = TIME2[5::12*4]
axs.set_xticks(ticks=xtick_location)
axs.set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
plt.legend(fontsize=16)
plt.tight_layout()
if 0:
    plt.savefig(w_path+Dir_pth+'/ppt/'+save_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+Dir_pth+'/'+save_name,bbox_inches='tight')
plt.show()









