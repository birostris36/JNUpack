# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 11:40:07 2023

@author: shjo9
"""


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
# from eofs.standard import Eof 
NN=10
pth='D:/HEAT/DATA/'
ncname='EN4_OHC_GLOBAL_c14_700m_1980_2023.nc'
w_path='D:/HEAT/EOF_H/'
# Dir_pth='EOF_ohc_700m'

# try:
#     os.mkdir(w_path+Dir_pth)
#     os.mkdir(w_path+Dir_pth+'/ppt')
# except:
#     pass

# CMAP=plt.get_cmap('RdYlBu_r',15)
CNN=14
fac=100
Mylim=[-0.04*fac,0.04*fac]
My_levels=np.linspace(Mylim[0], Mylim[-1], CNN+1,endpoint=True)
Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
 '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']      
Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)
CMAP = ListedColormap(Mycmap(
    np.linspace(0, 1, len(My_levels)-1,endpoint=True)) )

OHC = xr.open_dataset(pth+ncname).\
    loc[dict(time=slice('1980-01','2023-12'),lat=slice(-80,-10))].OHC
    
OHC_2Y=OHC.rolling(time=12,center=True).mean()[6:-5]

#lon_m,lat_m=np.meshgrid(OHC.lon,OHC.lat)
# coslat=np.cos(np.deg2rad(lat_m)) 
# wgts=np.sqrt(coslat)[:,np.newaxis] 
# solver=Eof(OHC.OHC,weights=wgts)

solver=Eof(OHC_2Y)
eofs = solver.eofs(neofs=NN, eofscaling=0)
pcs = solver.pcs(npcs=NN,pcscaling=0)

var_=solver.varianceFraction(NN)*100
var=var_/np.sum(var_)*100

def plot_pcs(time,time2,pc,t_name,w_path,save_name,fig_bool=True):
    Label_size = 18
    fig, axs = plt.subplots(1,1,figsize=(10,3.7),constrained_layout = True,
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

### Plot pcs ==================================================================
TIME= [str(i)[0:7] for i in pcs.time.values]
TIME2=[str(i)[2:4] for i in pcs.time.values]

def Plot_SO_Merc3(lonA,latA,MyDATA,t_name,CMAP,Mylim,My_levels,w_path,save_name,fig_bool=False):
    
    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,6),
                       subplot_kw={'projection': MERC},dpi=200)
    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 20}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular','fontstyle':'italic'})

    M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    # M=plt.pcolormesh(lonA, latA, MyDATA,
    #               transform=PC,cmap=CMAP)
    plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
 #   ax.set_extent([0, 360, -80, -24], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=24)
    cb.ax.tick_params(labelsize=19)
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+Dir_pth+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+Dir_pth+'/'+save_name,bbox_inches='tight')
    plt.show()
def plot_pcs(time,time2,pc,t_name,w_path,save_name,fig_bool=True):
    Label_size = 18
    fig, axs = plt.subplots(1,1,figsize=(10,3.7),constrained_layout = True,
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


import pandas as pd
import matplotlib.pyplot as plt

ORI=pd.read_csv('D:/HEAT/Signals/MEI_ori.csv',header=None)
MEI_index=ORI.values.flatten()

MEI_index[MEI_index<-100]=np.nan
MEI_index=MEI_index[:-8]

N=1
i=pcs.values.transpose()[N]
j=np.arange(1,11)[N]
n=var.values[N]
m=var_.values[N]
t_name='~2000m '+f'{j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'
plot_pcs(TIME,TIME2,-i,t_name,'w_path','save_name',fig_bool=False)



mei_2Y=pd.DataFrame({'mei':MEI_index}).rolling(12,center=True).mean()
normal_mei=mei_2Y/np.max(mei_2Y)

i_2Y=pd.DataFrame({'i':i}).rolling(12,center=True).mean()
normal_i_2Y=i_2Y/np.max(i_2Y)

normal_i=i/np.max(i)


np.corrcoef(normal_mei[6:-5].values.reshape(-1),normal_i.reshape(-1))

Label_size = 18
fig, axs = plt.subplots(1,1,figsize=(10,3.7),constrained_layout = True,
                    dpi=200)
f1 = axs.plot(TIME,normal_i, label='pc 3 mode',color='k',linewidth=2,zorder=0)
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







### Plot eof ==================================================================
    

    
N=3
i=eofs[N].values*fac
j=np.arange(1,11)[N]
n=var.values[N]
m=var_.values[N]
# t_name='~2000m '+f'{j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'
t_name='~700m '+f'{j:02d}'+' mode '+f'{n:.1f}'+'% ('+f'{m:.1f}'+'%)'
i[i>Mylim[-1]]=Mylim[-1]
i[i<Mylim[0]]=Mylim[0]
Plot_SO_Merc3(eofs.lon,eofs.lat,-i,t_name,CMAP,Mylim,My_levels,'w_path','save_name',fig_bool=False)











