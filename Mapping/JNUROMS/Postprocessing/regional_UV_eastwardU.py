# -*- coding: utf-8 -*-
"""
Created on Mon May  1 18:05:04 2023

@author: shjo9
"""



import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date, date2num,MFDataset
import os
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
sys.path.append('D:/OneDrive/JNUpack/')
from Mapping.Tools import d_modules as mm
from scipy.interpolate import interp2d, griddata
import matplotlib.ticker as mticker
import numpy.ma as ma
from matplotlib import ticker

AVGS_path='G:/Models/TK0525ED_CLM/'
save_path='D:/OneDrive/JNUpack/Mapping/Samples_figs/'

lat_rng,lon_rng=[-80,-60],[0,359]

AVGS=[AVGS_path+i for i in os.listdir(AVGS_path) if i.endswith('.nc')]

Sample_Data=Dataset(AVGS[-1])
lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
lon_co=np.where((lon_rho[0,:]>=lon_rng[0])&(lon_rho[0,:]<=lon_rng[-1]))[0]

lon_rho,lat_rho=np.meshgrid(lon_rho[0,lon_co],lat_rho[lat_co,0])

ZETA=Sample_Data['zeta'][0,lat_co,lon_co]

U=Sample_Data['u_eastward'][0,-1,lat_co,lon_co]
V=Sample_Data['v_northward'][0,-1,lat_co,lon_co]

lat_new = np.arange(lat_rng[0],lat_rng[-1], 2)
lon_new = np.arange(lon_rng[0],lon_rng[-1], 2)

lon_new_m,lat_new_m=np.meshgrid(lon_new,lat_new)

mask=np.invert(U.mask)

U_=griddata((lon_rho[mask].flatten(),lat_rho[mask].flatten()),U[mask].flatten(),
              (lon_new_m.flatten(),lat_new_m.flatten()),
           method='linear',fill_value=np.nan)
V_=griddata((lon_rho[mask].flatten(),lat_rho[mask].flatten()),V[mask].flatten(),
              (lon_new_m.flatten(),lat_new_m.flatten()),
           method='linear',fill_value=np.nan)
U_re = U_.reshape(lon_new_m.shape)
V_re = V_.reshape(lon_new_m.shape)

U_nega = ma.masked_array(U_re, mask=U_re>0)
U_posi = ma.masked_array(U_re, mask=U_re<=0)
V_nega = ma.masked_array(V_re, mask=U_re>0)
V_posi = ma.masked_array(V_re, mask=U_re<=0)

lat_nega=ma.masked_array(lat_new_m, mask=U_re<=0)
lon_nega=ma.masked_array(lon_new_m, mask=U_re<=0)
lat_posi=ma.masked_array(lat_new_m, mask=U_re>0)
lon_posi=ma.masked_array(lon_new_m, mask=U_re>0)

Speed_posi=(U_posi**2+V_posi**2)**(1/2)
Speed_nega=(U_nega**2+V_nega**2)**(1/2)



# UV={'u':U_re,'v':V_re,'lat':lat_new_m,'lon':lon_new_m}
# Scalar={'u':U,'lon':lon_rho,'lat':lat_rho}




plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams['font.family'] = 'Arial'


SIZE=(11,4)
FS=14
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt
import cmocean


NN=16
#My_levels=np.linspace(data_lim[0],data_lim[-1],NN+1,endpoint=True)

My_levels=np.arange(data_lim[0],data_lim[-1]+0.2/2,0.01)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
Nega_cmap = ListedColormap(cmocean.cm.balance(np.linspace(0, .4, len(My_levels)+1,endpoint=True)))
Posi_cmap = ListedColormap(cmocean.cm.balance(np.linspace(.6, 1, len(My_levels)+1,endpoint=True)))

CMAP=ListedColormap(cmocean.cm.balance(np.linspace(0, 1, len(My_levels)+1,endpoint=True)))


# Nega_cmap = ListedColormap(plt.get_cmap('PuOr_r')(np.linspace(0, .49, len(My_levels)+1,endpoint=True)))
# Posi_cmap = ListedColormap(plt.get_cmap('PuOr_r')(np.linspace(.5, 1, len(My_levels)+1,endpoint=True)))

# CMAP=plt.get_cmap('PuOr')
# zeta_CMAP=plt.get_cmap('jet',15)
# Plot_SO_Merc2(lon_rho,lat_rho,UV,ZETA,My_levels,zeta_CMAP,data_lim,save_path,'Merc_zeta_sample',fig_bool=False)

lonA=lon_rho
latA=lat_rho
# MyDATA_uv=UV
# MyDATA_scalar=Scalar
Mylim=[-.6,.6]
# Speed=(MyDATA_uv['u']**2+MyDATA_uv['v']**2)**(1/2)
# =============================================================================
# def
# =============================================================================
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

 # NORM
# def Plot_SO_Merc2(MyDATA_uv,MyDATA_scalar,My_levels,CMAP,Mylim,w_path,save_name,fig_bool=False):

#     Spheric=ccrs.SouthPolarStereo(central_longitude=180.0,globe=None)
#     PC = ccrs.PlateCarree(central_longitude=180.0,globe=None)
#     MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
#     # PC=MERC
#     # Now we will create axes object having specific projection 

#     fig, ax = plt.subplots(1, 1, figsize=(8,10),
#                        subplot_kw={'projection': PC},dpi=200)

#     gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
#                       linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
#                           )
        
#     gl.xlabels_top,gl.ylabels_right = False,False
#     gl.xlabel_style = gl.ylabel_style = {"size" : 18}
    
#     # ax.set_xticks(np.arange(1,359,10), crs=PC)
    
#     gl.xlocator = mticker.FixedLocator( np.arange(0,359,10) )
#     gl.xformatter = LONGITUDE_FORMATTER
#     # gl.xlocator = mticker.FixedLocator(np.arange(0,359,20))
    
#     # gl.xlocator = mticker.FixedLocator([0, 45, 180, 270, 359])
    
#     # To plot borders and coastlines, we can use cartopy feature
#     ax.add_feature(cf.COASTLINE.with_scale("50m"), lw=1,zorder=110)
#     ax.add_feature(cf.LAND.with_scale("50m"),color=[.75,.75,.75],zorder=100)

#     # q1 = plt.quiver(lon_posi,lat_posi,U_posi/Speed_posi,V_posi/Speed_posi,U_posi,
#     #     scale=36,headwidth=8.,headaxislength=10,headlength=13,color='k',
#     #     minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
#     #     pivot='mid',cmap=Posi_cmap,angles='xy')
#     # plt.clim(0,Mylim[-1])
    
#     # q2 = plt.quiver(lon_nega,lat_nega,U_nega/Speed_nega,V_nega/Speed_nega,U_nega,
#     #     scale=36,headwidth=8.,headaxislength=10,headlength=13,color='k',
#     #     minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
#     #     pivot='mid',cmap=Nega_cmap,angles='xy')
#     # plt.clim(Mylim[0],0)
    
    
#     q1 = plt.quiver(lon_posi,lat_posi,U_posi/Speed_posi,V_posi/Speed_posi,U_posi,
#         scale=30, width=0.0025, headwidth=8, headaxislength=3,color='k',
#         edgecolor='k',alpha=1.,transform=ccrs.PlateCarree(),zorder=1,minlength=1,
#         pivot='mid',cmap=Posi_cmap,angles='xy')
#     plt.clim(0,Mylim[-1])
    
#     q2 = plt.quiver(lon_nega,lat_nega,U_nega/Speed_nega,V_nega/Speed_nega,U_nega,
#         scale=30, width=0.0025, headwidth=8, headaxislength=3,color='k',
#         edgecolor='k',alpha=1.,transform=ccrs.PlateCarree(),zorder=1,minlength=1,
#         pivot='mid',cmap=Nega_cmap,angles='xy')
#     plt.clim(Mylim[0],0)
    
#     # M=plt.contourf(lonA,latA,ZETA,cmap=CMAP,levels=My_levels,transform=PC,zorder=0)
#     # M=plt.contourf(lonA,latA,ZETA,cmap=CMAP,transform=ccrs.PlateCarree(),zorder=0)

#     # M=plt.pcolormesh(MyDATA_scalar['lon'],MyDATA_scalar['lat'],MyDATA_scalar['u'],\
#     #                  cmap=CMAP,transform=PC,zorder=0)

#     #plt.pcolormesh(lonA, latA, MyDATA,
#     #              transform=PC,cmap=CMAP)
    
#     # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
#     ax.set_extent([lon_rng[0], lon_rng[-1], lat_rng[0], lat_rng[-1]], crs=ccrs.PlateCarree())
#     ax.tick_params(axis='both', which='major', labelsize=28)

#     divider = make_axes_locatable(ax)
#     ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

#     fig.add_axes(ax_cb)
#     cb=plt.colorbar(q1,extend='both',pad=0.01,cax=ax_cb)
#     cb.set_label(label='m', weight='regular',fontsize=28)
#     cb.ax.tick_params(labelsize=19)
    
#     plt.tight_layout()
#     if fig_bool:
#         plt.savefig(w_path+'ppt/'+save_name,
#                 facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
#         plt.savefig(w_path+save_name)
#     plt.show()

# =============================================================================
# Cartopy Mercator
# =============================================================================


def Plot_SO_Merc2(MyDATA_uv,MyDATA_scalar,My_levels,CMAP,Mylim,w_path,save_name,fig_bool=False):

    Spheric=ccrs.SouthPolarStereo(central_longitude=180.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=180.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # PC=MERC
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(8,10),
                       subplot_kw={'projection': PC},dpi=200)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
        
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 14}
    
    # ax.set_xticks(np.arange(1,359,10), crs=PC)
    
    gl.xlocator = mticker.FixedLocator( np.arange(0,359,20) )
    gl.xformatter = LONGITUDE_FORMATTER
    # gl.xlocator = mticker.FixedLocator(np.arange(0,359,20))
    
    # gl.xlocator = mticker.FixedLocator([0, 45, 180, 270, 359])
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("50m"), lw=1,zorder=110)
    ax.add_feature(cf.LAND.with_scale("50m"),color=[.75,.75,.75],zorder=100)

    q1 = plt.quiver(lon_posi,lat_posi,U_posi,V_posi,U_posi,
        scale=8,headwidth=8.,headaxislength=10,headlength=13,color='k',
        minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
        pivot='mid',cmap=Posi_cmap,angles='xy')
    plt.clim(0,Mylim[-1])
    
    q2 = plt.quiver(lon_nega,lat_nega,U_nega,V_nega,U_nega,
        scale=8,headwidth=8.,headaxislength=10,headlength=13,color='k',
        minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
        pivot='mid',cmap=Nega_cmap,angles='xy')
    plt.clim(Mylim[0],0)
    
    # plt.quiverkey(q1,123.8,-33.5,.5,"$5\cdot10^{-1} m/s$",coordinates='data',color='k',
    #             labelpos='E',alpha=1,labelcolor='k',fontproperties={'size':12},
    #             labelsep=0.13)
  
    # q1 = plt.quiver(lon_posi,lat_posi,U_posi,V_posi,U_posi,
    #     scale=8, width=0.0025, headwidth=8, headaxislength=3,color='k',
    #     edgecolor='k',alpha=1.,transform=ccrs.PlateCarree(),zorder=1,minlength=1,
    #     pivot='mid',cmap=Posi_cmap,angles='xy')
    # plt.clim(0,Mylim[-1])
    
    # q2 = plt.quiver(lon_nega,lat_nega,U_nega,V_nega,U_nega,
    #     scale=8, width=0.0025, headwidth=8, headaxislength=3,color='k',
    #     edgecolor='k',alpha=1.,transform=ccrs.PlateCarree(),zorder=1,minlength=1,
    #     pivot='mid',cmap=Nega_cmap,angles='xy')
    # plt.clim(Mylim[0],0)
    
    # M=plt.contourf(lonA,latA,ZETA,cmap=CMAP,levels=My_levels,transform=PC,zorder=0)
    # M=plt.contourf(lonA,latA,ZETA,cmap=CMAP,transform=ccrs.PlateCarree(),zorder=0)

    # M=plt.pcolormesh(MyDATA_scalar['lon'],MyDATA_scalar['lat'],MyDATA_scalar['u'],\
    #                  cmap=CMAP,transform=PC,zorder=0)

    #plt.pcolormesh(lonA, latA, MyDATA,
    #              transform=PC,cmap=CMAP)
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
    ax.set_extent([lon_rng[0], lon_rng[-1], lat_rng[0], lat_rng[-1]], crs=ccrs.PlateCarree())
    ax.tick_params(axis='both', which='major', labelsize=14)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(q1,extend='both',pad=0.01,cax=ax_cb)
    # cb.set_label(label='m', weight='regular',fontsize=14)
    # cb.ax.tick_params(labelsize=14)
    
    
    # cax2 = divider.new_horizontal(size="5%", pad=0.7, pack_start=True,axes_class=plt.Axes)
    cax2 = divider.new_horizontal(size="5%", pad=0.5,axes_class=plt.Axes)
    fig.add_axes(cax2)
    cb2 = fig.colorbar(q2,extend='both', cax=cax2)
    cb2.set_label(label='Eastward vel (m/s)', fontweight='bold',\
                  fontstyle='italic',fontsize=14,labelpad=13,fontname='Arial')
    cb2.ax.yaxis.set_ticks_position('right')
    
    
    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+save_name)
    plt.show()
