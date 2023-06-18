# -*- coding: utf-8 -*-
"""
Created on Tue May  2 08:04:44 2023

@author: shjo9
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May  2 06:47:47 2023

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
from tqdm import tqdm
fig_bool=0


data_path='G:/SODA/'
w_path='D:/OneDrive/base142/Model_results/TK0525ED_CLM/Surface_UV_uvel_soda/'
save_name='Surface_SODA_UV_uvel'

lat_rng,lon_rng=[-80,-23],[110,300]

OGCM=[data_path+i for i in os.listdir(data_path) if i.endswith('.nc')]

Sample_Data=Dataset(OGCM[-1])
lat_rho,lon_rho=Sample_Data['yt_ocean'][:],Sample_Data['xt_ocean'][:]
lat_co=np.where((lat_rho[:]>=lat_rng[0])&(lat_rho[:]<=lat_rng[-1]))[0]
lon_co=np.where((lon_rho[:]>=lon_rng[0])&(lon_rho[:]<=lon_rng[-1]))[0]

lon_rho,lat_rho=np.meshgrid(lon_rho[lon_co],lat_rho[lat_co])

lat_new = np.arange(lat_rng[0],lat_rng[-1], 2)
lon_new = np.arange(lon_rng[0],lon_rng[-1], 2)
lon_new_m,lat_new_m=np.meshgrid(lon_new,lat_new)

Zeta,U_re,V_re=[],[],[]
for i in tqdm(range(len(OGCM))):
    
    tmp_zeta=Dataset(OGCM[i])['ssh'][:,lat_co,lon_co]
    
    tmp_u=Dataset(OGCM[i])['u'][:,0,lat_co,lon_co]
    tmp_v=Dataset(OGCM[i])['v'][:,0,lat_co,lon_co]
    mask=np.invert(tmp_u[0].mask)
    
    tmp_mz=np.mean(tmp_zeta,axis=0)
    
    tmp_mu=np.mean(tmp_u,axis=0)
    tmp_mv=np.mean(tmp_v,axis=0)
    tmp_U_=griddata((lon_rho[mask].flatten(),lat_rho[mask].flatten()),tmp_mu[mask].flatten(),
                  (lon_new_m.flatten(),lat_new_m.flatten()),
               method='linear',fill_value=np.nan)
    tmp_V_=griddata((lon_rho[mask].flatten(),lat_rho[mask].flatten()),tmp_mv[mask].flatten(),
                  (lon_new_m.flatten(),lat_new_m.flatten()),
               method='linear',fill_value=np.nan)
    tmp_U_re = tmp_U_.reshape(lon_new_m.shape)
    tmp_V_re = tmp_V_.reshape(lon_new_m.shape)
    U_re.append(tmp_U_re); V_re.append(tmp_V_re); Zeta.append(tmp_mz)
U_re=np.array(U_re); V_re=np.array(V_re);
    

U_nega = ma.masked_array(U_re, mask=U_re>0)
U_posi = ma.masked_array(U_re, mask=U_re<=0)
V_nega = ma.masked_array(V_re, mask=U_re>0)
V_posi = ma.masked_array(V_re, mask=U_re<=0)

lat_nega=ma.masked_array(lat_new_m, mask=U_re[0]<=0)
lon_nega=ma.masked_array(lon_new_m, mask=U_re[0]<=0)
lat_posi=ma.masked_array(lat_new_m, mask=U_re[0]>0)
lon_posi=ma.masked_array(lon_new_m, mask=U_re[0]>0)

Speed_posi=(U_posi**2+V_posi**2)**(1/2)
Speed_nega=(U_nega**2+V_nega**2)**(1/2)



# UV={'u':U_re,'v':V_re,'lat':lat_new_m,'lon':lon_new_m}
# Scalar={'u':U,'lon':lon_rho,'lat':lat_rho}
import pandas as pd
My_date=pd.date_range('1980-01','2018-01',freq='Y').strftime('%Y')


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

data_lim1=[-.05,.05]
data_lim2=[-2.5,1.5]

NN=16
#My_levels=np.linspace(data_lim[0],data_lim[-1],NN+1,endpoint=True)

My_levels=np.arange(data_lim1[0],data_lim1[-1]+0.2/2,0.01)
#zeta_CMAP = ListedColormap(cmocean.cm.balance(np.linspace(0, 1, NN,endpoint=True)))
Nega_cmap = ListedColormap(cmocean.cm.balance(np.linspace(0, .4, len(My_levels)+1,endpoint=True)))
Posi_cmap = ListedColormap(cmocean.cm.balance(np.linspace(.6, 1, len(My_levels)+1,endpoint=True)))




zeta_levels=np.arange(data_lim2[0],data_lim2[-1]+0.2/2,0.2)
CMAP=ListedColormap(cmocean.cm.balance(np.linspace(0, 1, len(zeta_levels)+1,endpoint=True)))
CMAP=ListedColormap(plt.get_cmap('jet')(np.linspace(0, 1, len(zeta_levels)+1,endpoint=True)))


# Nega_cmap = ListedColormap(plt.get_cmap('PuOr_r')(np.linspace(0, .49, len(My_levels)+1,endpoint=True)))
# Posi_cmap = ListedColormap(plt.get_cmap('PuOr_r')(np.linspace(.5, 1, len(My_levels)+1,endpoint=True)))

# CMAP=plt.get_cmap('PuOr')
# zeta_CMAP=plt.get_cmap('jet',15)
# Plot_SO_Merc2(lon_rho,lat_rho,UV,ZETA,My_levels,zeta_CMAP,data_lim,save_path,'Merc_zeta_sample',fig_bool=False)

lonA=lon_rho
latA=lat_rho
# MyDATA_uv=UV
# MyDATA_scalar=Scalar
Mylim=data_lim1
Mylim=[-.6,.6]
Label_size=10
# Speed=(MyDATA_uv['u']**2+MyDATA_uv['v']**2)**(1/2)
# =============================================================================
# def
# =============================================================================
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
t=0
for z,i,j,ii,jj,tt in zip(Zeta,U_posi,V_posi,U_nega,V_nega,My_date): 
    t+=1
    Title_name='Years: '+tt+f' (+{t-1:02d}~{t:02d})'

    Spheric=ccrs.SouthPolarStereo(central_longitude=180.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=180.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # PC=MERC
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(8,11),
                       subplot_kw={'projection': PC},dpi=200)
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,y_inline=False,x_inline=False,rotate_labels=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    ax.set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})

    gl.xlabels_top,gl.ylabels_right = False,False
    # gl.xlabel_style = gl.ylabel_style = {"size" : 14}
    
    # ax.set_xticks(np.arange(1,359,10), crs=PC)
    
    # gl.xlocator = mticker.FixedLocator( np.arange(0,359,20) )
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.xlocator = mticker.FixedLocator(np.arange(0,359,20))
    
    # gl.xlocator = mticker.FixedLocator([0, 45, 180, 270, 359])
    # ax.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("50m"), lw=1,zorder=110)
    ax.add_feature(cf.LAND.with_scale("50m"),color=[.75,.75,.75],zorder=100)

    q1 = plt.quiver(lon_posi,lat_posi,i,j,i,
        scale=8,headwidth=8.,headaxislength=10,headlength=13,color='k',
        minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
        pivot='mid',cmap=Posi_cmap,angles='xy')
    plt.clim(0,Mylim[-1])
    
    q2 = plt.quiver(lon_nega,lat_nega,ii,jj,ii,
        scale=8,headwidth=8.,headaxislength=10,headlength=13,color='k',
        minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=ccrs.PlateCarree(),zorder=1,
        pivot='mid',cmap=Nega_cmap,angles='xy')
    plt.clim(Mylim[0],0)
    
        
    plt.quiverkey(q1,-58.,-19.,.5,"$5\cdot10^{-1} m/s$",coordinates='data',color='r',
                labelpos='E',alpha=1.,labelcolor='k',fontproperties={'size':Label_size},
                labelsep=0.13,transform=ccrs.PlateCarree(),zorder=3)
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
    
    # m1=plt.contourf(lon_rho,lat_rho,z,cmap=CMAP,levels=zeta_levels,transform=ccrs.PlateCarree(),zorder=0,alpha=.5)
    # m1=plt.contour(lon_rho,lat_rho,z,colors='k',levels=zeta_levels,transform=ccrs.PlateCarree(),zorder=0,alpha=1)

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
                  fontstyle='italic',fontsize=Label_size,labelpad=13,fontname='Arial')
    cb2.ax.yaxis.set_ticks_position('right')
    
    # cax2 = divider.new_horizontal(size="5%", pad=0.7, pack_start=True,axes_class=plt.Axes)
    # cax3 = divider.new_horizontal(size="5%", pad=0.5,axes_class=plt.Axes,pack_start=True)
    # fig.add_axes(cax3)
    # cb3 = fig.colorbar(m1,extend='both', cax=cax3)
    # cb3.set_label(label='m', fontweight='bold',\
    #               fontstyle='italic',fontsize=14,labelpad=13,fontname='Arial')
    # cb3.ax.yaxis.set_ticks_position('left')
    
    # plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'ppt/'+save_name+'_'+tt,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+save_name+'_'+tt,bbox_inches='tight')
    plt.show()
