import matplotlib
#matplotlib.use('agg')
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
import pandas as pd
import pickle

w_path='D:/HEAT/figs/Composite/'

coord_pth='D:/HEAT/DATA/EOFs/'
era_pth='/home/shjo/Downloads/'

### load coords ===============================================================
with open(coord_pth+'iap_coord.pickle', 'rb') as f:
    iap_coords = pickle.load(f)
with open(coord_pth+'gec_coord.pickle', 'rb') as f:
    gec_coords = pickle.load(f)
iap_posi_1std,iap_posi_05std=iap_coords['iap_posi_1std'],iap_coords['iap_posi_05std']
iap_nega_1std,iap_nega_05std=iap_coords['iap_nega_1std'],iap_coords['iap_nega_05std']
gec_posi_1std,gec_posi_05std=gec_coords['gec_posi_1std'],gec_coords['gec_posi_05std']
gec_nega_1std,gec_nega_05std=gec_coords['gec_nega_1std'],gec_coords['gec_nega_05std']

gec_posi_1std=gec_posi_1std.where( (gec_posi_1std>='1993-01') & (gec_posi_1std<='2016-12') ).dropna()
gec_nega_1std=gec_nega_1std.where( (gec_nega_1std>='1993-01') & (gec_nega_1std<='2016-12') ).dropna()
gec_posi_05std=gec_posi_05std.where( (gec_posi_05std>='1993-01') & (gec_posi_05std<='2016-12') ).dropna()
gec_nega_05std=gec_nega_05std.where( (gec_nega_05std>='1993-01') & (gec_nega_05std<='2016-12') ).dropna()

iap_posi_1std=iap_posi_1std.where( (iap_posi_1std>='1993-01') & (iap_posi_1std<='2016-12') ).dropna()
iap_nega_1std=iap_nega_1std.where( (iap_nega_1std>='1993-01') & (iap_nega_1std<='2016-12') ).dropna()
iap_posi_05std=iap_posi_05std.where( (iap_posi_05std>='1993-01') & (iap_posi_05std<='2016-12') ).dropna()
iap_nega_05std=iap_nega_05std.where( (iap_nega_05std>='1993-01') & (iap_nega_05std<='2016-12') ).dropna()

posi_coord,nega_coord=iap_posi_1std,iap_nega_1std

### load nc ===================================================================
nc_SSR=xr.open_dataset(era_pth+'ERA5_.nc').loc[dict(expver=1,latitude=slice(40,-80))].ssr/86400

TIME=pd.date_range('1993-01',periods=len(nc_EKE.time),freq='m')

nc_SSR=nc_SSR-nc_SSR.mean(dim='time')

### Composite ===================================================================
SSR_1p,SSR_1n=nc_SSR.sel(time=posi_coord.strftime('%Y-%m')).mean(dim='time'),nc_SSR.sel(time=nega_coord.strftime('%Y-%m')).mean(dim='time')

### define  ===================================================================
def Plot_SO_Merc3(lonA,latA,MyDATA,t_name,CMAP,Mylim,My_levels,FMT,w_path,save_name,fig_bool=False):

    a,b=[220, 220],[-65,-45]
    c,d=[260, 260],[-65,-45]
    e,f=[220, 260],[-65,-65]
    g,h=[220, 260],[-45,-45]

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

    ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5)
    ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5)
    ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5)
    ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5)

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
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter(FMT))
    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name,bbox_inches='tight')
    plt.show()

Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
 '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']     
Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)

lon_era,lat_era=np.meshgrid(SSR_1p.longitude,SSR_1n.latitude)

CNN=16

# ohclim=[-1.6,1.6] ; msllim=[-500,500] ; sstlim=[-1.6,1.6]
# adtlim=[-.16,.16] ; ekelim=[-150,150] ; uvlim=[-.06,.06]

ssrlim=[0,500]

ssr_levels=np.linspace(ssrlim[0], ssrlim[-1], CNN+1,endpoint=True)

CMAP_msl = ListedColormap(Mycmap(
    np.linspace(0, 1, len(ssr_levels)-1,endpoint=True)) )


Plot_SO_Merc3(lon_era,lat_era,SSR_1p,'Composite-1std P-SSR',CMAP_msl,ssrlim,ssr_levels,'%.0f',w_path,'Composite_1std_P_MSL',fig_bool=True)
Plot_SO_Merc3(lon_era,lat_era,SSR_1n,'Composite-1std N-SSR',CMAP_msl,ssrlim,ssr_levels,'%.0f',w_path,'Composite_1std_N_MSL',fig_bool=True)

'''

a,b=[220, 220],[-65,-45]
c,d=[260, 260],[-65,-45]
e,f=[220, 260],[-65,-65]
g,h=[220, 260],[-45,-45]

Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
lon__,lat__=np.meshgrid(ADT_1p_.longitude,ADT_1p_.latitude)
# Now we will create axes object having specific projection 

fig, ax = plt.subplots(1, 1, figsize=(7,6.5),
                    subplot_kw={'projection': PC},dpi=200)
gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                    linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                        )
gl.xlabels_top,gl.ylabels_right = False,False
gl.xlabel_style = gl.ylabel_style = {"size" : 16}

# To plot borders and coastlines, we can use cartopy feature
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
ax.set_title('Composite-1std P-SLA & UV',loc='right',fontdict={'fontsize':16,'fontweight':'regular','fontstyle':'italic'})

ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5)

M=plt.contourf(lon__,lat__,ADT_1p_,cmap=CMAP_adt,levels=adt_levels,transform=PC)
# M=plt.pcolormesh(lonA, latA, MyDATA,
#               transform=PC,cmap=CMAP)
plt.clim(adtlim[0],adtlim[-1])

q1 = ax.quiver(lon_uv_re,lat_uv_re,U_1p_re,V_1p_re,
    scale=3.0,headwidth=8.,headaxislength=8,headlength=11,color='k',
    minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=PC,zorder=1,
    pivot='mid',angles='xy')
#q1 = ax.quiver(lon_uv_re,lat_uv_re,U_1p_re,V_1p_re,
#    scale=3.0,headwidth=6.,headaxislength=5,headlength=6,color='k',width=0.0025,
#    edgecolor='k',alpha=1.,transform=PC,zorder=1,
#    pivot='mid',angles='xy')

# crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
#   ax.set_extent([0, 360, -80, -24], crs=PC)
ax.tick_params(axis='both', which='major', labelsize=14)

divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
cb.set_label(label='', weight='regular',fontsize=14)
cb.ax.tick_params(labelsize=14)
cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.tight_layout()
if 1:
    plt.savefig(w_path+'/ppt/'+'Composite_1std_P_UV',
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+'/'+'Composite_1std_P_UV',bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(7,6.5),
                    subplot_kw={'projection': PC},dpi=200)
gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                    linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                        )
gl.xlabels_top,gl.ylabels_right = False,False
gl.xlabel_style = gl.ylabel_style = {"size" : 16}

# To plot borders and coastlines, we can use cartopy feature
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
ax.set_title('Composite-1std N-SLA & UV',loc='right',fontdict={'fontsize':16,'fontweight':'regular','fontstyle':'italic'})

ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5)

M=plt.contourf(lon__,lat__,ADT_1n_,cmap=CMAP_adt,levels=adt_levels,transform=PC)
# M=plt.pcolormesh(lonA, latA, MyDATA,
#               transform=PC,cmap=CMAP)
plt.clim(adtlim[0],adtlim[-1])

q1 = ax.quiver(lon_uv_re,lat_uv_re,U_1n_re,V_1n_re,
    scale=3.0,headwidth=8.,headaxislength=8,headlength=11,color='k',
    minlength=1,edgecolor='k',minshaft=1.3,alpha=1.,transform=PC,zorder=1,
    pivot='mid',angles='xy')

# crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
#   ax.set_extent([0, 360, -80, -24], crs=PC)
ax.tick_params(axis='both', which='major', labelsize=14)

divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

fig.add_axes(ax_cb)
cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
cb.set_label(label='', weight='regular',fontsize=14)
cb.ax.tick_params(labelsize=14)
cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.tight_layout()
if 1:
    plt.savefig(w_path+'/ppt/'+'Composite_1std_N_UV',
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+'/'+'Composite_1std_N_UV',bbox_inches='tight')
plt.show()
'''