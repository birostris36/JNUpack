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
era_pth='G:/ERA5_monthly_85/'

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
nc_MSL=xr.open_dataset(era_pth+'ERA5_SLP.nc').loc[dict(expver=1,latitude=slice(40,-80))].msl
nc_SST=xr.open_dataset(era_pth+'ERA5_SST.nc').loc[dict(expver=1,latitude=slice(40,-80))].sst
nc_EKE=xr.open_mfdataset('G:/AVISO/eke/*.nc').loc[dict(latitude=slice(-80,40))].eke
nc_ADT=xr.open_mfdataset('G:/AVISO/madt_h/*.nc').loc[dict(latitude=slice(-80,40))].adt
nc_U=xr.open_mfdataset('G:/AVISO/madt_u/*.nc').loc[dict(latitude=slice(-80,40))].ugos
nc_V=xr.open_mfdataset('G:/AVISO/madt_v/*.nc').loc[dict(latitude=slice(-80,40))].vgos
nc_GECOHC=xr.open_dataset('D:/HEAT/DATA/GECCO_OHC_SO_c14_700m_1980_2018.nc').loc[dict(lat=slice(-80,40))].OHC
#nc_GECOHC=xr.open_dataset('D:/HEAT/DATA/GECCO_OHC_SO_c14_700m_1980_2018.nc').loc[dict(lat=slice(-80,40))].OHC700

TIME=pd.date_range('1993-01',periods=len(nc_EKE.time),freq='m')

nc_EKE=nc_EKE.assign_coords({'TIME':('time',TIME )})
nc_EKE=nc_EKE.swap_dims({"time":"TIME"})
nc_ADT=nc_ADT.assign_coords({'TIME':('time',TIME )})
nc_ADT=nc_ADT.swap_dims({"time":"TIME"})
nc_U=nc_U.assign_coords({'TIME':('time',TIME )})
nc_U=nc_U.swap_dims({"time":"TIME"})
nc_V=nc_V.assign_coords({'TIME':('time',TIME )})
nc_V=nc_V.swap_dims({"time":"TIME"})

nc_SPEED=(nc_U**2+nc_V**2)**(1/2)

nc_MSL=nc_MSL-nc_MSL.mean(dim='time')
nc_SST=nc_SST-nc_SST.mean(dim='time')
nc_EKE=nc_EKE-nc_EKE.mean(dim='TIME')
nc_SPEED=nc_SPEED-nc_SPEED.mean(dim='TIME')
nc_ADT=nc_ADT-nc_ADT.mean(dim='TIME')
nc_U=nc_U-nc_U.mean(dim='TIME')
nc_V=nc_V-nc_V.mean(dim='TIME')
nc_GECOHC=nc_GECOHC-nc_GECOHC.mean(dim='time')

### Composite ===================================================================
MSL_1p,MSL_1n=nc_MSL.sel(time=posi_coord.strftime('%Y-%m')).mean(dim='time'),nc_MSL.sel(time=nega_coord.strftime('%Y-%m')).mean(dim='time')
SST_1p,SST_1n=nc_SST.sel(time=posi_coord.strftime('%Y-%m')).mean(dim='time'),nc_SST.sel(time=nega_coord.strftime('%Y-%m')).mean(dim='time')
GECOHC_1p,GECOHC_1n=nc_GECOHC.sel(time=posi_coord).mean(dim='time'),nc_GECOHC.sel(time=nega_coord).mean(dim='time')

EKE_1p,EKE_1n=nc_EKE.sel(TIME=posi_coord).mean(dim='TIME'),nc_EKE.sel(TIME=nega_coord).mean(dim='TIME')
ADT_1p,ADT_1n=nc_ADT.sel(TIME=posi_coord).mean(dim='TIME'),nc_ADT.sel(TIME=nega_coord).mean(dim='TIME')
U_1p,U_1n=nc_U.sel(TIME=posi_coord).mean(dim='TIME'),nc_U.sel(TIME=nega_coord).mean(dim='TIME')
V_1p,V_1n=nc_V.sel(TIME=posi_coord).mean(dim='TIME'),nc_V.sel(TIME=nega_coord).mean(dim='TIME')

SPEED_1p,SPEED_1n=nc_SPEED.sel(TIME=posi_coord).mean(dim='TIME'),nc_SPEED.sel(TIME=nega_coord).mean(dim='TIME')

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

lon_aviso,lat_aviso=np.meshgrid(ADT_1p.longitude,ADT_1p.latitude)
lon_era,lat_era=np.meshgrid(SST_1p.longitude,SST_1p.latitude)
lon_gec,lat_gec=np.meshgrid(GECOHC_1p.lon,GECOHC_1p.lat)

CNN=16

ohclim=[-1.6,1.6] ; msllim=[-500,500] ; sstlim=[-1.6,1.6]
adtlim=[-.16,.16] ; ekelim=[-150,150] ; uvlim=[-.06,.06]

ohc_levels=np.linspace(ohclim[0], ohclim[-1], CNN+1,endpoint=True)
msl_levels=np.linspace(msllim[0], msllim[-1], CNN+1,endpoint=True)
sst_levels=np.linspace(sstlim[0], sstlim[-1], CNN+1,endpoint=True)
adt_levels=np.linspace(adtlim[0], adtlim[-1], CNN+1,endpoint=True)
eke_levels=np.linspace(ekelim[0], ekelim[-1], CNN+1,endpoint=True)
uv_levels= np.linspace(uvlim[0], uvlim[-1], CNN+1,endpoint=True)

EKE_1p,EKE_1n=EKE_1p.values,EKE_1n.values
EKE_1p[EKE_1p>=eke_levels[-1]]=eke_levels[-1]
EKE_1p[EKE_1p<=eke_levels[0]]=eke_levels[-0]
EKE_1n[EKE_1n>=eke_levels[-1]]=eke_levels[-1]
EKE_1n[EKE_1n<=eke_levels[0]]=eke_levels[-0]

SPEED_1p,SPEED_1n=SPEED_1p.values,SPEED_1n.values
SPEED_1p[SPEED_1p>=uv_levels[-1]]=uv_levels[-1]
SPEED_1p[SPEED_1p<=uv_levels[0]]=uv_levels[-0]
SPEED_1n[SPEED_1n>=uv_levels[-1]]=uv_levels[-1]
SPEED_1n[SPEED_1n<=uv_levels[0]]=uv_levels[-0]

CMAP_sst = ListedColormap(Mycmap(
    np.linspace(0, 1, len(sst_levels)-1,endpoint=True)) )
CMAP_msl = ListedColormap(Mycmap(
    np.linspace(0, 1, len(msl_levels)-1,endpoint=True)) )
CMAP_adt = ListedColormap(Mycmap(
    np.linspace(0, 1, len(adt_levels)-1,endpoint=True)) )
CMAP_ohc = ListedColormap(Mycmap(
    np.linspace(0, 1, len(ohc_levels)-1,endpoint=True)) )
CMAP_eke = ListedColormap(Mycmap(
    np.linspace(0, 1, len(eke_levels)-1,endpoint=True)) )
CMAP_uv = ListedColormap(Mycmap(
    np.linspace(0, 1, len(uv_levels)-1,endpoint=True)) )

### SST
Plot_SO_Merc3(lon_era,lat_era,SST_1p,'Composite-1std P-SST',CMAP_sst,sstlim,sst_levels,'%.1f',w_path,'Composite_1std_P_SST',fig_bool=True)
Plot_SO_Merc3(lon_era,lat_era,SST_1n,'Composite-1std N-SST',CMAP_sst,sstlim,sst_levels,'%.1f',w_path,'Composite_1std_N_SST',fig_bool=True)

### MSL
Plot_SO_Merc3(lon_era,lat_era,MSL_1p,'Composite-1std P-MSL',CMAP_msl,msllim,msl_levels,'%.0f',w_path,'Composite_1std_P_MSL',fig_bool=True)
Plot_SO_Merc3(lon_era,lat_era,MSL_1n,'Composite-1std N-MSL',CMAP_msl,msllim,msl_levels,'%.0f',w_path,'Composite_1std_N_MSL',fig_bool=True)

### ADT
Plot_SO_Merc3(lon_aviso,lat_aviso,ADT_1p,'Composite-1std P-SLA',CMAP_adt,adtlim,adt_levels,'%.2f',w_path,'Composite_1std_P_SLA',fig_bool=True)
Plot_SO_Merc3(lon_aviso,lat_aviso,ADT_1n,'Composite-1std N-SLA',CMAP_adt,adtlim,adt_levels,'%.2f',w_path,'Composite_1std_N_SLA',fig_bool=True)

Plot_SO_Merc3(lon_gec,lat_gec,GECOHC_1p*10**(-9),'Composite-1std P-OHC',CMAP_ohc,ohclim,ohc_levels,'%.2f',w_path,'Composite_1std_P_OHC',fig_bool=True)
Plot_SO_Merc3(lon_gec,lat_gec,GECOHC_1n*10**(-9),'Composite-1std N-OHC',CMAP_ohc,ohclim,ohc_levels,'%.2f',w_path,'Composite_1std_N_OHC',fig_bool=True)

Plot_SO_Merc3(lon_aviso,lat_aviso,EKE_1p,'Composite-1std P-EKE',CMAP_eke,ekelim,eke_levels,'%.0f',w_path,'Composite_1std_P_EKE',fig_bool=True)
Plot_SO_Merc3(lon_aviso,lat_aviso,EKE_1n,'Composite-1std N-EKE',CMAP_eke,ekelim,eke_levels,'%.0f',w_path,'Composite_1std_N_EKE',fig_bool=True)

Plot_SO_Merc3(lon_aviso,lat_aviso,SPEED_1p,'Composite-1std P-Speed',CMAP_uv,uvlim,uv_levels,'%.2f',w_path,'Composite_1std_P_SPEED',fig_bool=True)
Plot_SO_Merc3(lon_aviso,lat_aviso,SPEED_1n,'Composite-1std N-Speed',CMAP_uv,uvlim,uv_levels,'%.2f',w_path,'Composite_1std_N_SPEED',fig_bool=True)


### PLot vector ====================================================================================
U_1p_,V_1p_=U_1p.loc[dict(latitude=slice(-65-10,-45+10),longitude=slice(220-10,260+10))],V_1p.loc[dict(latitude=slice(-65-10,-45+10),longitude=slice(220-10,260+10))]
U_1n_,V_1n_=U_1n.loc[dict(latitude=slice(-65-10,-45+10),longitude=slice(220-10,260+10))],V_1n.loc[dict(latitude=slice(-65-10,-45+10),longitude=slice(220-10,260+10))]
ADT_1p_=ADT_1p.loc[dict(latitude=slice(-65-10,-45+10),longitude=slice(220-10,260+10))]
ADT_1n_=ADT_1n.loc[dict(latitude=slice(-65-10,-45+10),longitude=slice(220-10,260+10))]
lon_uv,lat_uv=np.meshgrid(U_1p_.longitude,V_1p_.latitude)

from scipy.interpolate import griddata

lon_re,lat_re=np.arange(U_1p_.longitude[0],U_1p_.longitude[-1],.5),np.arange(U_1p_.latitude[0],U_1p_.latitude[-1],.5)
lon_uv_re,lat_uv_re=np.meshgrid(lon_re,lat_re)

U_1p_re=griddata( (lon_uv.flatten(),lat_uv.flatten()),U_1p_.values.flatten(),(lon_uv_re,lat_uv_re),'linear' )
V_1p_re=griddata( (lon_uv.flatten(),lat_uv.flatten()),V_1p_.values.flatten(),(lon_uv_re,lat_uv_re),'linear' )

U_1n_re=griddata( (lon_uv.flatten(),lat_uv.flatten()),U_1n_.values.flatten(),(lon_uv_re,lat_uv_re),'linear' )
V_1n_re=griddata( (lon_uv.flatten(),lat_uv.flatten()),V_1n_.values.flatten(),(lon_uv_re,lat_uv_re),'linear' )

U_1p_re[U_1p_re< (U_1p_re-U_1p_re.std())/2 ]=np.nan ; V_1p_re[V_1p_re< (V_1p_re-V_1p_re.std())/2 ]=np.nan
U_1n_re[U_1n_re< (U_1n_re-U_1n_re.std())/2 ]=np.nan ; V_1n_re[V_1n_re< (V_1n_re-V_1n_re.std())/2 ]=np.nan

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