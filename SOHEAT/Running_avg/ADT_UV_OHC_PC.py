import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.path as mpath
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle 
import pandas as pd
mpl.use('agg')

aviso_pth='E:/_data/AVISO/'
coord_pth='E:/HEAT/DATA/EOFs/'

Sample=xr.open_dataset(aviso_pth+'madt_h/dt_global_allsat_madt_h_y2002_m04.nc')

### load coords ===============================================================
with open(coord_pth+'iap_coord.pickle', 'rb') as f:
    iap_coords = pickle.load(f)
with open(coord_pth+'gec_coord.pickle', 'rb') as f:
    gec_coords = pickle.load(f)
iap_posi_1std,iap_posi_05std=iap_coords['iap_posi_1std'],iap_coords['iap_posi_05std']
iap_nega_1std,iap_nega_05std=iap_coords['iap_nega_1std'],iap_coords['iap_nega_05std']
gec_posi_1std,gec_posi_05std=gec_coords['gec_posi_1std'],gec_coords['gec_posi_05std']
gec_nega_1std,gec_nega_05std=gec_coords['gec_nega_1std'],gec_coords['gec_nega_05std']

gec_posi_1std=gec_posi_1std.where( (gec_posi_1std>='1980-01') & (gec_posi_1std<='2016-12') ).dropna()
gec_nega_1std=gec_nega_1std.where( (gec_nega_1std>='1980-01') & (gec_nega_1std<='2016-12') ).dropna()
gec_posi_05std=gec_posi_05std.where( (gec_posi_05std>='1980-01') & (gec_posi_05std<='2016-12') ).dropna()
gec_nega_05std=gec_nega_05std.where( (gec_nega_05std>='1980-01') & (gec_nega_05std<='2016-12') ).dropna()

iap_posi_1std=iap_posi_1std.where( (iap_posi_1std>='1980-01') & (iap_posi_1std<='2016-12') ).dropna()
iap_nega_1std=iap_nega_1std.where( (iap_nega_1std>='1980-01') & (iap_nega_1std<='2016-12') ).dropna()
iap_posi_05std=iap_posi_05std.where( (iap_posi_05std>='1980-01') & (iap_posi_05std<='2016-12') ).dropna()
iap_nega_05std=iap_nega_05std.where( (iap_nega_05std>='1980-01') & (iap_nega_05std<='2016-12') ).dropna()

posi_coord,nega_coord=iap_posi_05std,iap_nega_05std

lon_rng,lat_rng=[220-35,260+40],[-65-10,-45+30]

nc_ADT=xr.open_mfdataset(aviso_pth+'madt_h/*.nc').loc[dict(longitude=slice(lon_rng[0],lon_rng[-1]),latitude=slice(lat_rng[0],lat_rng[-1]))].adt
nc_U=xr.open_mfdataset(aviso_pth+'madt_u/*.nc').loc[dict(longitude=slice(lon_rng[0],lon_rng[-1]),latitude=slice(lat_rng[0],lat_rng[-1]))].ugos
nc_V=xr.open_mfdataset(aviso_pth+'madt_v/*.nc').loc[dict(longitude=slice(lon_rng[0],lon_rng[-1]),latitude=slice(lat_rng[0],lat_rng[-1]))].vgos
nc_GECOHC=xr.open_dataset('E:/HEAT/DATA/GECCO_OHC_SO_c14_700m_1980_2018.nc').loc[dict(lon=slice(lon_rng[0]-360,lon_rng[-1]-360),lat=slice(lat_rng[0],lat_rng[-1]))].OHC
# nc_GECOHC=xr.open_dataset('E:/HEAT/DATA/GECCO_OHC_SO_c14_700m_1980_2018.nc').loc[dict(lat=slice(lat_rng[0],lat_rng[-1]))].OHC

### detrend ======================================================================

NC_adt=nc_ADT.assign_coords({'TT':('time',range(len(nc_ADT.time)))})
NC_u  =nc_U.assign_coords({'TT':('time',range(len(nc_ADT.time)))})
NC_v  =nc_V.assign_coords({'TT':('time',range(len(nc_ADT.time)))})

NC_adt=NC_adt.swap_dims({"time":"TT"})
NC_u=NC_u.swap_dims({"time":"TT"})
NC_v=NC_v.swap_dims({"time":"TT"})

NC_s_adt=NC_adt.polyfit(dim='TT',deg=1,skipna=True)
NC_s_u=NC_u.polyfit(dim='TT',deg=1,skipna=True)
NC_s_v=NC_v.polyfit(dim='TT',deg=1,skipna=True)

fit_adt = xr.polyval(NC_adt.TT, NC_s_adt.polyfit_coefficients)
fit_u   = xr.polyval(NC_adt.TT, NC_s_u.polyfit_coefficients)
fit_v   = xr.polyval(NC_adt.TT, NC_s_v.polyfit_coefficients)

# Coef=NC_s.polyfit_coefficients[0]
# Coef_var=Coef.values*12 # (m/year)

adt_dt=NC_adt-fit_adt
adt_dt=adt_dt.swap_dims({"TT":"time"})
u_dt=NC_u-fit_u
u_dt=u_dt.swap_dims({"TT":"time"})
v_dt=NC_v-fit_v
v_dt=v_dt.swap_dims({"TT":"time"})

EKE = ((u_dt-u_dt.mean(dim='time'))**2+(v_dt-v_dt.mean(dim='time'))**2)/2

YY=12
ADT_dt_nY=adt_dt.rolling(time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]
U_dt_nY=u_dt.rolling(time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]
V_dt_nY=v_dt.rolling(time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]
GECOHC_nY=nc_GECOHC.rolling(time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]
EKE_1Y= EKE.rolling(time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]

EKE_1Y_a = EKE_1Y-EKE_1Y.mean(dim='time')

# print(np.nanmax(EKE_1Y_a[0]))
# print(np.nanmin(EKE_1Y_a[0]))

GECOHC_nY_a=GECOHC_nY-GECOHC_nY.mean(dim='time')

### Figure configure ================================

Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
 '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']     
Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)

lon_aviso,lat_aviso=np.meshgrid(ADT_dt_nY.longitude,ADT_dt_nY.latitude)
# lon_era,lat_era=np.meshgrid(SST_1p.longitude,SST_1p.latitude)
lon_gec,lat_gec=np.meshgrid(GECOHC_nY.lon,GECOHC_nY.lat)

print(lon_gec.shape)
print(GECOHC_nY.shape)


CNN=16

ohclim=[-2.,2.] ; adtlim=[-.16,.16] ;  uvlim=[-.06,.06]; ekelim=[-.03,.03]

ohc_levels=np.linspace(ohclim[0], ohclim[-1], CNN+1,endpoint=True)
adt_levels=np.linspace(adtlim[0], adtlim[-1], CNN+1,endpoint=True)
uv_levels= np.linspace(uvlim[0], uvlim[-1], CNN+1,endpoint=True)
eke_levels= np.linspace(ekelim[0], ekelim[-1], CNN+1+4,endpoint=True)

CMAP_adt = ListedColormap(Mycmap(
    np.linspace(0, 1, len(adt_levels)-1,endpoint=True)) )
CMAP_ohc = ListedColormap(Mycmap(
    np.linspace(0, 1, len(ohc_levels)-1,endpoint=True)) )
CMAP_uv = ListedColormap(Mycmap(
    np.linspace(0, 1, len(uv_levels)-1,endpoint=True)) )
CMAP_eke = ListedColormap(Mycmap(
    np.linspace(0, 1, len(eke_levels)-1,endpoint=True)) )

a,b=[220, 220],[-65,-45] ; c,d=[260, 260],[-65,-45]
e,f=[220, 260],[-65,-65] ; g,h=[220, 260],[-45,-45]

FMT='%.2f'

w_path='E:/_tmp/RMnT/GECOHC/'

Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)

COMPO_p=iap_posi_1std.strftime('%Y-%m')
COMPO_n=iap_nega_1std.strftime('%Y-%m')

### Plot figure ================================================================
DATA=GECOHC_nY_a.values*10**(-9); LEVELS=ohc_levels; CMAP=CMAP_ohc ; My_lim=ohclim
lon_m,lat_m= lon_gec,lat_gec ;#lon_aviso,lat_avios
# Index_nm=pd.Index([str(i)[:7] for i in ADT_dt_nY.time.values])
Index_nm=pd.Index([str(i)[:7] for i in GECOHC_nY_a.time.values])

for i,j in zip(DATA,Index_nm): # ADT_dt_nY EKE_1Y_a
    t_name='GECOHC 1Y '+j
    s_name='GECOHC 1Y_'+j.replace('-','_')

    if sum(COMPO_p == j):
        T_clr='r'
    elif sum(COMPO_n == j):
        T_clr='b'
    else:
        T_clr='k'
    
    i[i<=My_lim[0]]=My_lim[0] ; i[i>=My_lim[-1]]=My_lim[-1] ;
    
    fig, ax = plt.subplots(1, 1, figsize=(12.5,6),
                        subplot_kw={'projection': PC},dpi=200)
    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                        linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                            )
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 20}
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',color=T_clr,fontdict={'fontsize':24,'fontweight':'regular','fontstyle':'italic'})
    ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5)
    ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5)
    ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5)
    ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5)

    M=plt.contourf(lon_m,lat_m,i,cmap=CMAP,levels=LEVELS,transform=PC)
    # M=plt.pcolormesh(lon_aviso, lat_aviso, i,
    #            transform=PC,cmap=CMAP)
    # plt.clim(My_lim[0],My_lim[-1])

    ax.set_extent([lon_rng[0], lon_rng[-1], lat_rng[0], lat_rng[-1]], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=24)
    cb.ax.tick_params(labelsize=19)
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter(FMT))
    plt.tight_layout()
    if 1:
        plt.savefig(w_path+'/ppt/'+s_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+s_name,bbox_inches='tight')
    plt.show()

