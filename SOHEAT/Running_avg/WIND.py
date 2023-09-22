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
import sys
sys.path.append('C:/Users/shjo9/Bridge/JNUpack/JNUROMS/Tools/')
from Manta_WindStress import ra_windstrcurl
import cmocean

era_pth='E:/_data/ERA5_monthly_85/'
coord_pth='E:/HEAT/DATA/EOFs/'

# Sample=xr.open_dataset(era_pth+'madt_h/dt_global_allsat_madt_h_y2002_m04.nc')

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

posi_coord,nega_coord=iap_posi_05std,iap_nega_05std

lon_rng,lat_rng=[220-35,260+40],[-65-10,-45+30]

nc_U10=xr.open_dataset(era_pth+'ERA5_monthly_wind_1980_2021.nc').loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]),wind_time=slice('1993-01','2020-12'))].Uwind
nc_V10=xr.open_dataset(era_pth+'ERA5_monthly_wind_1980_2021.nc').loc[dict(lat=slice(lat_rng[0],lat_rng[-1]),lon=slice(lon_rng[0],lon_rng[-1]),wind_time=slice('1993-01','2020-12'))].Vwind
wind_lat,wind_lon=nc_U10.lat.values,nc_U10.lon.values

### Running anverage =========================================================
YY=12
U10_nY=nc_U10.rolling(wind_time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]
V10_nY=nc_V10.rolling(wind_time=YY,center=True).mean()[int(YY/2):int(-YY/2-1)]

u10_,v10_=U10_nY.values,V10_nY.values
t_,at_,on_=u10_.shape

wsc=np.load('E:/_tmp/RMnT/wsc.npy')*10**8
# wsc=np.zeros([t_,at_,on_])
# for i,j,n in zip(u10_,v10_,range(t_)):
#     wsc[n]=ra_windstrcurl(wind_lat,wind_lon,i,j)
#     print(n)
# np.save('E:/_tmp/RMnT/wsc.npy',wsc)
# raise

U10_nY_a=U10_nY-U10_nY.mean(dim='wind_time')
V10_nY_a=V10_nY-V10_nY.mean(dim='wind_time')

wsc_a=wsc-np.nanmean(wsc,axis=0)
# print(np.nanmin(wsc_a)*10**7)
# print(np.nanmax(wsc_a)*10**7)
# raise


### Figure configure ================================

# Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
#  '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']     
# Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)


lon_m,lat_m=np.meshgrid(wind_lon,wind_lat)


CNN=16

uvlim=[-.06,.06]; wsclim=[-5,5]

wsc_levels=np.linspace(wsclim[0], wsclim[-1], CNN+1,endpoint=True)
MyCmap = ListedColormap(cmocean.cm.curl(np.linspace(0, 1, len(wsc_levels)+1,endpoint=True)))
CMAP_wsc = ListedColormap(MyCmap(
    np.linspace(0, 1, len(wsc_levels)-1,endpoint=True)) )


a,b=[220, 220],[-65,-45] ; c,d=[260, 260],[-65,-45]
e,f=[220, 260],[-65,-65] ; g,h=[220, 260],[-45,-45]

FMT='%.2f'

w_path='E:/_tmp/RMnT/wsc/'

Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)

COMPO_p=iap_posi_1std.strftime('%Y-%m')
COMPO_n=iap_nega_1std.strftime('%Y-%m')
Index_nm=pd.Index([str(i)[:7] for i in V10_nY.wind_time.values])

### Plot figure ================================================================
DATA=wsc_a; LEVELS=wsc_levels; CMAP=CMAP_wsc ; My_lim=wsclim
for i,j in zip(DATA,Index_nm): # ADT_dt_nY EKE_1Y_a
    t_name='WSC 1Y '+j
    s_name='WSC 1Y_'+j.replace('-','_')

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

