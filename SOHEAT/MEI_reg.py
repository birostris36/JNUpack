# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 12:41:48 2023

@author: shjo9
"""
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy
pth='D:/HEAT/DATA/'
ncname='GECCO_OHC_SO_c14_700m_1980_2018.nc'
ORI=pd.read_csv('D:/HEAT/Signals/MEI_ori.csv',header=None)

CNN=16
fac=10**9
Mylim=[-1*fac,1*fac]
My_levels=np.linspace(Mylim[0], Mylim[-1], CNN+1,endpoint=True)
Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
 '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']      
Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)
CMAP = ListedColormap(Mycmap(
    np.linspace(0, 1, len(My_levels)-1,endpoint=True)) )

### Read data index ===========================================================
MEI_index=ORI.values.flatten()
MEI_index[MEI_index<-100]=np.nan
MEI_index=MEI_index[:-8]

mei_2Y=pd.DataFrame({'mei':MEI_index}).rolling(12,center=True).mean()

OHC = xr.open_dataset(pth+ncname).\
    loc[dict(time=slice('1980-01','2023-12'),lat=slice(-80,-10))]
    
OHC_2Y=OHC.rolling(time=12,center=True).mean()    

### define ====================================================================
def linearRegress4Cube(sig,dataset,Slicing_date,method=1):
    '''
    Inputs :
        sig --> DataSeries with index
        dataset --> numpy array (3D [t,~,~])
        Slicing_date --> list / Slicing date ex) ['1993-01', '2018-01']

    -------
    Outputs :
        Coef :
        p_values :

    '''
    import pandas as pd
    from tqdm import tqdm
    import numpy as np
    from scipy import io
    import xarray as xr
    import statsmodels.api as sm
    from sklearn.linear_model import LinearRegression
    # from sklearn.datasets import fetch_20newsgroups_vectorized
    # from sklearn.feature_selection import chi2


    _,idx,idy = dataset.shape
    Coef,Intercept, p_values = np.zeros([idx,idy]),np.zeros([idx,idy]),np.zeros([idx,idy])
    for i in tqdm(range(idx)):
        for j in range(idy):
            if np.isnan(dataset[:,i,j].mean()):
                # print('dsadsadsadsadsadas')
                Coef[i,j] = np.nan
                Intercept[i,j] = np.nan
                continue

            tmp_dataset_sig = dataset[:,i,j]
            tmp_set = pd.DataFrame({'tmp_dataset_sig':tmp_dataset_sig,'sig':sig.values.reshape(-1)},index = sig.index)
            tmp = tmp_set.loc[Slicing_date[0]:Slicing_date[1]]
            X, X_data = tmp.sig.values.reshape(-1,1), tmp.tmp_dataset_sig.values.reshape(-1,1)

            if method == 'sm' :
                mod = sm.OLS(X_data,X)
                fii = mod.fit()
                Coef[i,j] = fii.summary2().tables[1]['Coef.'].x1        # print(line_fitter.coef_[1])
                p_values[i,j] = fii.summary2().tables[1]['P>|t|']

            else:
                line_fitter = LinearRegression()
                line_fitter.fit(X,X_data)
                Coef[i,j] = line_fitter.coef_
                Intercept[i,j] = line_fitter.intercept_

    return Coef, p_values

### Linear regression & p-values ==============================================
re=pd.DataFrame({'X':mei_2Y.values[6:-57].reshape(-1)},index=OHC_2Y.OHC[6:-5].time)
Coef,pvalue=linearRegress4Cube(re,OHC_2Y.OHC.values[6:-5],['1980-1','2018-12'],method='sm')
Coef,_=linearRegress4Cube(re,OHC_2Y.OHC.values[6:-5],['1980-1','2018-12'],method='1')

'''re=pd.DataFrame({'X':MEI_index.reshape(-1)},index=OHC.OHC.time)
Coef,pvalue=linearRegress4Cube(re,OHC.OHC.values,['1980-1','2023-12'],method='sm')
Coef,_=linearRegress4Cube(re,OHC.OHC.values,['1980-1','2023-12'],method='1')'''

# plt.pcolor(Coef,vmin=-2*10**9,vmax=2*10**9,cmap=CMAP)
# plt.colorbar()


p_values=deepcopy(pvalue)
p_values[p_values==0]=np.nan
for i in range(p_values.shape[0]):
    for j in range(p_values.shape[-1]):
        if np.abs(p_values[i,j]) <= 0.05:
            p_values[i,j]=1
        else:
            p_values[i,j]=np.nan


### PLot ======================================================================
Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
MERC=ccrs.Mercator(central_longitude=180.0,globe=None)

lon_m,lat_m=np.meshgrid(OHC.lon.values,OHC.lat.values)

t_name='~700m OHC regressed to MEI'

MyDATA=deepcopy(Coef)
MyDATA[MyDATA>Mylim[-1]]=Mylim[-1]
MyDATA[MyDATA<Mylim[0]]=Mylim[0]

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

M=plt.contourf(lon_m,lat_m,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
plt.contourf(lon_m,lat_m,p_values,3,hatches=['.'],transform=PC,alpha=0)

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
# cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.tight_layout()
if 0:
    plt.savefig(w_path+Dir_pth+'/ppt/'+save_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+Dir_pth+'/'+save_name,bbox_inches='tight')
plt.show()







