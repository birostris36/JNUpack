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
# mpl.use('agg')

aviso_pth='E:/_data/AVISO/'
coord_pth='E:/HEAT/DATA/EOFs/'


# lon_rng,lat_rng=[220-35,260+80],[-65-10,-45+30]
lon_rng,lat_rng=[220,260],[-65,-45]

nc_GECOHC=xr.open_dataset('E:/HEAT/DATA/GECCO_OHC_SO_c14_700m_1980_2018.nc').loc[dict(lon=slice(lon_rng[0]-360,lon_rng[-1]-360),lat=slice(lat_rng[0],lat_rng[-1]))].OHC
# nc_GECOHC=xr.open_dataset('E:/HEAT/DATA/GECCO_OHC_SO_c14_700m_1980_2018.nc').loc[dict(lat=slice(lat_rng[0],lat_rng[-1]))].OHC
# nc_GECOHC[0].plot()
# plt.show()
target_GECOHC_ts=nc_GECOHC.mean(dim=['lat','lon']).values
np.save('E:/_tmp/RMnT/target_GECOHC_ts_1980_2018.npy',target_GECOHC_ts)

from scipy import io
io.savemat('E:/_tmp/RMnT/target_GECOHC_ts_1980_2018.mat',{'OHC':target_GECOHC_ts}) 
### detrend ======================================================================


