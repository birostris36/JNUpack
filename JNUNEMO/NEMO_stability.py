import os; import sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from libs.stb_master import Stb
from libs.plot import plot_2d,plot_2df
from netCDF4 import Dataset
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

rpth='E:/_tmp/EAST_NX_result/'
nc_name='EAST_NX_1d_20220112_20220119_grid_T.nc'
wpth=rpth+'tmp/'

### mkdir ====================================
try:
    os.mkdir(rpth+'tmp')
    os.mkdir(rpth+'tmp/ppt')
except:
    pass
    

My_stb=Stb(rpth=rpth,wpth=wpth)
My_stb.pro_stat1()
My_stb.pro_stat2()

### Plot output.abort =========================
# NC=Dataset(rpth+'output.abort.nc')
# Header=NC.variables.keys()
# lon,lat=NC['nav_lon'][:],NC['nav_lat'][:]

Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba','#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8','#e3f4fb',\
 '#f2f9e3','#fcf0b4','#fddb81','#fdc152','#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27','#ae191f']      
Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=22)
# CMAP = ListedColormap(Mycmap)

# ARGS={'lon':lon,'lat':lat,
#     'data':NC['sossheig'][0],
#     'Mylim':[-1.3,1.3],
#     'CMAP':Mycmap,
#     'Title_name':'Abort ssh',
#     'wpth':wpth,
#     'wname':'abort_ssh',
#     'fig_bool':True
#     }

# plot_2d(ARGS)


### Plot grid T ==================================================
zeta_lim=[-1.,1.]
zeta_levels=np.arange(zeta_lim[0],zeta_lim[-1]+0.2/2,0.2)
zeta_CMAP = ListedColormap(LinearSegmentedColormap.from_list('',Mycolorlist,N=256)(
    np.linspace(0, 1, len(zeta_levels)+1,endpoint=True)) )

sst_lim=[22,28]
sst_levels=np.arange(sst_lim[0],sst_lim[-1]+0.2/2,0.2)
sst_CMAP = ListedColormap(plt.get_cmap('jet')(
    np.linspace(0, 1, len(sst_levels)+1,endpoint=True)) )

NCT=Dataset(rpth+nc_name)
HeaderT=NCT.variables.keys()
lon,lat=NCT['nav_lon'][:],NCT['nav_lat'][:]

try:
    os.mkdir(wpth+'sst')
    os.mkdir(wpth+'ssh')
    os.mkdir(wpth+'sst/ppt/')
    os.mkdir(wpth+'ssh/ppt/')
except :
    pass

i=-1
while i<7-1:
    i+=1
    
    tmp_ssh,tmp_sst=NCT['sossheig'][i],NCT['sosstsst'][i]
    
    tmp_ssh[tmp_ssh<zeta_lim[0]]=zeta_lim[0]
    tmp_ssh[tmp_ssh>zeta_lim[-1]]=zeta_lim[-1]
    
    tmp_sst[tmp_sst<sst_lim[0]]=sst_lim[0]
    tmp_sst[tmp_sst>=sst_lim[-1]]=sst_lim[-1]

    ARGS_ssh={'lon':lon,'lat':lat,
        'data':tmp_ssh,
        'Mylim':zeta_lim,
        'levels':zeta_levels,
        'CMAP':zeta_CMAP,
        'Title_name':'ssh (time step: '+f'{i:02d})',
        'wpth':wpth+'ssh/',
        'wname':'ssh_'+f'{i:02d}',
        'fig_bool':True
        }

    ARGS_sst={'lon':lon,'lat':lat,
        'data':tmp_sst,
        'Mylim':sst_lim,
        'levels':sst_levels,
        'CMAP':sst_CMAP,
        'Title_name':'sst (time step: '+f'{i:02d})',
        'wpth':wpth+'sst/',
        'wname':'sst_'+f'{i:02d}',
        'fig_bool':True
        }
    
    plot_2df(ARGS_ssh)
    plot_2d(ARGS_sst)
    






