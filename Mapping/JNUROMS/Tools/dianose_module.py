# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:50:12 2023

@author: shjo9
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 15 11:11:48 2023

@author: shjo9
"""


import sys 
sys.path.append('D:/OneDrive/JNUpack/')
import Tools.JNUROMS as jr
from Mapping.Tools import d_modules as mm
import numpy as np
from netCDF4 import Dataset,MFDataset,date2num,num2date
import os
import xarray as xr
import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from scipy.interpolate import interp2d, griddata
import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime as dt
import cmocean

Avg_pth='G:/Models/TK0525ED_CLM/'
Log_npth='G:/TEST/LogTK25EM_DQD.txt'
save_pth='G:/TEST/'
Grd_npth='G:/MODEL_DATA/Grd/Grd_SO_05d_sponge.nc'
SODA_pth=''
fig_bool=0

plt.rcParams["font.weight"] = "regular"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.linewidth'] = 1.
# plt.rcParams['axes.grid'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams["font.family"] = 'Arial'
mpl.rcParams['axes.unicode_minus'] = False

# Weighted monthly mean
def season_mean(ds, calendar="standard"):
    
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = ds['ocean_time'].dt.days_in_month

    # Calculate the weights by grouping by 'time.season'
    weights = (
        month_length.groupby("ocean_time.season") / month_length.groupby("ocean_time.season").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby("ocean_time.season").sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby("ocean_time.season").sum(dim="ocean_time")
def yearly_mean(ds, calendar="standard"):
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = ds['ocean_time'].dt.days_in_month

    # Calculate the weights by grouping by 'time.season'
    weights = (
        month_length.groupby("ocean_time.year") / month_length.groupby("ocean_time.year").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    # np.testing.assert_allclose(weights.groupby("ocean_time.year").sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby("ocean_time.year").sum(dim="ocean_time")

# Evaluates model stability
def Stability01():
    global Avg_pth, Log_pth, save_pth,fig_bool
    
    with open(Log_pth) as f:
        Model_Log=f.readlines()
        Logs1000 = Model_Log[:1000]

    def get_model_stability(Model_Log):
        i=0;j=0; stid=10**5; 
        while i<=len(Model_Log):
            if i>=len(Model_Log):
                break
            elif Model_Log[i].find('rx0') !=-1:
                rx0=Model_Log[i].split(':')[-1].strip().replace(' ','')
                rx1=Model_Log[i+1].split(':')[-1].strip().replace(' ','')
                rx=rx0+'\n'+rx1
            if Model_Log[i].find('STEP')!=-1:
                header_name = [ii.strip() for ii in Model_Log[i].split(' ') if len(ii)]
                [header_name.append(ii.strip()) for ii in Model_Log[i+1].split('  ') if len(ii)]
                PD = pd.DataFrame({},columns=header_name)
                stid=i
            tmp1=Model_Log[i].lstrip()
            if i>stid and len(tmp1) and tmp1[0].isnumeric() and not i>len(Model_Log)-300:
                A= [ii.strip() for ii in tmp1.split(' ') if len(ii)]
                [A.append(ii.strip()) for ii in Model_Log[i+1].lstrip().split(' ') if len(ii)]
                PD.loc[j] = A
                i+=2; j+=1
            i+=1
        for i in PD.columns:
            try:
                PD[i] = PD[i].astype(float)
            except :
                pass
        return PD,rx
    
    PD,rx=get_model_stability(Model_Log)
    
    PD['NET_VOLUME'].plot()
    PD.columns
    
    RX=rx.split('\n')
    RX1=RX[0][:16]
    RX2=RX[-1][:-7]
    
    Title_name='Topo ('+RX1+' / '+RX2+')'
    Model_Times1 = PD['YYYY-MM-DD']
    Model_Times2 = [i[2:4] for i in Model_Times1]
    Label_size = 25
    
    fig, axs = plt.subplots(4,1,figsize=(11,8.5),constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': [1,1, 1.,1]},dpi=200)
    f1 = axs[0].plot(Model_Times1,PD['KINETIC_ENRG'], label='KINETIC_ENRG',color='k',linewidth=2,zorder=0)
    # axs[0].plot(np.array(Model_Times1)[t_co],f1_z2(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[0].tick_params(axis='y', labelsize=Label_size)
    # axs[0].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    xtick_location = Model_Times1[5::12]
    xtick_labels = Model_Times2[5::12]
    axs[0].set_xticks(ticks=xtick_location)
    axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
    # axs[0].grid(axis='x', alpha=.3,linestyle='-.',color='k')
    # axs[0].set_ylim(Zeta2d.data.mean()-Zeta2d.data.std()*2.5,\
    #                 Zeta2d.data.mean()+Zeta2d.data.std()*7)# axs[0].set_yticks(ticks=np.arange(18,23,1))
    axs[0].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[0].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
    axs[0].legend(fontsize=18)
    #! Fig2 
    f1 = axs[1].plot(Model_Times1,PD['POTEN_ENRG'], label='POTEN_ENRG',color='k',linewidth=2,zorder=0)
    # axs[1].plot(np.array(Model_Times1)[t_co],f1_m2(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[1].tick_params(axis='y', labelsize=Label_size)
    xtick_location = Model_Times1[5::12*10]
    xtick_labels =Model_Times2[5::12*10]
    axs[1].set_xticks(ticks=xtick_location)
    axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    # axs[1].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    # axs[1].set_ylim(Momentum2d.data.mean()-Momentum2d.data.std()*2.5,\
    #                 Momentum2d.data.mean()+Momentum2d.data.std()*3.)
    # axs[1].set_yticks(ticks=np.arange(18,23,1))
    axs[1].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[1].legend(fontsize=18,loc=4)
    
    f1 = axs[2].plot(Model_Times1,PD['TOTAL_ENRG'], label='TOTAL_ENRG',color='k',linewidth=2,zorder=0)
    # axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[2].tick_params(axis='y', labelsize=Label_size)
    xtick_location = Model_Times1[5::12*40]
    xtick_labels =Model_Times2[5::12*40]
    axs[2].set_xticks(ticks=xtick_location)
    axs[2].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    # axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    # axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
    #                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
    # axs[1].set_yticks(ticks=np.arange(18,23,1))
    axs[2].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[2].legend(fontsize=18,loc=4)
    
    f1 = axs[3].plot(Model_Times1,PD['NET_VOLUME'], label='NET_VOLUME',color='k',linewidth=2,zorder=0)
    # axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[3].tick_params(axis='y', labelsize=Label_size)
    xtick_location = Model_Times1[5::12*60]
    xtick_labels =Model_Times2[5::12*60]
    axs[3].set_xticks(ticks=xtick_location)
    axs[3].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    # axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    # axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
    #                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
    # axs[1].set_yticks(ticks=np.arange(18,23,1))
    axs[3].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[3].legend(fontsize=18,loc=4)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(save_pth+'ppt/'+'Model_momentum_logs',
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(save_pth+'Model_momentum_logs')
    plt.show()
def Stability02():
    global Avg_pth, Log_npth, save_pth,fig_bool

    ncs=[Avg_pth+i for i in os.listdir(Avg_pth)]

    SampleNC=Dataset(ncs[0])
    
    # Variables=SampleNC.variables.keys()
    
    PN,PM=SampleNC['pn'][:]*10**3,SampleNC['pm'][:]*10**3 # km
    CELL_size=1/PN*1/PM # m**2
    
    # ZETA
    ZETA=MFDataset(ncs)['zeta']
    ZETA_values=ZETA[:]
    ZETA_area=np.zeros_like(ZETA_values)
    n=0
    for i in ZETA_values:
        ZETA_area[n]=i*CELL_size
        n+=1
    
    Total_size=np.sum(CELL_size)

    ZETA_timeseries=np.sum(ZETA_area,axis=2).sum(axis=1)/Total_size
        
    # AICE
    AICE=MFDataset(ncs)['aice']
    AICE_values=AICE[:]
    AICE_area=np.zeros_like(AICE_values)
    n=0
    for i in AICE_values:
        AICE_area[n]=i*CELL_size
        n+=1
    
    AICE_timeseries=np.sum(AICE_area,axis=2).sum(axis=1)/Total_size
    
    t=np.arange(len(AICE_timeseries))
    
    AICE_fp1=np.polyfit(t[200:],AICE_timeseries[200:],1)
    AICE_trend=np.polyval(AICE_fp1,t[200:])
    
    # SHFLUX 
    SHFLUX=MFDataset(ncs)['shflux']
    SHFLUX_values=SHFLUX[:]
    SHFLUX_area=np.zeros_like(SHFLUX_values)
    n=0
    for i in SHFLUX_values:
        SHFLUX_area[n]=i*CELL_size
        n+=1
    
    SHFLUX_timeseries=np.sum(SHFLUX_area,axis=2).sum(axis=1)/Total_size
    
    
    t=range(len(SHFLUX_timeseries))
    
    SHFLUX_fp1=np.polyfit(t,SHFLUX_timeseries,1)
    SHFLUX_trend=np.polyval(SHFLUX_fp1,t)
    
    plt.plot(SHFLUX_timeseries)
    plt.plot(SHFLUX_trend)
    
    # SST
    SST=MFDataset(ncs)['temp'][:,-1,:,:]
    SST_values=SST[:]
    SST_area=np.zeros_like(SST_values)
    n=0
    for i in SST_values:
        SST_area[n]=i*CELL_size
        n+=1
    
    SST_timeseries=np.sum(SST_area,axis=2).sum(axis=1)/Total_size
    
    t=range(len(SST_timeseries))
    
    SST_fp1=np.polyfit(t,SST_timeseries,1)
    SST_trend=np.polyval(SST_fp1,t)
    
    
    Model_Times1 = t
    Model_Times2=[str(i/12)[:2] for i in t]

    Label_size = 25
    fig, axs = plt.subplots(4,1,figsize=(11,9.),constrained_layout = True,
                            sharex=True,gridspec_kw={'height_ratios': [1, 1.,1,1]},dpi=200)
    f1 = axs[0].plot(Model_Times1,ZETA_timeseries, label='ZETA',color='k',linewidth=2,zorder=0)
    # axs[0].plot(np.array(Model_Times1)[t_co],f1_z2(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[0].tick_params(axis='y', labelsize=Label_size)
    # axs[0].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    xtick_location = Model_Times1[5::12]
    xtick_labels = Model_Times2[5::12]
    axs[0].set_xticks(ticks=xtick_location)
    axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
    # axs[0].grid(axis='x', alpha=.3,linestyle='-.',color='k')
    # axs[0].set_ylim(Zeta2d.data.mean()-Zeta2d.data.std()*2.5,\
    #                 Zeta2d.data.mean()+Zeta2d.data.std()*7)# axs[0].set_yticks(ticks=np.arange(18,23,1))
    axs[0].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[0].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
    axs[0].legend(fontsize=18)
    axs[0].set_ylim(-0.4,0.3)
    #! Fig2 
    f1 = axs[1].plot(Model_Times1,AICE_timeseries, label='AICE',color='k',linewidth=2,zorder=0)
    # axs[1].plot(np.array(Model_Times1)[t_co],f1_m2(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[1].tick_params(axis='y', labelsize=Label_size)
    xtick_location = Model_Times1[5::12*10]
    xtick_labels =Model_Times2[5::12*10]
    axs[1].set_xticks(ticks=xtick_location)
    axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    # axs[1].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    # axs[1].set_ylim(Momentum2d.data.mean()-Momentum2d.data.std()*2.5,\
    #                 Momentum2d.data.mean()+Momentum2d.data.std()*3.)
    # axs[1].set_yticks(ticks=np.arange(18,23,1))
    axs[1].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[1].legend(fontsize=18,loc=4)
    
    f1 = axs[2].plot(Model_Times1,SHFLUX_timeseries, label='SHFLX',color='k',linewidth=2,zorder=0)
    # axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[2].tick_params(axis='y', labelsize=Label_size)
    xtick_location = Model_Times1[0::12*10]
    xtick_labels =Model_Times2[0::12*10]
    axs[2].set_xticks(ticks=xtick_location)
    axs[2].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    # axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    # axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
    #                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
    # axs[1].set_yticks(ticks=np.arange(18,23,1))
    axs[2].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[2].legend(fontsize=18,loc=4)
    axs[2].set_ylim(-100,100)
    
    f1 = axs[3].plot(Model_Times1,SST_timeseries, label='SST',color='k',linewidth=2,zorder=0)
    # axs[2].plot(np.array(Model_Times1)[t_co],f1_m3(tmp_t),color='r',linewidth=3,linestyle='dashed')
    axs[3].tick_params(axis='y', labelsize=Label_size)
    xtick_location = Model_Times1[0::12*10]
    xtick_labels =Model_Times2[0::12*10]
    axs[3].set_xticks(ticks=xtick_location)
    axs[3].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    # axs[2].set_xlim(Model_Times1.values[0],Model_Times1.values[-1])
    # axs[2].set_ylim(Momentum3d.data.mean()-Momentum3d.data.std()*2.5,\
    #                 Momentum3d.data.mean()+Momentum3d.data.std()*3.)
    # axs[1].set_yticks(ticks=np.arange(18,23,1))
    axs[3].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[3].legend(fontsize=18,loc=4)
    axs[3].set_ylim(7,12.6)
    
    plt.tight_layout()
    if fig_bool:
        plt.savefig(save_pth+'ppt/'+'Model_stability2',
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(save_pth+'Model_stability2')
    plt.show()

# Plot Surface results
def data_drift(data_nm,lat_rng,My_levels,cmap,data_lim,**kargs):
    global Avg_pth, Log_pth, save_pth,fig_bool

    AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]
    
    Sample_Data=Dataset(AVGS[0])
    lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
    lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
        
    Sample=xr.open_dataset(AVGS[0])
    # Sample.s_rho.values
    
    if [i for i in Sample[data_nm].coords].count('s_rho'):
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(s_rho=Sample.s_rho.values[-1],\
                                                        eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]
    else:
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]

    if kargs['mean']=='ann':
        data=yearly_mean(data_).rename({'year':'ocean_time'})
    elif kargs['mean']=='season':
        data=season_mean(data_).rename({'season':'ocean_time'})
    elif kargs['mean']=='monthly':
        data=data_.resample(ocean_time='1MS').mean()
    elif kargs['mean']=='monthly_clm':
        data=data_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
    elif kargs['mean']=='total':
        data=data_.mean(dim='ocean_time',keepdims=True)
        
    for i in data:
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)
        s_name_S='Spherical_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
        
        s_name_M=s_name_S.replace('Spherical','Merc')
        
        mm.Plot_SO_Spherical2(data.lon_rho.values,data.lat_rho.values,i.squeeze().values,\
                              t_name,My_levels,cmap,data_lim,save_pth,s_name_S,fig_bool)
            
        mm.Plot_SO_Merc2(lon_rho,lat_rho,i,t_name,My_levels,cmap,data_lim,\
                          save_pth,s_name_M,fig_bool)
            
# Auger Temp vertical section
def Auger_temp_section(data_nm,My_levels,cmap,**kargs):
    global Avg_pth,Grd_npth, save_pth,fig_bool
    
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams["font.weight"] = "regular"

    save_name='Auger_TK0525ED_CLM'
    
    AVGS=np.sort([Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')])
    
    # Process Times
    t_rng=[kargs['st'],kargs['ed']]
    OGCM_TIMES=MFDataset(Avg_pth+'*nc')['ocean_time']
    units=OGCM_TIMES.units
    print(11)
    OGCM_times=num2date(OGCM_TIMES[:],units)
    Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
    Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
    TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
    Avg_co=list(set(TIMES_co//12))
    print(22)

    My_date=list(np.sort(list(set([i.strftime('%Y') for i in OGCM_times[TIMES_co]]))))
    
    # Get My Grid info
    ncG=Dataset(Grd_npth)
    # lonG,latG=ncG['lon_rho'][:],ncG['lat_rho'][:]
    # angle,topo,mask=ncG['angle'][:],ncG['h'][:],ncG['mask_rho'][:]
    ncG.close()
    
    # atG,onG=lonG.shape
    # cosa,sina=np.cos(angle),np.sin(angle)
    
    # Get OGCM Grid info
    Sample_ncO=Dataset(AVGS[0])
    # lonO,latO=Sample_ncO['lon_rho'][:],Sample_ncO['lat_rho'][:]
    # depthO=Sample_ncO['s_rho'][:]
    Tcline=Sample_ncO['hc'][:]
    Sample_ncO.close()
    
    data_lim=[-2.1,15]
    My_levels=np.linspace(data_lim[0],data_lim[-1],15)
    MyCmap = ListedColormap(cmocean.cm.thermal(np.linspace(0, 1, len(My_levels)+1,endpoint=True)))
    cmap=MyCmap
    # =============================================================================
    # Process Times
    # =============================================================================
    x_rng=[-66.5,-42.5]
    VAR1,VAR2=[],[]
    for i in Avg_co:
        ncA,ncG=Dataset(AVGS[i]),Dataset(Grd_npth)
        for n in [0,10,11]:
            X,Z,tmp1=jr.get_section(ncG,ncA,'temp',[140],x_rng,tindx=n)
            # _,_,tmp2=jr.get_section(ncG,ncA,'rho',[60],[-70,-34])
            VAR2.append(tmp1)
        tmp_JFD=np.mean(np.array(VAR2),0)
        VAR2=[]
        ncA.close(); ncG.close()
        VAR1.append(tmp_JFD);
    VAR1=np.array(VAR1)

    Label_size=14
    t=0; data_lim=[-2.1,15]
    for i,j in zip(VAR1,My_date):
        t+=1
        Title_name='Years: '+j+f' (+{t-1:02d}~{t:02d})'
        # Figures
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        axs[0].set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
        im0=axs[0].contour(X,Z,i,colors='k',levels=[1,3])
        im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(X,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-Tcline,0)
        axs[0].set_xlim(-66.5,-42.5)
        im3=axs[1].contour(X,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
        im3.collections[1].set_linestyle('dashed')
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(X,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-5000,-Tcline)
        axs[1].set_xlim(-66.5,-42.5)
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if fig_bool:
            plt.savefig(save_pth+'ppt/'+save_name+'_'+j.replace('-','_'),
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(save_pth+save_name+'_'+j.replace('-','_'),bbox_inches='tight')
        plt.show()

# zonal averaged vertical section
def zonal_data_drift(data_nm,cmap,data_lim,**kargs):
    global Avg_pth,Grd_npth, save_pth,fig_bool
    
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams["font.weight"] = "regular"

    save_name='Zonal_temp_average_section'
    
    t_rng=[kargs['st'],kargs['ed']]
    
    # Read Grd
    TOPO=xr.open_dataset(Grd_npth).h.mean(dim='xi_rho')
    
    
    AVGS=np.sort([Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')])
    NC=xr.open_mfdataset(AVGS)
    DATA=NC[data_nm].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]
    zeta=NC['zeta'].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]

    zonal_m=DATA.mean(dim='xi_rho')
    zonal_zeta_m=zeta.mean(dim='xi_rho')
    
    lat=NC.lat_rho.values[:,0]
    
    if kargs['mean']=='ann':
        data=yearly_mean(zonal_m).rename({'year':'ocean_time'})
        d_zeta=yearly_mean(zonal_zeta_m).rename({'year':'ocean_time'})
    elif kargs['mean']=='season':
        data=season_mean(zonal_m).rename({'season':'ocean_time'})
        d_zeta=season_mean(zonal_zeta_m).rename({'season':'ocean_time'})
    elif kargs['mean']=='monthly':
        data=zonal_m.resample(ocean_time='1MS').mean()
        d_zeta=zonal_zeta_m.resample(ocean_time='1MS').mean()
    elif kargs['mean']=='monthly_clm':
        data=zonal_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        d_zeta=zonal_zeta_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
    elif kargs['mean']=='total':
        data=zonal_m.mean(dim='ocean_time',keepdims=True)
        d_zeta=zonal_zeta_m.mean(dim='ocean_time',keepdims=True)

    for i in d_zeta.values:
        Z=jr.zlevs(NC['Vtransform'].values[0], NC['Vstretching'].values[0],NC['theta_s'].values[0],\
               NC['theta_b'].values[0], NC['Tcline'].values[0], NC.s_rho.shape[0],1, TOPO.values, i)
    
    Label_size=12
    xtick_location = np.linspace(lat[0], lat[-1],6)
    xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]
    

    lat_m,z_m=np.meshgrid(lat,Z[:,0])

    for i in data:
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)
        
        s_name_S='Zonal_mean_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')

        # Figures
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        # fig.subplots_adjust(wspace=0, hspace=0)
        axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
        im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[1,3])
        im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(lat_m,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-NC['Tcline'].values[0],0)
        axs[0].set_xlim(-80,-23.5)
        im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
        im3.collections[1].set_linestyle('dashed')
        axs[0].set_xticks(ticks=xtick_location)
        axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(lat_m,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-5000,-NC['Tcline'].values[0])
        axs[1].set_xlim(-80,-23.5)
        axs[1].set_xticks(ticks=xtick_location)
        axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if fig_bool:
            plt.savefig(save_pth+'ppt/'+s_name_S,
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(save_pth+s_name_S,bbox_inches='tight')
        plt.show()

# zonal data diff soda
def zonal_data_diff_Soda(data_nm,cmap,cmap1,data_lim,data_lim1,**kargs):
    global Avg_pth,Grd_npth, save_pth,fig_bool, SODA_pth
    
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams["font.weight"] = "regular"


    t_rng=[kargs['st'],kargs['ed']]

    
    if data_nm=='zeta':
        data_soda_nm='ssh'
    elif data_nm=='u_eastward':
        data_soda_nm='u'
        # data_lim1=[-0.3,0.3]
    elif data_nm=='v_northward':
        data_soda_nm='v'
        # data_lim1=[-0.1,0.1]
    elif data_nm=='temp':
        data_soda_nm='temp'
        # data_lim1=[-3.5,3.5]
    elif data_nm=='salt':
        data_soda_nm='salt'
        # data_lim1=[-1,1.]


    # Proceses SODA
    SODA=xr.open_mfdataset(SODA_pth+'*.nc')[data_soda_nm].\
        loc[dict(yt_ocean=slice(-80,-23.5),time=slice(t_rng[0],t_rng[-1]))].rename({'time':'ocean_time'})
    zonal_soda_m=SODA.mean(dim='xt_ocean')
    SODA_lat=SODA.yt_ocean.values
    SODA_Z=SODA.st_ocean.values
    SODA_lat_m,SODA_Z_m=np.meshgrid(SODA_lat,SODA_Z)
    
    save_name='Zonal_temp_average_section'
    
    # Read Grd
    TOPO=xr.open_dataset(Grd_npth).h.mean(dim='xi_rho')
    
    
    AVGS=np.sort([Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')])
    NC=xr.open_mfdataset(AVGS)
    DATA=NC[data_nm].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]
    zeta=NC['zeta'].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]

    zonal_m=DATA.mean(dim='xi_rho')
    zonal_zeta_m=zeta.mean(dim='xi_rho')
    
    lat=NC.lat_rho.values[:,0]
    
    if kargs['mean']=='ann':
        data=yearly_mean(zonal_m).rename({'year':'ocean_time'})
        data_soda=yearly_mean(zonal_soda_m).rename({'year':'ocean_time'})
        d_zeta=yearly_mean(zonal_zeta_m).rename({'year':'ocean_time'})
    elif kargs['mean']=='season':
        data=season_mean(zonal_m).rename({'season':'ocean_time'})
        data_soda=season_mean(zonal_soda_m).rename({'season':'ocean_time'})
        d_zeta=season_mean(zonal_zeta_m).rename({'season':'ocean_time'})
    elif kargs['mean']=='monthly':
        data=zonal_m.resample(ocean_time='1MS').mean()
        data_soda=zonal_soda_m.resample(ocean_time='1MS').mean()
        d_zeta=zonal_zeta_m.resample(ocean_time='1MS').mean()
    elif kargs['mean']=='monthly_clm':
        data=zonal_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        data_soda=zonal_soda_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        d_zeta=zonal_zeta_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
    elif kargs['mean']=='total':
        data=zonal_m.mean(dim='ocean_time',keepdims=True)
        data_soda=zonal_soda_m.mean(dim='ocean_time',keepdims=True)
        d_zeta=zonal_zeta_m.mean(dim='ocean_time',keepdims=True)

    for i in d_zeta.values:
        Z=jr.zlevs(NC['Vtransform'].values[0], NC['Vstretching'].values[0],NC['theta_s'].values[0],\
               NC['theta_b'].values[0], NC['Tcline'].values[0], NC.s_rho.shape[0],1, TOPO.values, i)
    
    Label_size=12
    xtick_location = np.linspace(lat[0], lat[-1],6)
    xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]
    

    lat_m,z_m=np.meshgrid(lat,Z[:,0])
    
    
    # Interpolates SODA grid to model grid
    from scipy.interpolate import interp2d, griddata

 
    for i,j in zip(data,data_soda):
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)
        
        
        tmp_soda_=griddata( (SODA_lat_m.flatten(),-SODA_Z_m.flatten()),j.values.flatten(),
                (lat_m.flatten(),Z.flatten() ),
            method='linear',fill_value=np.nan)
        
        soda_re = tmp_soda_.reshape(lat_m.shape)

        model_soda=i-soda_re
        
        s_name_S1='Zonal_mean_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+'Model_SODA'+'_'+kargs['ed'].replace('-','')
        s_name_S='Zonal_mean_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+'SODA'+'_'+kargs['ed'].replace('-','')

        # Figure 1) Model - SODA  
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        # fig.subplots_adjust(wspace=0, hspace=0)
        axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
        # im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[1,3])
        # im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(lat_m,Z,model_soda,cmap=cmap1,vmin=data_lim1[0],vmax=data_lim1[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-NC['Tcline'].values[0],0)
        axs[0].set_xlim(-80,-23.5)
        # im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
        # im3.collections[1].set_linestyle('dashed')
        axs[0].set_xticks(ticks=xtick_location)
        axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(lat_m,Z,model_soda,cmap=cmap1,vmin=data_lim1[0],vmax=data_lim1[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-5000,-NC['Tcline'].values[0])
        axs[1].set_xlim(-80,-23.5)
        axs[1].set_xticks(ticks=xtick_location)
        axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if fig_bool:
            plt.savefig(save_pth+'ppt/'+s_name_S,
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(save_pth+s_name_S,bbox_inches='tight')
        plt.show()
        
        # Figure 2) SODA  
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        # fig.subplots_adjust(wspace=0, hspace=0)
        axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
        # im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[1,3])
        # im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(SODA_lat_m,-SODA_Z_m,j,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-NC['Tcline'].values[0],0)
        axs[0].set_xlim(-80,-23.5)
        # im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[1,3])
        # im3.collections[1].set_linestyle('dashed')
        axs[0].set_xticks(ticks=xtick_location)
        axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(SODA_lat_m,-SODA_Z_m,j,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-6000,-NC['Tcline'].values[0])
        axs[1].set_xlim(-80,-23.5)
        axs[1].set_xticks(ticks=xtick_location)
        axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if fig_bool:
            plt.savefig(save_pth+'ppt/'+s_name_S,
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(save_pth+s_name_S,bbox_inches='tight')
        plt.show()      

# Surface data soda diff
def Surface_data_Soda_diff(data_nm,lat_rng,cmap,My_levels,My_levels1,data_lim,data_lim1,**kargs):
    global Avg_pth,Grd_npth, save_pth,fig_bool, SODA_pth
    
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams["font.weight"] = "regular"

    # t_rng=[kargs['st'],kargs['ed']]
    
    if data_nm=='zeta':
        data_soda_nm='ssh'
    elif data_nm=='u_eastward':
        data_soda_nm='u'
    elif data_nm=='v_northward':
        data_soda_nm='v'
    elif data_nm=='temp':
        data_soda_nm='temp'
    elif data_nm=='salt':
        data_soda_nm='salt'
        
    
    print('!!! Code must be chaged for expansions !!!')
    print('!!! lon_rho,lat_rho=np.meshgrid(lon_rho[0,lat_co],lat_rho[lat_co,0]) !!!')

    # lon_rho,lat_rho=np.meshgrid(lon_rho[0,lat_co],lat_rho[lat_co,0])
        
    SODA_tmp=[SODA_pth+i for i in os.listdir(SODA_pth) if i.endswith('.nc')][0]
    AVGS=[Avg_pth+i for i in os.listdir(Avg_pth) if i.endswith('.nc')]
    Sample=xr.open_dataset(AVGS[0]); Sample_soda=xr.open_dataset(SODA_tmp)
        
    Sample_Data=Dataset(AVGS[0])
    lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
    lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])

    
    if [i for i in Sample[data_nm].coords].count('s_rho'): 
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(s_rho=Sample.s_rho.values[-1],\
                                                        eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]
        SODA_=xr.open_mfdataset(SODA_pth+'*.nc')[data_soda_nm].\
            loc[dict(st_ocean=Sample_soda.st_ocean.values[0],\
                     yt_ocean=slice(lat_rng[0],lat_rng[-1]),time=slice(kargs['st'],kargs['ed']))]\
                .rename({'time':'ocean_time'})
    else:
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,\
                                                        ocean_time=slice(kargs['st'],kargs['ed']))]
        SODA_=xr.open_mfdataset(SODA_pth+'*.nc')[data_soda_nm].\
            loc[dict(time=slice(kargs['st'],kargs['ed']))].rename({'time':'ocean_time'})

    if kargs['mean']=='ann':
        data=yearly_mean(data_).rename({'year':'ocean_time'})
        SODA=yearly_mean(SODA_).rename({'year':'ocean_time'})
    elif kargs['mean']=='season':
        data=season_mean(data_).rename({'season':'ocean_time'})
        SODA=season_mean(SODA_).rename({'season':'ocean_time'})
    elif kargs['mean']=='monthly':
        data=data_.resample(ocean_time='1MS').mean()
        SODA=SODA_.resample(ocean_time='1MS').mean()
    elif kargs['mean']=='monthly_clm':
        data=data_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        SODA=SODA_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
    elif kargs['mean']=='total':
        data=data_.mean(dim='ocean_time',keepdims=True)
        SODA=SODA_.mean(dim='ocean_time',keepdims=True)


    lat_soda,lon_soda=SODA.yt_ocean.values,SODA.xt_ocean.values
    lon_soda_m, lat_soda_m=np.meshgrid(lon_soda,lat_soda)

    for i,j in zip(data,SODA):
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)

        tmp_soda_=griddata( (lon_soda_m.flatten(),lat_soda_m.flatten()),j.values.flatten(),
                (lon_rho.flatten(),lat_rho.flatten() ),
            method='linear',fill_value=np.nan)
        
        soda_re = tmp_soda_.reshape(lon_rho.shape)
    
        model_soda=i-soda_re
    
        s_name_S1='Spherical_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+'Model_SODA'+'_'+kargs['ed'].replace('-','')
        s_name_S='Spherical_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                kargs['st'].replace('-','')+'_'+'SODA'+'_'+kargs['ed'].replace('-','')
    
        s_name_M=s_name_S.replace('Spherical','Merc')
        s_name_M1=s_name_S1.replace('Spherical','Merc')

        # figure SODA
        mm.Plot_SO_Spherical2(lon_soda_m,lat_soda_m,j.values,\
                              t_name,My_levels,cmap,data_lim,save_pth,s_name_S,fig_bool)
        mm.Plot_SO_Merc2(lon_soda_m,lat_soda_m,j,t_name,My_levels,cmap,data_lim,\
                          save_pth,s_name_M,fig_bool)
        # figure SODA
        mm.Plot_SO_Spherical2(lon_rho,lat_rho,model_soda.values,\
                              t_name,My_levels1,cmap,data_lim1,save_pth,s_name_S1,fig_bool)
        mm.Plot_SO_Merc2(lon_rho,lat_rho,model_soda,t_name,My_levels1,cmap,data_lim1,\
                          save_pth,s_name_M1,fig_bool)
            











