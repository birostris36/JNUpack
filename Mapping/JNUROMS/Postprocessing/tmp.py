def zonal_data_diff_Soda(data_nm,cmap,data_lim,**kargs):
    global Avg_pth,Grd_npth, save_pth,fig_bool, SODA_pth
    
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.rcParams["font.weight"] = "regular"


    # Proceses SODA
    SODA=xr.open_mfdataset(SODA_pth+'*.nc')[data_nm].loc[dict(time=slice(t_rng[0],t_rng[-1]))].rename({'time':'ocean_time'})
    zonal_soda_m=SODA.mean(dim='xt_ocean')
    SODA_lat=SODA.yt_ocean.values
    SODA_Z=SODA.st_ocean.values
    SODA_lat_m,SODA_Z_m=np.meshgrid(SODA_lat,SODA_Z)
    

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

    
    for i in data_soda:
        tmp_i=i.values
        tmp_U_=griddata( (SODA_lat_m.flatten(),SODA_Z_m.flatten()),i.flatten(),
                    (lat_m.flatten(),z_m.flatten() ),
                method='linear',fill_value=np.nan)
    

    for i,j in zip(data,data_soda):
        if kargs['mean']=='monthly':
            t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
        else:
            t_name=str(i.ocean_time.values)
        
        tmp_soda_=griddata( (SODA_lat_m.flatten(),SODA_Z_m.flatten()),i.values.flatten(),
                (lat_m.flatten(),z_m.flatten() ),
            method='linear',fill_value=np.nan)
        
        soda_re = tmp_soda_.reshape(lat_m.shape)

        
        
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