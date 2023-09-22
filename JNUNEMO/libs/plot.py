import matplotlib as mpl
mpl.use('Agg') #Generates figures in Backend
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np

def plot_subs(time1,time2,pd_data,Title_name,wpth,wname,fig_bool=0):

    Label_size = 18
    lbls=pd_data.columns
    lbl_N=len(lbls)
    Figsize=(11,2.1*lbl_N)
    Height_ratios=list(np.ones(lbl_N))
    xtick_location = time1[::12*10]
    xtick_labels = time2[::12*10]
    
    fig, axs = plt.subplots(lbl_N,1,figsize=Figsize,constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': Height_ratios},dpi=200)

    for n,L in zip(range(lbl_N), lbls):
        axs[n].plot(time1,pd_data[L], label=L,color='k',linewidth=2,zorder=0)

        axs[n].tick_params(axis='y', labelsize=Label_size)     
        axs[n].set_xticks(ticks=xtick_location)
        axs[n].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        axs[n].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
        axs[n].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
        axs[n].legend(fontsize=12,loc='upper right')
        
    if fig_bool:
        plt.savefig(wpth+'ppt/'+wname,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(wpth+wname)
    plt.show()
    
def plot_2d(kwargs):
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    # MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                    subplot_kw={'projection': PC})
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cf.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(kwargs['Title_name'],loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                    linewidth=.6, color='k', alpha=0.45, linestyle='-.')
    gl.rotate_labels=False
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
    
    # M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    M=plt.pcolormesh(kwargs['lon'], kwargs['lat'], kwargs['data'],
                transform=PC,cmap=kwargs['CMAP'],vmin=kwargs['Mylim'][0],vmax=kwargs['Mylim'][-1])
    # ax.set_extent([0, 360, -80, -24], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)
    if kwargs['fig_bool']:
        plt.savefig(kwargs['wpth']+'/ppt/'+kwargs['wname'],
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(kwargs['wpth']+'/'+kwargs['wname'],bbox_inches='tight')
    plt.show()
    
def plot_2df(kwargs):
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    # MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                    subplot_kw={'projection': PC})
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cf.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(kwargs['Title_name'],loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                    linewidth=.6, color='k', alpha=0.45, linestyle='-.')
    gl.rotate_labels=False
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
    
    M=plt.contourf(kwargs['lon'], kwargs['lat'], kwargs['data'], cmap=kwargs['CMAP'],\
        levels=kwargs['levels'],transform=PC)
    #M=plt.pcolormesh(kwargs['lon'], kwargs['lat'], kwargs['data'],
    #            transform=PC,cmap=kwargs['CMAP'],vmin=kwargs['Mylim'][0],vmax=kwargs['Mylim'][-1])
    ax.set_extent([kwargs['lon'][0,0], kwargs['lon'][-1,-1],\
        kwargs['lat'][0,0], kwargs['lat'][-1,-1]], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)
    if kwargs['fig_bool']:
        plt.savefig(kwargs['wpth']+'/ppt/'+kwargs['wname'],
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(kwargs['wpth']+'/'+kwargs['wname'],bbox_inches='tight')
    plt.show()