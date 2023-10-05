import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import os
import sys
import numpy as np

plt.rcParams["font.family"] = 'Arial'

def dta_colr(myNm):
    # print(myNm)
    if myNm == 'ERA5':
        myClr='C0'
    elif myNm=='ERSSTv5':
        myClr='C1'
    elif myNm=='HadlSST':
        myClr='C2'
    elif myNm=='OISST':
        myClr='C3'
    elif myNm=='ECCO':
        myClr='C4' 
    elif myNm=='EN4':
        myClr='C5'
    elif myNm=='GECCO':
        myClr='C6'
    elif myNm=='IAP':
        myClr='C7'
    elif myNm=='ISHII':
        myClr='C8'
    else:
        raise
    return myClr
# mySetting={
#     'figsize': (10,3.7),
#     'mylabel': [],
#     'Label_size':18,
#     'title_loc':'right',
#     # 'xdomain_label_freq':[5::12*4],
#     'fontParams':'Arial'
# }

def recall_myCMAP(Cname):
    if Cname=='myblc2':
        Mycolorlist=['#1b2c62','#1f4181','#2455a1','#3877ba',\
            '#529bd2','#71b8e4','#91d2f2','#b2e0fa','#cbebf8',\
                '#e3f4fb','#f2f9e3','#fcf0b4','#fddb81','#fdc152',\
                    '#fca12f','#f8822b','#ef5e29','#e03b28','#cc1e27',\
                        '#ae191f']     
    else: 
        pass
    Mycmap = LinearSegmentedColormap.from_list('',Mycolorlist,N=256)
    return Mycmap

def myClrbr(myCname,myLIM,N):
    myCMAP=recall_myCMAP(myCname)
    mylevels=np.linspace(myLIM[0], myLIM[-1], N+1,endpoint=True)
    CMAP = ListedColormap(myCMAP(
    np.linspace(0, 1, len(mylevels)-1,endpoint=True)) )
    return CMAP,mylevels

class figmaster:
    
    def __init__(self,mySetting):
        self.figsize=mySetting['figsize']
        self.mylabel=mySetting['mylabel']
        self.Label_size=mySetting['Label_size']
        # self.xdomain_label=mySetting['xdomain_label_freq']
        self.wpth=mySetting['wpth']
        self.title_loc=mySetting['title_loc']

        
    def create_loc(wpth):
        try:
            os.mkdir(wpth)
        except:
            print('!!! Directory already exits !!!')
            ans=input('!!! Delete it ? (y/n) !!!')
            if ans=='y':
                os.rmdir(wpth)
                os.mkdir(wpth)
            else:
                pass
                
    def myCrtpy_cyl(self,LAT,LON,DATA,CMAP,LEVELS):
        PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
        PC180 = ccrs.PlateCarree(central_longitude=180.0,globe=None)
        fig, ax = plt.subplots(1, 1, figsize=(12.5,6),
                            subplot_kw={'projection': PC180},dpi=200)
        gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                            linewidth=.6, color='k', alpha=0.45, linestyle='-.')
        gl.xlabels_top,gl.ylabels_right = False,False
        gl.xlabel_style = gl.ylabel_style = {"size" : 20}
        ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
        ax.add_feature(cf.LAND,color=[.75,.75,.75],zorder=100)
        ax.set_title(self.ttl[0],loc=self.ttl[-1],fontdict={'fontsize':32,'fontweight':'regular'})

        M=plt.contourf(LON,LAT,DATA,cmap=CMAP,levels=LEVELS,transform=PC)

        ax.set_extent([LON[0][0], LON[0][-1], LAT[0][0], LAT[-1][0]], crs=PC)
        ax.tick_params(axis='both', which='major', labelsize=28)
        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
        cb.set_label(label='', weight='regular',fontsize=24)
        cb.ax.tick_params(labelsize=19)
        # cb.ax.yaxis.set_major_formatter(FormatStrFormatter(FMT))
        plt.tight_layout()
        if 1:
            # plt.savefig(w_path+'/ppt/'+s_name,
            #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(w_path+'/'+s_name,bbox_inches='tight')
        plt.show()
        
    
    # Now we will create axes object having specific projection 
    def myCrtpy_sph_pcolor(self,LAT,LON,DATA,CMAP,myName):
        Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
        PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
        fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                        subplot_kw={'projection': Spheric})
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
        ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
        ax.set_title(myName,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

        gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                        linewidth=.6, color='k', alpha=0.45, linestyle='-.')
        gl.rotate_labels=False
        gl.xlabels_top,gl.ylabels_right = True,True
        gl.xlabel_style = gl.ylabel_style = {"size" : 26}
        
        # M=plt.contourf(LON,LAT,DATA,cmap=CMAP,levels=LEVELS,transform=PC)
        # M=plt.contourf(LON,LAT,DATA,cmap=CMAP,transform=PC,vmin=-3.5,vmax=3.5)

        M=plt.pcolormesh(LON,LAT,DATA,cmap=CMAP,transform=PC,vmin=-3.5,vmax=3.5)

        ax.set_extent([LON[0][0], LON[0][-1], LAT[0][0], LAT[-1][0]], crs=PC)
        
        ax.tick_params(axis='both', which='major', labelsize=28)

        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
        cb.set_label(label='', weight='regular',fontsize=28)
        cb.ax.tick_params(labelsize=19)
        plt.tight_layout()
        if 1:
            myName.replace(' ','_')
            # plt.savefig(w_path+'/ppt/'+save_name,
            #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.wpth+'/'+myName)
        plt.show()
        
    def myCrtpy_sph2(self,LAT,LON,DATA,HATCH,CMAP,LEVELS,myName):
        Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
        PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
        fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                        subplot_kw={'projection': Spheric})
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
        ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
        ax.set_title(myName,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

        gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                        linewidth=.6, color='k', alpha=0.45, linestyle='-.')
        gl.rotate_labels=False
        gl.xlabels_top,gl.ylabels_right = True,True
        gl.xlabel_style = gl.ylabel_style = {"size" : 26}
        
        plt.contourf(LON,LAT,HATCH,levels=LEVELS,colors='none',hatches='.',transform=PC,zorder=2)
        M=plt.contourf(LON,LAT,DATA,cmap=CMAP,levels=LEVELS,transform=PC,zorder=0)
        
        ax.set_extent([LON[0][0], LON[0][-1], LAT[0][0], LAT[-1][0]], crs=PC)
        
        ax.tick_params(axis='both', which='major', labelsize=28)

        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
        cb.set_label(label='', weight='regular',fontsize=28)
        cb.ax.tick_params(labelsize=19)
        plt.tight_layout()
        if 1:
            myName.replace(' ','_')
            # plt.savefig(w_path+'/ppt/'+save_name,
            #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.wpth+'/'+myName.replace(' ','_'))
        plt.show()

    def myCrtpy_sph2_box(self,LAT,LON,DATA,HATCH,CMAP,LEVELS,myName,lat_rng,lon_rng):
        
        a,b=[lon_rng[0], lon_rng[0]],[lat_rng[0],lat_rng[-1]]
        c,d=[lon_rng[-1], lon_rng[-1]],[lat_rng[0],lat_rng[-1]]
        e,f=[lon_rng[0], lon_rng[-1]],[lat_rng[0],lat_rng[0]]
        g,h=[lon_rng[0], lon_rng[-1]],[lat_rng[-1],lat_rng[-1]]
        
        Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
        PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
        fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                        subplot_kw={'projection': Spheric})
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
        ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
        ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
        ax.set_title(myName,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

        gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                        linewidth=.6, color='k', alpha=0.45, linestyle='-.')
        gl.rotate_labels=False
        gl.xlabels_top,gl.ylabels_right = True,True
        gl.xlabel_style = gl.ylabel_style = {"size" : 26}
        
        ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
        ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
        ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
        ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5,zorder=200)
        
        plt.contourf(LON,LAT,HATCH,levels=LEVELS,colors='none',hatches='.',transform=PC,zorder=2)
        M=plt.contourf(LON,LAT,DATA,cmap=CMAP,levels=LEVELS,transform=PC,zorder=0)
        
        ax.set_extent([LON[0][0], LON[0][-1], LAT[0][0], LAT[-1][0]], crs=PC)
        
        ax.tick_params(axis='both', which='major', labelsize=28)

        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
        cb.set_label(label='', weight='regular',fontsize=28)
        cb.ax.tick_params(labelsize=19)
        plt.tight_layout()
        if 1:
            myName.replace(' ','_')
            # plt.savefig(w_path+'/ppt/'+save_name,
            #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.wpth+'/'+myName.replace(' ','_'))
        plt.show()

    def Zonal_mean_ver1(self,data1d,latR,myName,fig_bool=True):
        Label_size = 18
        fig, axs = plt.subplots(1,1,figsize=(5,6.5),constrained_layout = True,
                            dpi=200)
        f1 = axs.plot(data1d,latR, label='',color='k',linewidth=2,zorder=0)
        axs.set_title(myName,loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'italic'})
        axs.tick_params(axis='both', labelsize=Label_size)
        axs.grid(axis='x',linestyle='-.')
        # xtick_location = time[5::12*4]
        # xtick_labels = time2[5::12*4]
        # axs.set_xticks(ticks=xtick_location)
        # axs.set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
        axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
        axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
        plt.tight_layout()
        if fig_bool:
            myName.replace(' ','_')
            # plt.savefig(wnpth'/ppt/'+save_name,
            #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.wpth+'/'+myName.replace(' ','_'),bbox_inches='tight')
        plt.show()





def plot_pcs(time,time2,pc,t_name,w_path,save_name,fig_bool=True):
    Label_size = 18
    fig, axs = plt.subplots(1,1,figsize=(10,3.7),constrained_layout = True,
                        dpi=200)
    f1 = axs.plot(time,pc, label='KINETIC_ENRG',color='k',linewidth=2,zorder=0)
    axs.set_title(t_name,loc='right',fontdict={'fontsize':20,'fontweight':'regular','fontstyle':'italic'})
    axs.tick_params(axis='both', labelsize=Label_size)
    axs.grid(axis='x',linestyle='-.')
    xtick_location = time[5::12*4]
    xtick_labels = time2[5::12*4]
    axs.set_xticks(ticks=xtick_location)
    axs.set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    axs.tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs.tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
    plt.tight_layout()
    if fig_bool:
        # plt.savefig(wnpth'/ppt/'+save_name,
        #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(self.wpth+'/'+save_name,bbox_inches='tight')
    plt.show()


def plot_trend(self,y_est,y_err):
    ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2,color='r')
    
    
    
    
    
    
    
    


