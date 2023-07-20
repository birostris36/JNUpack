# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:06:13 2023

@author: shjo9
"""



import cartopy.crs as ccrs
import cartopy.feature as cf
import datetime as dt
import cmocean
import sys 
sys.path.append('D:/JNUpack/')
sys.path.append('D:/JNUpack/JNUROMS/')
import Tools.JNUROMS as jr
import Tools.Inputs as ti
# from Mapping.Tools import d_modules as mm
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
# from pptx import Presentation # 라이브러리 
# from pptx.util import Inches,Cm, Pt # 사진, 표등을 그리기 위해
# from pptx.enum.text import PP_ALIGN
import cartopy
import matplotlib.path as mpath
from copy import deepcopy


def Plot_SO_Spherical2(lonA,latA,MyDATA,t_name,My_levels,CMAP,Mylim,w_path,save_name,fig_bool=False):

    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                       subplot_kw={'projection': Spheric})
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.rotate_labels=False
    gl.xlabels_top,gl.ylabels_right = True,True
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
     
    M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    #M=plt.pcolormesh(lonA, latA, MyDATA,
    #              transform=PC,cmap=CMAP,vmin=Mylim[0],vmax=Mylim[-1])
    plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
#    ax.set_extent([0, 360, -80, -24], crs=PC)
    
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()
def Plot_SO_Spherical3(lonA,latA,MyDATA,t_name,CMAP,my_lim,w_path,save_name,fig_bool=False):

    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,11),
                       subplot_kw={'projection': Spheric})
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':32,'fontweight':'regular'})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.rotate_labels=False
    gl.xlabels_top,gl.ylabels_right = True,True
    gl.xlabel_style = gl.ylabel_style = {"size" : 26}
     
    #M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    M=plt.pcolormesh(lonA, latA, MyDATA,
                  transform=PC,cmap=CMAP,vmin=my_lim[0],vmax=my_lim[-1])
    #plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
    # ax.set_extent([0, 360, -80, -24], crs=PC)
    
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=1., axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.08,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()
    
def Plot_SO_Merc2(lonA,latA,MyDATA,t_name,My_levels,CMAP,Mylim,w_path,save_name,fig_bool=False):
    
    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,4),
                       subplot_kw={'projection': MERC})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 24}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular'})

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
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()
     
def Plot_SO_Merc3(lonA,latA,MyDATA,t_name,CMAP,Mylim,w_path,save_name,fig_bool=False):
    
    Spheric=ccrs.SouthPolarStereo(central_longitude=0.0,globe=None)
    PC = ccrs.PlateCarree(central_longitude=0.0,globe=None)
    MERC=ccrs.Mercator(central_longitude=180.0,globe=None)
    
    # Now we will create axes object having specific projection 

    fig, ax = plt.subplots(1, 1, figsize=(12.5,4),
                       subplot_kw={'projection': MERC})

    gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                      linewidth=.6, color='k', alpha=0.45, linestyle='-.',\
                          )
    gl.xlabels_top,gl.ylabels_right = False,False
    gl.xlabel_style = gl.ylabel_style = {"size" : 24}
    
    # To plot borders and coastlines, we can use cartopy feature
    ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
    ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
    ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular'})

    # M=plt.contourf(lonA,latA,MyDATA,cmap=CMAP,levels=My_levels,transform=PC)
    M=plt.pcolormesh(lonA, latA, MyDATA,
                  transform=PC,cmap=CMAP)
    plt.clim(Mylim[0],Mylim[-1])
    
    # crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
 #   ax.set_extent([0, 360, -80, -24], crs=PC)
    ax.tick_params(axis='both', which='major', labelsize=28)

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=.1, axes_class=plt.Axes)

    fig.add_axes(ax_cb)
    cb=plt.colorbar(M,extend='both',pad=0.01,cax=ax_cb)
    cb.set_label(label='', weight='regular',fontsize=28)
    cb.ax.tick_params(labelsize=19)

    plt.tight_layout()
    if fig_bool:
        plt.savefig(w_path+'/ppt/'+save_name,
                facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(w_path+'/'+save_name)
    plt.show()

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

class Diag():
    
    def __init__(self,Grd_npth,Avg_pth,Log_npth,sv_pth,SODA_pth):
        
        if sv_pth[-1]!='/':
            sv_pth=sv_pth+'/'
        
        try:
            sv_pth=sv_pth+Log_npth.split('/')[-1].split('.')[0].replace('Log_','')+'/'
        except:
            sv_pth=Log_npth.split('/')[-1].split('.')[0]+'/'
        
        try:
            os.mkdir(sv_pth)
            os.mkdir(sv_pth+'ppt/')
        except:
            pass
        
        self.grd=Grd_npth
        self.avg=Avg_pth
        self.sv=sv_pth
        self.log=Log_npth
        self.soda=SODA_pth
        
        # Define figparams
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

    def check_inputs(self): #Log_npth,Avg_pth,save_pth
    
        DATA=[self.avg+i for i in os.listdir(self.avg)]
        Parameters=ti.Log_Manager(self.log).get_parameters(NN=12)
        SM=ti.Sample_Manager(DATA[0])
        CPP=SM.get_cpp_Sample(8)
        NLM_LBC=SM.get_LBC_Sample()
        Files=SM.get_inputs_files()
        
        prs = Presentation() # 파워포인트 객체 선언
        
        ### !!! Info !!! ###
        bullet_slide_layout = prs.slide_layouts[1] # 1 : 제목 + 내용 슬라이드
        slide = prs.slides.add_slide(bullet_slide_layout) # 기존에 있던 슬라이드에 추가
        
        # 제목
        title_shape = slide.placeholders[0] 
        title_shape.text = 'Adding a Bullet Slide'
        
        # 내용
        body_shape = slide.placeholders[1]
        tf = body_shape.text_frame
        tf.text = 'Find the bullet slide layout'
        
        # 단락 추가
        p = tf.add_paragraph()
        p.text = 'Use _TextFrame.text for first bullet'
        p.level = 1  # 1 : 들여쓰기 레벨
        
        # 단락 추가 
        p = tf.add_paragraph()
        p.text = 'Use _TextFrame.add_paragraph() for subsequent bullets'
        p.level = 2  # 2 : 들여쓰기 레벨
        
        ### !!! Figure 01 !!! ###
        # img_path = 'D:/Working_hub/OneDrive/Projects/ROMS/Figs/'
        
        blank_slide_layout = prs.slide_layouts[6] # 6 : 제목/내용이 없는 '빈' 슬라이드
        slide = prs.slides.add_slide(blank_slide_layout)
        
        left = top = Inches(2.5)
        width = height = Inches(7)
        # width, hegith가 없을 경우 원본 사이즈로 
        # pic = slide.shapes.add_picture(img_path+'Model_momentum.png', left, top, width=width,height=height)
        
        left = Inches(3)
        width = Inches(6)
        height = Inches(3)
        # pic = slide.shapes.add_picture(img_path+'AVISO_adt.png', left, top, width=width,height=height)
        
        
        ### !!! CPP Table !!! ###
        
        title_only_slide_layout = prs.slide_layouts[5]
        slide = prs.slides.add_slide(title_only_slide_layout)
        shapes = slide.shapes
        
        title_shape = slide.placeholders[0] 
        title_shape.text = 'Activated CPP options'
        title_shape.text_frame.paragraphs[0].font.bold = True
        
        tmp_L,tmp_M = CPP.shape
        tmp = CPP.values
        
        rows = tmp_L
        cols = tmp_M
        left = top = Inches(0*2.54)
        width = Inches(3.5*2.54)
        height = Inches(1.5*2.54)
        
        Shape = shapes.add_table(rows, cols, left, top, width, height)
        tbl =  Shape._element.graphic.graphicData.tbl
        style_id = '{616DA210-FB5B-4158-B5E0-FEB733F419BA}'
        tbl[0][-1].text = style_id    
        table = Shape.table
        # set column widths
        # table.columns.width = Inches(3.5*2.54)
        # table.columns.height = Inches(1.5*2.54)
        
        
        for i in range(tmp_L):
            for j in range(tmp_M):
                table.cell(i,j).text=tmp[i][j]
                table.cell(i,j).text_frame.paragraphs[0].font.size = Pt(12)
                table.cell(i,j).text_frame.paragraphs[0].font.bold = False
                table.cell(i,j).text_frame.paragraphs[0].alignment= PP_ALIGN.CENTER
        
        ### !!! LBC Table !!! ###
        title_only_slide_layout2 = prs.slide_layouts[5]
        slide2 = prs.slides.add_slide(title_only_slide_layout2)
        shapes2 = slide2.shapes
        
        title_shape2 = slide2.placeholders[0] 
        title_shape2.text = 'NLB_LBC'
        title_shape2.text_frame.paragraphs[0].font.bold = True
        
        tmp_L,tmp_M = NLM_LBC.values.shape
        tmp = NLM_LBC.values
        
        rows,cols = tmp_L,tmp_M  
        left = top = Inches(0*2.54)
        width = Inches(1.5*2.54)
        height = Inches(1.*2.54)
        
        Shape2 = shapes2.add_table(rows, cols, left, top, width, height)
        tbl2 =  Shape2._element.graphic.graphicData.tbl
        style_id2 = '{616DA210-FB5B-4158-B5E0-FEB733F419BA}'
        tbl2[0][-1].text = style_id    
        table2 = Shape2.table
        # set column widths
        # table.columns.width = Inches(3.5*2.54)
        # table.columns.height = Inches(1.5*2.54)
        
        for i in range(tmp_L):
            for j in range(tmp_M):
                table2.cell(i,j).text=tmp[i][j]
                table2.cell(i,j).text_frame.paragraphs[0].font.size = Pt(12)
                table2.cell(i,j).text_frame.paragraphs[0].alignment= PP_ALIGN.CENTER
                if j==0:
                    table2.cell(i,j).text_frame.paragraphs[0].font.bold = True
        
        
        ### !!! Parameters !!! ###
        title_only_slide_layout3 = prs.slide_layouts[5]
        slide3 = prs.slides.add_slide(title_only_slide_layout3)
        shapes3 = slide3.shapes
        
        title_shape3 = slide3.placeholders[0] 
        title_shape3.text = 'Parameters'
        title_shape3.text_frame.paragraphs[0].font.bold = True
        
        tmp_L,tmp_M = np.shape(Parameters)
        
        rows,cols = tmp_L,tmp_M  
        left = top = Inches(0*2.54)
        width = Inches(3.*2.54)
        height = Inches(1.*2.54)
        
        Shape3 = shapes3.add_table(rows, cols, left, top, width, height)
        tbl3 =  Shape3._element.graphic.graphicData.tbl
        style_id2 = '{616DA210-FB5B-4158-B5E0-FEB733F419BA}'
        tbl3[0][-1].text = style_id    
        table3 = Shape3.table
        # set column widths
        # table.columns.width = Inches(3.5*2.54)
        # table.columns.height = Inches(1.5*2.54)
        
        for i in range(tmp_L):
            for j in range(tmp_M):
                table3.cell(i,j).text=Parameters[i][j]
                table3.cell(i,j).text_frame.paragraphs[0].font.size = Pt(11)
                table3.cell(i,j).text_frame.paragraphs[0].alignment= PP_ALIGN.CENTER
                table3.cell(i,j).text_frame.paragraphs[0].font.bold = False
                if not j%2:
                    table3.cell(i,j).text_frame.paragraphs[0].font.bold = True
        
        
        ### !!! Input files !!! ###
        title_only_slide_layout4 = prs.slide_layouts[5]
        slide4 = prs.slides.add_slide(title_only_slide_layout4)
        shapes4 = slide4.shapes
        
        title_shape4 = slide4.placeholders[0] 
        title_shape4.text = 'Input files'
        title_shape4.text_frame.paragraphs[0].font.bold = True
        
        tmp_L,tmp_M = np.shape(Files)
        
        rows,cols = tmp_L,tmp_M  
        left = top = Inches(0*2.54)
        width = Inches(3.*2.54)
        height = Inches(1.*2.54)
        
        Shape4 = shapes4.add_table(rows, cols, left, top, width, height)
        tbl4 =  Shape4._element.graphic.graphicData.tbl
        style_id2 = '{616DA210-FB5B-4158-B5E0-FEB733F419BA}'
        tbl4[0][-1].text = style_id    
        table4 = Shape4.table
        # set column widths
        # table.columns.width = Inches(3.5*2.54)
        # table.columns.height = Inches(1.5*2.54)
        
        for i in range(tmp_L):
            for j in range(tmp_M):
                table4.cell(i,j).text=Files[i][j]
                table4.cell(i,j).text_frame.paragraphs[0].font.size = Pt(11)
                table4.cell(i,j).text_frame.paragraphs[0].alignment= PP_ALIGN.CENTER
                table4.cell(i,j).text_frame.paragraphs[0].font.bold = False
                if not j%2:
                    table4.cell(i,j).text_frame.paragraphs[0].font.bold = True
        
        ######################################
        print(self.sv+self.log.split('/')[-1].split('.')[0]+'.pptx')
        prs.save(self.sv+self.log.split('/')[-1].split('.')[0]+'.pptx')

    # Evaluates model stability
    def Stability01(self):
        
        with open(self.log) as f:
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
                if i>stid and len(tmp1) and tmp1[0].isnumeric() and not i>len(Model_Log)-1000:
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
        
        # PD['NET_VOLUME'].plot()
        # PD.columns
        
        RX=rx.split('\n')
        RX1=RX[0][:16]
        RX2=RX[-1][:-7]
        
        Title_name='Topo ('+RX1+' / '+RX2+')'
        Model_Times1 = PD['YYYY-MM-DD']
        Model_Times2 = [i[2:4] for i in Model_Times1]
        Label_size = 18
        
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
        axs[0].legend(fontsize=12)
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
        axs[1].legend(fontsize=12)
        
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
        axs[2].legend(fontsize=12)
        
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
        axs[3].legend(fontsize=12)
    
        plt.tight_layout()
        if True:
            plt.savefig(self.sv+'ppt/'+'Model_momentum_logs',
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.sv+'Model_momentum_logs')
        plt.show()
        
    def Stability02(self):
    
        ncs=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
    
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
        
        #AICE_fp1=np.polyfit(t[200:],AICE_timeseries[200:],1)
        #AICE_trend=np.polyval(AICE_fp1,t[200:])
        
        # 
        Ross_lat_rng,Ross_lon_rng=[-73,-60], [180,215]
        Wedd_lat_rng,Wedd_lon_rng=[-73,-60], [300,360]
        
        Sample_Data=Dataset(ncs[0])
        lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
        Ross_lat_co=np.where((lat_rho[:,0]>=Ross_lat_rng[0])&(lat_rho[:,0]<=Ross_lat_rng[-1]))[0]
        Ross_lon_co=np.where((lon_rho[0,:]>=Ross_lon_rng[0])&(lon_rho[0,:]<=Ross_lon_rng[-1]))[0]
        
        Wedd_lat_co=np.where((lat_rho[:,0]>=Wedd_lat_rng[0])&(lat_rho[:,0]<=Wedd_lat_rng[-1]))[0]
        Wedd_lon_co=np.where((lon_rho[0,:]>=Wedd_lon_rng[0])&(lon_rho[0,:]<=Wedd_lon_rng[-1]))[0]
        
        ROSS=xr.open_mfdataset(self.avg+'*.nc').zeta.loc[dict(xi_rho=Ross_lon_co,eta_rho=Ross_lat_co)].mean(dim=['eta_rho','xi_rho']).values
        WEDD=xr.open_mfdataset(self.avg+'*.nc').zeta.loc[dict(xi_rho=Wedd_lon_co,eta_rho=Wedd_lat_co)].mean(dim=['eta_rho','xi_rho']).values
    
        #SST_fp1=np.polyfit(t,SST_timeseries,1)
        #SST_trend=np.polyval(SST_fp1,t)
        
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
        
        f1 = axs[2].plot(Model_Times1,ROSS, label='ROSS zeta',color='k',linewidth=2,zorder=0)
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
        # axs[2].set_ylim(-100,100)
        
        f1 = axs[3].plot(Model_Times1,WEDD, label='WEDD zeta',color='k',linewidth=2,zorder=0)
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
        # axs[3].set_ylim(7,12.6)
        plt.tight_layout()
        if True:
            plt.savefig(self.sv+'ppt/'+'Model_stability2',
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.sv+'Model_stability2')
        plt.show()

    def Stability03(self):
        # lon_rng,lat_rng --> For transport (Sv)        
       
        AVGS=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
    
        Sample_Data=Dataset(AVGS[0])
        lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
        
        def Cal_Sv(lat_rng,lon_rng):
            dl_x=np.min(np.diff(lon_rho))
            lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
            lon_co=np.where((lon_rho[0,:]>lon_rng-dl_x/2)&(lon_rho[0,:]<lon_rng+dl_x/2))[0]
            if len(lon_co) !=1  :
                print('!!!!!!!!!!!!!!!!!!')
                raise
            # lon_rho,lat_rho=np.meshgrid(lon_rho[0,lon_co],lat_rho[lat_co,0])
            NC=xr.open_mfdataset(AVGS[0])
            # Transport of Drake passage (Sv)
            U=xr.open_mfdataset(AVGS)['u_eastward'].loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze()
            ZETA=xr.open_mfdataset(AVGS)['zeta'].loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze().values
            TOPO=xr.open_dataset(self.grd).h.loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze().values
            PN=xr.open_dataset(self.grd).pn.loc[dict(xi_rho=lon_co,eta_rho=lat_co)].squeeze().values
            T,th,at=U.shape
            X_dist=np.tile(1/PN,th).reshape([th,at])
            Z=jr.zlevs(NC['Vtransform'].values, NC['Vstretching'].values,NC['theta_s'].values,\
                   NC['theta_b'].values, NC['Tcline'].values,  U.s_rho.shape[0] ,5, TOPO, ZETA[-1])
            # Z_=jr.zlevs(NC['Vtransform'].values, NC['Vstretching'].values,NC['theta_s'].values,\
            #        NC['theta_b'].values, NC['Tcline'].values,  U.s_rho.shape[0] ,1, TOPO, ZETA[-1])
            # lat_m,_=np.meshgrid(lat_rho,Z_[:,0])
            delta_Z=Z[1:,:]-Z[:-1,:]
            CELL=delta_Z*X_dist #(Unit: m*m=m**2)
            CELL=CELL*10**-6
            SV=[]
            for i in U.values:
                SV_sum=np.nansum(CELL*i)    
                SV.append(SV_sum)
            SV=np.array(SV)
            return SV
        
        lat_DRK_rng,lon_DRK_rng=[-65,-55],294
        lat_SAF_rng,lon_SAF_rng=[-80,-24],24.5
        lat_SAU_rng,lon_SAU_rng=[-80,-24],129

        DRK_Sv=Cal_Sv(lat_DRK_rng,lon_DRK_rng)
        SAF_Sv=Cal_Sv(lat_SAF_rng,lon_SAF_rng)
        SAU_Sv=Cal_Sv(lat_SAU_rng,lon_SAU_rng)

        Model_Times1 = range(len(SAU_Sv))
        Model_Times2=[str(i/12)[:2] for i in Model_Times1]
    
        Label_size = 25
        fig, axs = plt.subplots(3,1,figsize=(11,6.5),constrained_layout = True,
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.,1]},dpi=200)
        f1 = axs[0].plot(Model_Times1,DRK_Sv, label='Drake passage (66W)',color='k',linewidth=2,zorder=0)
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
        axs[0].set_ylim(90,170)
        #! Fig2 
        f1 = axs[1].plot(Model_Times1,SAF_Sv, label='South of Africa (24.5E)',color='k',linewidth=2,zorder=0)
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
        axs[1].set_ylim(90,170)        
        
        f1 = axs[2].plot(Model_Times1,SAU_Sv, label='South of Australia (129E)',color='k',linewidth=2,zorder=0)
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
        axs[2].legend(fontsize=18,loc='upper right')
        axs[2].set_ylim(90,170)        
        plt.tight_layout()
        if True:
            plt.savefig(self.sv+'ppt/'+'Major_transports',
                    facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            plt.savefig(self.sv+'Major_transports')
        plt.show()

    # Plot Surface results
    def data_drift(self,data_nm,lat_rng,My_levels,cmap,data_lim,**kargs):
    
        AVGS=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
        
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
    
        # if kargs['mean']=='ann':
        #     data=yearly_mean(data_).rename({'year':'ocean_time'})
        # elif kargs['mean']=='season':
        #     data=season_mean(data_).rename({'season':'ocean_time'})
        # elif kargs['mean']=='monthly':
        #     data=data_.resample(ocean_time='1MS').mean()
        # elif kargs['mean']=='monthly_clm':
        #     data=data_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        # elif kargs['mean']=='total':
        #     data=data_.mean(dim='ocean_time',keepdims=True)
            
    # =============================================================================
        if kargs['mean']=='ann':
            data=data_.groupby('ocean_time.year').mean().rename({'year':'ocean_time'})
        elif kargs['mean']=='season':
            data=data_.groupby('ocean_time.season').mean().rename({'season':'ocean_time'})
        elif kargs['mean']=='monthly':
            data=data_.resample(ocean_time='1MS').mean()
        elif kargs['mean']=='monthly_clm':
            data=data_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        elif kargs['mean']=='total':
            data=data_.mean(dim='ocean_time',keepdims=True)
    # =============================================================================
        try:
            os.mkdir(self.sv+'Surface_mean_'+kargs['mean']+'_'+data_nm)
            os.mkdir(self.sv+'Surface_mean_'+kargs['mean']+'_'+data_nm+'/ppt')
        except:
            pass
                
        for i in data:
            # if kargs['mean']=='monthly':
    
            #     t_name=pd.to_datetime(i.ocean_time.values).strftime('%Y-%m')
            # else:
            #     t_name=str(i.ocean_time.values)
            t_name=str(i.ocean_time.values)[:7]
                
                
            s_name_S='Spherical_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                    kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
            
            s_name_M=s_name_S.replace('Spherical','Merc')
            
            Plot_SO_Spherical2(data.lon_rho.values,data.lat_rho.values,i.squeeze().values,\
                                  t_name,My_levels,cmap,data_lim,\
                                  self.sv+'Surface_mean_'+kargs['mean']+'_'+data_nm\
                                  ,s_name_S,True)
                
    #        Plot_SO_Merc2(lon_rho,lat_rho,i,t_name,My_levels,cmap,data_lim,\
    #                          save_pth+'Surface_mean_'+kargs['mean']+'_'+data_nm\
    #                          ,s_name_M,fig_bool)
                
    # Auger Temp vertical section
    def Auger_temp_section(self,data_nm,My_levels,cmap,has_year_zero=False,**kargs):
        
        plt.rcParams['contour.negative_linestyle'] = 'solid'
        plt.rcParams["font.weight"] = "regular"
    
        save_name='Auger_temp'
        
        AVGS=np.sort([self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')])
        
        # Process Times
        t_rng=[kargs['st'],kargs['ed']]
        OGCM_TIMES=MFDataset(self.avg+'*nc')['ocean_time']
        units=OGCM_TIMES.units
    
    
        if has_year_zero:
            OGCM_times=num2date(OGCM_TIMES[:],units,has_year_zero=has_year_zero)
            Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
            Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
            OGCM_times= np.array([dt.datetime(int(i.strftime().split(' ')[0].split('-')[0]),\
                         int(i.strftime().split(' ')[0].split('-')[1]),\
                             int(i.strftime().split(' ')[0].split('-')[-1])) for i in OGCM_times])
            TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
    
        else:
            OGCM_times=num2date(OGCM_TIMES[:],units)
            Tst=dt.datetime(int(t_rng[0].split('-')[0]), int(t_rng[0].split('-')[1]),1)
            Ted=dt.datetime(int(t_rng[1].split('-')[0]), int(t_rng[1].split('-')[1]),31)
            TIMES_co=np.where( (OGCM_times>=Tst)&(OGCM_times<=Ted) )[0]
    
        Avg_co=list(set(TIMES_co//12))
    
        My_date=list(np.sort(list(set([i.strftime('%Y') for i in OGCM_times[TIMES_co]]))))
        
        # Get My Grid info
        ncG=Dataset(self.grd)
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
        
        try:
            os.mkdir(self.sv+save_name) 
            os.mkdir(self.sv+save_name+'/ppt') 
        except:
            pass
        
        # =============================================================================
        # Process Times
        # =============================================================================
        x_rng=[-66.5,-42.5]
        VAR1,VAR2=[],[]
        for i in Avg_co:
            ncA,ncG=Dataset(AVGS[i]),Dataset(self.grd)
            for n in [0,10,11]:
                X,Z,tmp1=jr.get_section(ncG,ncA,'temp',[140],x_rng,tindx=n)
                # _,_,tmp2=jr.get_section(ncG,ncA,'rho',[60],[-70,-34])
                VAR2.append(tmp1)
            tmp_JFD=np.mean(np.array(VAR2),0)
            VAR2=[]
            ncA.close(); ncG.close()
            VAR1.append(tmp_JFD);
        VAR1=np.array(VAR1)
    
        Label_size=18
        t=0; data_lim=[-2.1,15]
        for i,j in zip(VAR1,My_date):
            t+=1
            Title_name='Years: '+j+f' (+{t-1:02d}~{t:02d})'
            # Figures
            fig, axs = plt.subplots(2,1,figsize=(10.5,6),
                                    sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
            axs[0].set_title(Title_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
            im0=axs[0].contour(X,Z,i,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
            axs[0].clabel(im0, inline=1, fontsize=14)
    
            # im0.collections[1].set_linestyle('dashed')
            im1=axs[0].pcolor(X,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
            axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
            axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
            axs[0].set_ylim(-Tcline,0)
            axs[0].set_xlim(-66.5,-42.5)
            im3=axs[1].contour(X,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
            axs[1].clabel(im0, inline=1, fontsize=14)
            # im3.collections[1].set_linestyle('dashed')
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
            if True:
                plt.savefig(self.sv+save_name+'/ppt/'+save_name+'_'+j.replace('-','_'),
                            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
                plt.savefig(self.sv+save_name+'/'+save_name+'_'+j.replace('-','_'),bbox_inches='tight')
            plt.show()

    # zonal averaged vertical section
    def zonal_data_drift(self,data_nm,cmap,data_lim,**kargs):
        
        plt.rcParams['contour.negative_linestyle'] = 'solid'
        plt.rcParams["font.weight"] = "regular"
    
        save_name='Zonal_temp_average_section'
        
        t_rng=[kargs['st'],kargs['ed']]
        
        # Read Grd
        TOPO=xr.open_dataset(self.grd).h.mean(dim='xi_rho')
        
        AVGS=np.sort([self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')])
        NC=xr.open_mfdataset(AVGS)
        DATA=NC[data_nm].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]
        zeta=NC['zeta'].loc[dict(ocean_time=slice(t_rng[0],t_rng[-1]))]
    
        zonal_m=DATA.mean(dim='xi_rho')
        zonal_zeta_m=zeta.mean(dim='xi_rho')
        
        lat=NC.lat_rho.values[:,0]
        
        # if kargs['mean']=='ann':
        #     data=yearly_mean(zonal_m).rename({'year':'ocean_time'})
        #     d_zeta=yearly_mean(zonal_zeta_m).rename({'year':'ocean_time'})
        # elif kargs['mean']=='season':
        #     data=season_mean(zonal_m).rename({'season':'ocean_time'})
        #     d_zeta=season_mean(zonal_zeta_m).rename({'season':'ocean_time'})
        # elif kargs['mean']=='monthly':
        #     data=zonal_m.resample(ocean_time='1MS').mean()
        #     d_zeta=zonal_zeta_m.resample(ocean_time='1MS').mean()
        # elif kargs['mean']=='monthly_clm':
        #     data=zonal_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        #     d_zeta=zonal_zeta_m.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        # elif kargs['mean']=='total':
        #     data=zonal_m.mean(dim='ocean_time',keepdims=True)
        #     d_zeta=zonal_zeta_m.mean(dim='ocean_time',keepdims=True)
        if kargs['mean']=='ann':
            data=zonal_m.groupby('ocean_time.year').mean().rename({'year':'ocean_time'})
            d_zeta=zonal_zeta_m.groupby('ocean_time.year').mean().rename({'year':'ocean_time'})
        elif kargs['mean']=='season':
            data=zonal_m.groupby('ocean_time.season').mean().rename({'season':'ocean_time'})
            d_zeta=zonal_zeta_m.groupby('ocean_time.season').mean().rename({'season':'ocean_time'})
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
        NC.close()
        Label_size=12
        xtick_location = np.linspace(lat[0], lat[-1],6)
        xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]
        
        lat_m,z_m=np.meshgrid(lat,Z[:,0])
        
        try:
            os.mkdir(self.sv+'Zonal_mean_'+kargs['mean']+'_'+data_nm) 
            os.mkdir(self.sv+'Zonal_mean_'+kargs['mean']+'_'+data_nm+'/ppt') 
        except:
            pass
    
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
          #  im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
          #  axs[0].clabel(im0, inline=1, fontsize=14)
          #  im0.collections[1].set_linestyle('dashed')
            im1=axs[0].pcolor(lat_m,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
            axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
            axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
            axs[0].set_ylim(-NC['Tcline'].values[0],0)
            axs[0].set_xlim(-80,-23.5)
          #  im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
        #    axs[1].clabel(im0, inline=1, fontsize=14)
            axs[0].set_xticks(ticks=xtick_location)
            axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
            axs[0].set_facecolor(color='#dddddd')
            
            # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
            im2=axs[1].pcolor(lat_m,Z,i,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
            axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
            axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
            axs[1].set_ylim(-5000,-NC['Tcline'].values[0])
            axs[1].set_xlim(-80,-23.5)
            axs[1].set_xticks(ticks=xtick_location)
            axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
            axs[1].set_facecolor(color='#dddddd')
            divider = make_axes_locatable(axs[1])
            cax = divider.append_axes("bottom", size="7%", pad=.35)
            cax.tick_params(labelsize=Label_size)
            cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
            h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
            if True:
                plt.savefig(self.sv+'Zonal_mean_'+kargs['mean']+'_'+data_nm+'/ppt/'+s_name_S,
                            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
                plt.savefig(self.sv+'Zonal_mean_'+kargs['mean']+'_'+data_nm+'/'+s_name_S,bbox_inches='tight')
            plt.show()
    
    # zonal data diff soda
    def zonal_data_diff_Soda(self,data_nm,cmap,cmap1,data_lim,data_lim1,**kargs):
        
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
        if (data_nm=='u_eastward') or (data_nm=='v_northward') :
            SODA=xr.open_mfdataset(self.soda+'*.nc')[data_soda_nm].\
                loc[dict(yu_ocean=slice(-80,-23.5),time=slice(t_rng[0],t_rng[-1]))].rename({'time':'ocean_time'})
            zonal_soda_m=SODA.mean(dim='xu_ocean')
            SODA_lat=SODA.yu_ocean.values
        else:
            SODA=xr.open_mfdataset(self.soda+'*.nc')[data_soda_nm].\
                loc[dict(yt_ocean=slice(-80,-23.5),time=slice(t_rng[0],t_rng[-1]))].rename({'time':'ocean_time'})
            zonal_soda_m=SODA.mean(dim='xt_ocean')
            SODA_lat=SODA.yt_ocean.values
            
        SODA_Z=SODA.st_ocean.values
        SODA_lat_m,SODA_Z_m=np.meshgrid(SODA_lat,SODA_Z)
        
        save_name='Zonal_temp_average_section'
        
        # Read Grd
        TOPO=xr.open_dataset(self.grd).h.mean(dim='xi_rho')
        
        
        AVGS=np.sort([self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')])
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
        
        
        try:
            [os.mkdir(self.sv+'Zonal_meanDiff_'+kargs['mean']+'_'+data_nm+'_'+i) for i in ['MODEL_SODA','SODA']]
            [os.mkdir(self.sv+'Zonal_meanDiff_'+kargs['mean']+'_'+data_nm+'_'+i+'/ppt') for i in ['MODEL_SODA','SODA']]
        except:
            pass
    
     
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
            
            s_name_S1='Zonal_meanDiff_'+data_nm+'_'+t_name.replace('-','')+'_'+\
                    kargs['st'].replace('-','')+'_'+'Model_SODA'+'_'+kargs['ed'].replace('-','')
            s_name_S='Zonal_meanDiff_'+data_nm+'_'+t_name.replace('-','')+'_'+\
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
            axs[0].set_facecolor(color='#dddddd')
            
            # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
            im2=axs[1].pcolor(lat_m,Z,model_soda,cmap=cmap1,vmin=data_lim1[0],vmax=data_lim1[-1])
            axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
            axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
            axs[1].set_ylim(-5000,-NC['Tcline'].values[0])
            axs[1].set_xlim(-80,-23.5)
            axs[1].set_xticks(ticks=xtick_location)
            axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
            axs[1].set_facecolor(color='#dddddd')
            divider = make_axes_locatable(axs[1])
            cax = divider.append_axes("bottom", size="7%", pad=.35)
            cax.tick_params(labelsize=Label_size)
            cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
            h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
            if True:
                plt.savefig(self.sv+'Zonal_meanDiff_'+kargs['mean']+'_'+data_nm+'_'+'MODEL_SODA/'+'ppt/'+s_name_S1,
                            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
                plt.savefig(self.sv+'Zonal_meanDiff_'+kargs['mean']+'_'+data_nm+'_'+'MODEL_SODA/'+s_name_S1,bbox_inches='tight')
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
            axs[0].set_facecolor(color='#dddddd')
            
            # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
            im2=axs[1].pcolor(SODA_lat_m,-SODA_Z_m,j,cmap=cmap,vmin=data_lim[0],vmax=data_lim[-1])
            axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
            axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
            axs[1].set_ylim(-6000,-NC['Tcline'].values[0])
            axs[1].set_xlim(-80,-23.5)
            axs[1].set_xticks(ticks=xtick_location)
            axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
            axs[1].set_facecolor(color='#dddddd')
            divider = make_axes_locatable(axs[1])
            cax = divider.append_axes("bottom", size="7%", pad=.35)
            cax.tick_params(labelsize=Label_size)
            cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
            h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
            if True:
                plt.savefig(self.sv+'Zonal_meanDiff_'+kargs['mean']+'_'+data_nm+'_'+'SODA/'+'ppt/'+s_name_S,
                            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
                plt.savefig(self.sv+'Zonal_meanDiff_'+kargs['mean']+'_'+data_nm+'_'+'SODA/'+s_name_S,bbox_inches='tight')
            plt.show()      

    # Surface data soda diff
    def Surface_data_Soda_diff(self,data_nm,lat_rng,cmap,cmap1,My_levels,My_levels1,data_lim,data_lim1,**kargs):
        
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
            
        SODA_tmp=[self.soda+i for i in os.listdir(self.soda) if i.endswith('.nc')][0]
        AVGS=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
        Sample=xr.open_dataset(AVGS[0]); Sample_soda=xr.open_dataset(SODA_tmp)
            
        Sample_Data=Dataset(AVGS[0])
        lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
        lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
        lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
    
        if [i for i in Sample[data_nm].coords].count('s_rho'): 
            data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(s_rho=Sample.s_rho.values[-1],\
                                                            eta_rho=lat_co,\
                                                            ocean_time=slice(kargs['st'],kargs['ed']))]
                
            if (data_nm=='u_eastward') or (data_nm=='v_northward') :
                SODA_=xr.open_mfdataset(self.soda+'*.nc')[data_soda_nm].\
                    loc[dict(st_ocean=Sample_soda.st_ocean.values[0],\
                             yu_ocean=slice(lat_rng[0],lat_rng[-1]),time=slice(kargs['st'],kargs['ed']))]\
                        .rename({'time':'ocean_time'})
            else:
                SODA_=xr.open_mfdataset(self.soda+'*.nc')[data_soda_nm].\
                    loc[dict(st_ocean=Sample_soda.st_ocean.values[0],\
                             yt_ocean=slice(lat_rng[0],lat_rng[-1]),time=slice(kargs['st'],kargs['ed']))]\
                        .rename({'time':'ocean_time'})
        else:
            data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,\
                                                            ocean_time=slice(kargs['st'],kargs['ed']))]
            SODA_=xr.open_mfdataset(self.soda+'*.nc')[data_soda_nm].\
                loc[dict(yt_ocean=slice(lat_rng[0],lat_rng[-1]),time=slice(kargs['st'],kargs['ed']))].rename({'time':'ocean_time'})
            
        # if kargs['mean']=='ann':
        #     data=yearly_mean(data_).rename({'year':'ocean_time'})
        #     SODA=yearly_mean(SODA_).rename({'year':'ocean_time'})
        # elif kargs['mean']=='season':
        #     data=season_mean(data_).rename({'season':'ocean_time'})
        #     SODA=season_mean(SODA_).rename({'season':'ocean_time'})
        # elif kargs['mean']=='monthly':
        #     data=data_.resample(ocean_time='1MS').mean()
        #     SODA=SODA_.resample(ocean_time='1MS').mean()
        # elif kargs['mean']=='monthly_clm':
        #     data=data_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        #     SODA=SODA_.groupby('ocean_time.month').mean().rename({'month':'ocean_time'})
        # elif kargs['mean']=='total':
        #     data=data_.mean(dim='ocean_time',keepdims=True)
        #     SODA=SODA_.mean(dim='ocean_time',keepdims=True)
            
        # =============================================================================
        if kargs['mean']=='ann':
            data=data_.groupby('ocean_time.year').mean()\
                .rename({'year':'ocean_time'})
            SODA=SODA_.groupby('ocean_time.year').mean()\
                .rename({'year':'ocean_time'})
    
        elif kargs['mean']=='season':
            data=data_.groupby('ocean_time.season').mean()\
                .rename({'season':'ocean_time'})
            SODA=SODA_.groupby('ocean_time.season').mean()\
                .rename({'season':'ocean_time'})
            
        elif kargs['mean']=='monthly':
            data=data_.resample(ocean_time='1MS').mean()
            SODA=SODA_.resample(ocean_time='1MS').mean()
        elif kargs['mean']=='monthly_clm':
            data=data_.groupby('ocean_time.month').mean()\
                .rename({'month':'ocean_time'})
            SODA=SODA_.groupby('ocean_time.month').mean()\
                .rename({'month':'ocean_time'})
        elif kargs['mean']=='total':
            data=data_.mean(dim='ocean_time',keepdims=True)
            SODA=SODA_.mean(dim='ocean_time',keepdims=True)
        # =============================================================================
        try:
            [os.mkdir(self.sv+'Surface_Diff_'+kargs['mean']+'_'+data_nm+'_'+i)\
             for i in ['MODEL_SODA','SODA']]
            [os.mkdir(self.sv+'Surface_Diff_'+kargs['mean']+'_'+data_nm+'_'+i+
                      '/ppt') for i in ['MODEL_SODA','SODA']]
        except:
            pass
        
        try:
            lat_soda,lon_soda=SODA.yt_ocean.values,SODA.xt_ocean.values
        except:
            lat_soda,lon_soda=SODA.yu_ocean.values,SODA.xu_ocean.values
    
        
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
            Plot_SO_Spherical2(lon_soda_m,lat_soda_m,j.values,\
                                  t_name,My_levels,cmap,data_lim,self.sv+\
                                  'Surface_Diff_'+kargs['mean']+'_'+data_nm+'_'+'SODA/',s_name_S,True)
            Plot_SO_Merc2(lon_soda_m,lat_soda_m,j,t_name,My_levels,cmap,data_lim,\
                              self.sv+'Surface_Diff_'+kargs['mean']+'_'+data_nm+'_'+'SODA/',s_name_M,True)
            # figure SODA
            Plot_SO_Spherical2(lon_rho,lat_rho,model_soda.values,\
                                  t_name,My_levels1,cmap1,data_lim1,self.sv+\
                                  'Surface_Diff_'+kargs['mean']+'_'+data_nm+'_'+'MODEL_SODA/',s_name_S1,True)
            Plot_SO_Merc2(lon_rho,lat_rho,model_soda,t_name,My_levels1,cmap1,data_lim1,\
                              self.sv+\
                              'Surface_Diff_'+kargs['mean']+'_'+data_nm+'_'+'MODEL_SODA/',s_name_M1,True)
                
    # Surface Linear trend
    def Surface_data_trend(self,data_nm,lat_rng,cmap,**kargs):
        
        AVGS=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
        
        Sample_Data=Dataset(AVGS[0])
        lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
        lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
        lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
            
        Sample=xr.open_dataset(AVGS[0])
        
        if [i for i in Sample[data_nm].coords].count('s_rho'):
            data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(s_rho=Sample.s_rho.values[-1],\
                                                            eta_rho=lat_co,\
                                                            ocean_time=slice(kargs['st'],kargs['ed']))]     
        else:
            data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,\
                                                            ocean_time=slice(kargs['st'],kargs['ed']))]
        # Linear trend
        data=data_.assign_coords({'TT':('ocean_time',range(len(data_.ocean_time)))})
        data=data.swap_dims({"ocean_time":"TT"})
        data_s=data.polyfit(dim='TT',deg=1,skipna=True)
        Coef=data_s.polyfit_coefficients[0]
        Coef_var=Coef.values*12 # (m/year)
        
        My_lim=np.nanmean(Coef_var)+np.nanstd(Coef_var)*1/2
        my_lim=[-My_lim,My_lim]
        
        t_name=data_nm.upper()+' trend : '+kargs['st'].replace('-','-')+' ~ '+kargs['ed'].replace('-','-')  
        t_name=''
        s_name_S='Surface_trend_'+data_nm+'_'+kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
    
        try:
            os.mkdir(self.sv+s_name_S)
            os.mkdir(self.sv+s_name_S+'/ppt')
        except:
            pass
        
        Plot_SO_Spherical3(data.lon_rho,data.lat_rho,Coef_var,\
                              t_name,cmap,my_lim,\
                              self.sv+s_name_S,s_name_S,True)
            
        # Plot_SO_Merc3(lon_rho,lat_rho,i,t_name,My_levels,cmap,data_lim,\
        #                   save_pth+'Surface_mean_'+kargs['mean']+'_'+data_nm\
        #                   ,s_name_M,fig_bool)

    def zonal_data_trend(self,data_nm,lat_rng,cmap,**kargs):
    
        AVGS=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
          
        TOPO=xr.open_dataset(self.grd).h.mean(dim='xi_rho').values
    
        Sample_Data=Dataset(AVGS[0])
        lat_rho,lon_rho=Sample_Data['lat_rho'][:],Sample_Data['lon_rho'][:]
        lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
        lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
            
        Sample=xr.open_dataset(AVGS[0])
        
        NC=xr.open_mfdataset(AVGS[0])
        data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,ocean_time=slice(kargs['st'],kargs['ed']))].mean(dim='xi_rho')
        
        data=data_.assign_coords({'TT':('ocean_time',range(len(data_.ocean_time)))})
        data=data.swap_dims({"ocean_time":"TT"})
        data_s=data.polyfit(dim='TT',deg=1,skipna=True)
        Coef=data_s.polyfit_coefficients[0]
        Coef_var=Coef.values*12 # (m/year)
        
        My_lim=np.nanmean(Coef_var)+np.nanstd(Coef_var)*3
        my_lim=[-My_lim,My_lim]
        
        Z=jr.zlevs(NC['Vtransform'].values, NC['Vstretching'].values,NC['theta_s'].values,\
               NC['theta_b'].values, NC['Tcline'].values,  data_.s_rho.shape[0] ,1, TOPO, np.zeros_like(TOPO))
    
        lat=NC.lat_rho.values[:,0]
        lat_m,z_m=np.meshgrid(lat,Z[:,0])
    
        t_name=data_nm.upper()+' trend : '+kargs['st'].replace('-','-')+' ~ '+kargs['ed'].replace('-','-')   
        s_name_S='Zonal_trend_'+data_nm+'_'+kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
    
        try:
            os.mkdir(self.sv+s_name_S)
            os.mkdir(self.sv+s_name_S+'/ppt')
        except:
            pass
    
        Label_size=12
        xtick_location = np.linspace(lat[0], lat[-1],6)
        xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]
            
        fig, axs = plt.subplots(2,1,figsize=(6,4),
                                sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
        # ax = plt.gca()
        # fig.subplots_adjust(wspace=0, hspace=0)
        axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
      #  im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
      #  axs[0].clabel(im0, inline=1, fontsize=14)
      #  im0.collections[1].set_linestyle('dashed')
        im1=axs[0].pcolor(lat_m,Z,Coef_var,cmap=cmap,vmin=my_lim[0],vmax=my_lim[-1])
        axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
        axs[0].set_ylim(-NC['Tcline'].values,0)
        axs[0].set_xlim(-80,-23.5)
      #  im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
    #    axs[1].clabel(im0, inline=1, fontsize=14)
        axs[0].set_xticks(ticks=xtick_location)
        axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        axs[0].set_facecolor(color='#dddddd')
        
        # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
        im2=axs[1].pcolor(lat_m,Z,Coef_var,cmap=cmap,vmin=my_lim[0],vmax=my_lim[-1])
        axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
        axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
        axs[1].set_ylim(-5000,-NC['Tcline'].values)
        axs[1].set_xlim(-80,-23.5)
        axs[1].set_xticks(ticks=xtick_location)
        axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
        axs[1].set_facecolor(color='#dddddd')
        
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("bottom", size="7%", pad=.35)
        cax.tick_params(labelsize=Label_size)
        cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
        h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
        if True:
            plt.savefig(self.sv+s_name_S+'/ppt/'+s_name_S,
                        facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
            print(self.sv+s_name_S+'/'+s_name_S)
            plt.savefig(self.sv+s_name_S+'/'+s_name_S,bbox_inches='tight')
        plt.show()

    # def Zonal_data_trend(self,data_nm,lat_rng,cmap,**kargs):
        
    #     TOPO=xr.open_dataset(self.grd).h.mean(dim='xi_rho')
    
    #     AVGS=[self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
        
    #     NC=Dataset(AVGS[0])
    #     lat_rho,lon_rho=NC['lat_rho'][:],NC['lon_rho'][:]
    #     lat_co=np.where((lat_rho[:,0]>=lat_rng[0])&(lat_rho[:,0]<=lat_rng[-1]))[0]
    #     lon_rho,lat_rho=np.meshgrid(lon_rho[0,:],lat_rho[lat_co,0])
                
    #     data_=xr.open_mfdataset(AVGS)[data_nm].loc[dict(eta_rho=lat_co,ocean_time=slice(kargs['st'],kargs['ed']))]\
    #         .mean(dim='xi_rho').squeeze()
    
    #     zeta=xr.open_mfdataset(AVGS[0])['zeta'].mean(dim='xi_rho').squeeze()
    #     lat=NC['lat_rho'][:,0]
    
    #     # Linear trend
    #     data=data_.assign_coords({'TT':('ocean_time',range(len(data_.ocean_time)))})
    #     data=data.swap_dims({"ocean_time":"TT"})
    #     data_s=data.polyfit(dim='TT',deg=1,skipna=True)
    #     Coef=data_s.polyfit_coefficients[0]
    #     Coef_var=Coef.values*12 # (m/year)
        
    #     My_lim=np.nanmean(Coef_var)+np.nanstd(Coef_var)*3.5
    #     my_lim=[-My_lim,My_lim]
        
    #     # zonal_zeta_m=zeta.mean(dim='xi_rho')
    
    #     zeros_zeta=np.zeros_like(zeta[0,:])
    #     Z=jr.zlevs(NC['Vtransform'][:], NC['Vstretching'][:],NC['theta_s'][:],\
    #            NC['theta_b'][:], NC['Tcline'][:], NC['s_rho'][:].shape[0],1, TOPO.values, zeros_zeta)
            
    #     lat_m,z_m=np.meshgrid(lat,Z[:,0])
    
    #     Label_size=12
    #     xtick_location = np.linspace(lat[0], lat[-1],6)
    #     xtick_labels = [f'{ii:0.1f}' for ii in xtick_location]    
    
    #     #CMAP_trend=plt.get_cmap('seismic',15)
        
    
    #     t_name=data_nm.upper()+' trend ['+kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')+']'
            
    #     s_name_S='Zonal_trend_'+data_nm+'_'+\
    #         kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
            
    #     try:
    #         os.mkdir(self.sv+s_name_S) 
    #         os.mkdir(self.sv+s_name_S+'/ppt') 
    #     except:
    #         pass
            
    #     # Figures
    #     fig, axs = plt.subplots(2,1,figsize=(6,4),
    #                             sharex=True,gridspec_kw={'height_ratios': [1, 1.3],'wspace':0, 'hspace':0.05},dpi=200)
    #     # fig.subplots_adjust(wspace=0, hspace=0)
    #     axs[0].set_title(t_name,loc='right',fontdict={'fontsize':Label_size,'fontweight':'regular'})
    #   #  im0=axs[0].contour(lat_m,Z,i,colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
    #   #  axs[0].clabel(im0, inline=1, fontsize=14)
    #   #  im0.collections[1].set_linestyle('dashed')
    #     im1=axs[0].pcolor(lat_m,Z,Coef_var,cmap=cmap,vmin=my_lim[0],vmax=my_lim[-1])
    #     axs[0].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    #     axs[0].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size, color='k',right=True)
    #     axs[0].set_ylim(-NC['Tcline'][0],0)
    #     axs[0].set_xlim(-80,-23.5)
    #   #  im3=axs[1].contour(lat_m,Z,i,vmin=data_lim[0],vmax=data_lim[-1],colors='k',levels=[-1.5,1.5,4.5,8,11],linestyle='-')
    # #    axs[1].clabel(im0, inline=1, fontsize=14)
    #     axs[0].set_xticks(ticks=xtick_location)
    #     axs[0].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
    #     axs[0].set_facecolor(color='#dddddd')
        
    #     # im4=axs[1].clabel(colors='k',CS=im3,inline=True,fmt='%1.f')
    #     im2=axs[1].pcolor(lat_m,Z,Coef_var,cmap=cmap,vmin=my_lim[0],vmax=my_lim[-1])
    #     axs[1].tick_params(axis='x', direction='in', length=4.5, pad=8, labelsize=Label_size, labelcolor='k', top=True)
    #     axs[1].tick_params(axis='y', direction='in', length=4.5, pad=8, labelsize=Label_size,  color='k',right=True)
    #     axs[1].set_ylim(-5000,-NC['Tcline'][-1])
    #     axs[1].set_xlim(-80,-23.5)
    #     axs[1].set_xticks(ticks=xtick_location)
    #     axs[1].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=.7)
    #     axs[1].set_facecolor(color='#dddddd')
    #     divider = make_axes_locatable(axs[1])
    #     cax = divider.append_axes("bottom", size="7%", pad=.35)
    #     cax.tick_params(labelsize=Label_size)
    #     cax.set_ylabel('',{'fontsize':Label_size,'fontweight':'bold','style':'italic'})
    #     h = fig.colorbar(im1, ax=axs[:],label='',cax=cax,orientation="horizontal",extend='both',aspect=50)
    #     if True:
    #         plt.savefig(self.sv+s_name_S+'/ppt/'+s_name_S,
    #                     facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    #         plt.savefig(self.sv+s_name_S+'/'+s_name_S,bbox_inches='tight')
    #     plt.show()

    def mk_OHC_nc(self,sv_dt=False):
        
        import gsw.density as gsw_d
        import gsw.conversions as gsw_c
        from gsw._wrapped_ufuncs import cp_t_exact
        
        NCs = [self.avg+i for i in os.listdir(self.avg) if i.endswith('.nc')]
        
        NC=Dataset(NCs[0])
        TOPO=Dataset(self.grd)['h'][:]
        tmp_zeta=np.zeros_like(TOPO)
        
        Z=jr.zlevs(NC['Vtransform'][:], NC['Vstretching'][:],NC['theta_s'][:],\
               NC['theta_b'][:], NC['Tcline'][:], NC['s_rho'].shape[0],1, TOPO[:], tmp_zeta)
            
        Z_w=jr.zlevs(NC['Vtransform'][:], NC['Vstretching'][:],NC['theta_s'][:],\
               NC['theta_b'][:], NC['Tcline'][:], NC['s_rho'].shape[0],5, TOPO[:], tmp_zeta)
               
        TEMP,SALT=np.flip(xr.open_mfdataset(NCs).temp,axis=1), np.flip(xr.open_mfdataset(NCs).salt,axis=1)
        DEPTH=-np.flip(Z,axis=0)
        DEPTH_w=-np.flip(Z_w,axis=0)
        dZ=DEPTH_w[1:,:,:]-DEPTH_w[:-1,:,:]
        
        fac_2000m=deepcopy(DEPTH)
        fac_2000m[fac_2000m>2000]=0
        fac_2000m[fac_2000m<=2000]=1

        # Potential temp --> Conservative temp
        CT=gsw_c.CT_from_pt(SALT,TEMP) #CT = gsw_CT_from_pt(SA,pt)
        
        rho = gsw_d.rho(SALT,CT,DEPTH)
        
        # Potential temperature --> In-situ temperature
        t=gsw_c.t_from_CT(SALT,CT,DEPTH)
        
        # Heat capacity
        CP = cp_t_exact(SALT,t,DEPTH)
        
        OHC_=CP*CT*rho*dZ
        OHC_=OHC_*fac_2000m

        # Integrates from ref depth (2000m) 
        OHC=OHC_.sum(dim='s_rho',skipna=False)
        OHC=OHC.rename('OHC')
    
        # d_t,d_at,d_on=OHC.shape
        
        ### !!!!!!!!! Only for squre grid !!!!!!!!!!!!! ###
        lon_coord,lat_coord = OHC.lon_rho.values[0,:], OHC.lat_rho.values[:,0] 
        OHC=OHC.assign_coords({'longitude':('xi_rho',lon_coord )})
        OHC=OHC.assign_coords({'latitude':('eta_rho',lat_coord )})
        OHC=OHC.swap_dims({"eta_rho":"latitude"}).swap_dims({"xi_rho":"longitude"})
        
        if sv_dt:
            try:
                OHC.to_netcdf(self.sv+self.log.split('/')[-1].split('.')[0].replace('Log','OHC')+'.nc')
            except:
                OHC.to_netcdf(self.sv+'OHC_'+self.log.split('/')[-1].split('.')[0]+'.nc')
        return OHC


    def OHC_trend(self,OCH,lat_rng, lon_rng,cmap,my_lim,**kargs):
        
        OCH=OCH.loc[dict(longitude=slice(lon_rng[0],lon_rng[-1]),latitude=slice(lat_rng[0],lat_rng[-1]),\
                ocean_time=slice(kargs['st'],kargs['ed']))]
        
        OCH=OCH.assign_coords({'TT':('ocean_time',range(len(OCH.ocean_time)))})
        OCH=OCH.swap_dims({"ocean_time":"TT"})
        OCH_s=OCH.polyfit(dim='TT',deg=1,skipna=True)
        Coef=OCH_s.polyfit_coefficients[0]
        Coef_var=Coef.values*12 # (m/year)
            
        t_name='OHC trend : '+kargs['st'].replace('-','-')+' ~ '+kargs['ed'].replace('-','-')  
        t_name=''
        s_name_S='OHC_trend_'+kargs['st'].replace('-','')+'_'+kargs['ed'].replace('-','')
        
        try:
            os.mkdir(self.sv+s_name_S)
            os.mkdir(self.sv+s_name_S+'/ppt')
        except:
            pass
        
        Plot_SO_Spherical3(OCH.lon_rho,OCH.lat_rho,Coef_var,\
                              t_name,cmap,my_lim,\
                              self.sv+s_name_S,s_name_S,True)
        Plot_SO_Merc3(OCH.lon_rho,OCH.lat_rho,Coef_var,\
                              t_name,cmap,my_lim,\
                              self.sv+s_name_S,s_name_S,True)
            












