# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:19:32 2023

@author: shjo9
"""
a,b=[180, 180],[-73,-60]
c,d=[215, 215],[-73,-60]
e,f=[180, 215],[-73,-73]
g,h=[180, 215],[-60,-60]

a1,b1=[300, 300],[-73,-60]
c1,d1=[360, 360],[-73,-60]
e1,f1=[300, 360],[-73,-73]
g1,h1=[300, 360],[-60,-60]

 # Now we will create axes object having specific projection 
fig, ax = plt.subplots(1, 1, figsize=(12.5,4),
                   subplot_kw={'projection': MERC})

gl = ax.gridlines(crs=PC, draw_labels=True,y_inline=False,x_inline=False,
                  linewidth=.6, color='k', alpha=0.45, linestyle='-.')
gl.xlabels_top,gl.ylabels_right = False,False
gl.xlabel_style = gl.ylabel_style = {"size" : 24}

# To plot borders and coastlines, we can use cartopy feature
ax.add_feature(cf.COASTLINE.with_scale("110m"), lw=1,zorder=110)
ax.add_feature(cartopy.feature.LAND,color=[.75,.75,.75],zorder=100)
# ax.set_title(t_name,loc='right',fontdict={'fontsize':24,'fontweight':'regular'})
# plt.plot([T_point1[0],T_point1[-1]],[T_point2[0], T_point2[-1]])
ax.plot(a,b,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(c,d,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(e,f,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(g,h,transform=PC,color='k',linestyle='--',linewidth=2.5)

ax.plot(a1,b1,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(c1,d1,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(e1,f1,transform=PC,color='k',linestyle='--',linewidth=2.5)
ax.plot(g1,h1,transform=PC,color='k',linestyle='--',linewidth=2.5)

# crs is PlateCarree -> we are explicitly telling axes, that we are creating bounds that are in degrees
ax.set_extent([0, 360, -80, -24], crs=PC)
ax.tick_params(axis='both', which='major', labelsize=28)
plt.tight_layout()
if 0:
    plt.savefig(w_path+'/ppt/'+save_name,
            facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(w_path+'/'+save_name)
plt.show()