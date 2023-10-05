from scipy import io
import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')

pth='C:/Users/shjo/OneDrive/mySO/SCP_EEMD/OHC700_1993_2020_230E250E_60S50S/EEMD_results/'
wpth='C:/Users/shjo/OneDrive/mySO/SCP_EEMD/OHC700_1993_2020_230E250E_60S50S/'

eemd_lst=[i for i in os.listdir(pth) if i.endswith('.mat')]

try:
    os.mkdir(wpth+'Figs/')
except:
    pass
wpth=wpth+'Figs/'


for i,nm in zip(eemd_lst,eemd_lst):
    tmp=nm.split('.')[0].replace('_',' ').split(' ')
    nm_=tmp[1]+' '+tmp[2]+' '+tmp[3]+' '+tmp[4]+' '+'\n'+tmp[5]+' '+tmp[6]
    sm=tmp[0]+'_'+tmp[1]+'_'+tmp[2]+'_'+tmp[3]+'_'+tmp[4]+'_'+tmp[5]+'_'+tmp[6]
    
    myEEMD=io.loadmat(pth+i)['tmp']
    myTime=pd.date_range('1993-01',periods=myEEMD.shape[0],freq='1M')
    myEEMD=myEEMD.transpose()
        
    lbl_N=np.shape(myEEMD)[0]
    Label_size = 24
    Figsize=(12,2.1*lbl_N)
    Height_ratios=list(np.ones(lbl_N))
    xtick_location = myTime[24::12*5]
    xtick_labels = myTime.strftime('%Y')[24::12*5]

    fig, axs = plt.subplots(lbl_N,1,figsize=Figsize,constrained_layout = True,
                        sharex=True,gridspec_kw={'height_ratios': Height_ratios},dpi=200)
    axs[0].set_title(nm_,loc='right',fontdict={'fontsize':32,'fontweight':'regular'}, pad=20)
    for i,n in zip(myEEMD,range(lbl_N)):
        axs[n].plot(myTime,i,color='k',linewidth=2,zorder=0)
        axs[n].tick_params(axis='y', labelsize=Label_size)     
        axs[n].set_xticks(ticks=xtick_location)
        axs[n].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
        axs[n].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
        axs[n].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
        # axs[n].legend(fontsize=12,loc='upper right')
        
    if 1:
        # plt.savefig(wpth+'ppt/'+wname,
        #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
        plt.savefig(wpth+sm)
    plt.show()             
                        
                            
                         














