import numpy as np
import matplotlib.pyplot as plt
from scipy import io
import pandas as pd

wpth='E:/_tmp/RMnT/'
wname='GECOHC_5_mod'

OHC=io.loadmat('E:/_tmp/RMnT/ohc_eemd.mat')

# print(OHC)
print(OHC['ohce'].shape)
time1=pd.date_range('1980-01','2019-01',freq='m').strftime('%Y-%m')
time2=pd.date_range('1980-01','2019-01',freq='m').strftime('%Y')

DATA=np.transpose(OHC['ohce'])

DATA=DATA[4:,:]

lbl_N=len(DATA)
print(lbl_N)
Label_size = 12
Figsize=(8,2.1*lbl_N)
Height_ratios=list(np.ones(lbl_N))
xtick_location = time1[::12*5]
xtick_labels = time2[::12*5]

fig, axs = plt.subplots(lbl_N,1,figsize=Figsize,constrained_layout = True,
                    sharex=True,gridspec_kw={'height_ratios': Height_ratios},dpi=200)

for i,n in zip(DATA,range(lbl_N)):
    axs[n].plot(time1,i,color='k',linewidth=2,zorder=0)
    axs[n].tick_params(axis='y', labelsize=Label_size)     
    axs[n].set_xticks(ticks=xtick_location)
    axs[n].set_xticklabels(xtick_labels, rotation=0, fontsize=Label_size, alpha=1)
    axs[n].tick_params(axis='x', direction='in', length=6, pad=8, labelsize=Label_size, labelcolor='k', top=True,width=1.)
    axs[n].tick_params(axis='y', direction='in', length=6, pad=8, labelsize=Label_size-3, width=1., color='k')
    axs[n].legend(fontsize=12,loc='upper right')
    
if 1:
    # plt.savefig(wpth+'ppt/'+wname,
    #         facecolor='none',edgecolor='none',bbox_inches='tight',transparent=True)
    plt.savefig(wpth+wname)
plt.show()

    
