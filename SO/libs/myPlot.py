import matplotlib.pyplot as plt
import os
import sys

plt.rcParams["font.family"] = 'Arial'


mySetting={
    'figsize': (10,3.7),
    'mylabel': [],
    'Label_size':18,
    'title':['test','right'],
    'xdomain_label_freq':[5::12*4],
    'fontParams':'Arial'
}



class figmaster:
    
    def __init__(self,mySetting):
        self.figsize=mySetting['figsize']
        self.mylabel=mySetting['mylabel']
        self.Label_size=mySetting['Label_size']
        self.xdomain_label=mySetting['xdomain_label_freq']
        self.title=mySetting['title']
        
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
        
        
        '''
            sbc = atm_to_nemo(pre_config=pre_config,
                      tools  = tools,
                      wrinc  = wrinc,
                      interp = interp,
                      agrif_mode=pre_config.agrif_mode)
        '''


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
        plt.savefig(wnpth+'/'+save_name,bbox_inches='tight')
    plt.show()


def plot_trend(y_est,):

    ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2,color='r')


import seaborn as sns
fig, ax = plt.subplots()

# plot manually calculated interval (std interval) --- the blue one
ax.plot(x, y_est, '-')

# plot seaborn calculated interval (std interval, i.e. when ci=68.27) --- the orange one
# sns.regplot(x=x, y=y, ci=68.27)
plt.show()
