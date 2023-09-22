
import sys, os
# append path above path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from libs.plot import plot_subs 
from netCDF4 import Dataset
# import numpy as np
class Stb:
    
    def __init__(self,rpth,wpth):
        self.log_name=rpth+'ocean.output'
        self.stat_name=rpth+'run.stat'
        self.statNC_name=rpth+'run.stat.nc'
        self.timestep = rpth+'time.step'
        self.setting_info=rpth+'setting_info'
        self.wpth=wpth
        
    def pro_stat1(self):
        
        with open(self.stat_name) as f:
            Model_Log=f.readlines()
            Model_Log = Model_Log
            
        CLMN=['it','ssh_max','u_max','s_min','s_max']
        LOGS=pd.DataFrame(data=None,columns=CLMN)
        for i in Model_Log:
            tmp=[j for j in i.split(' ') if len(j)!=0]
            tmp_list=[]
            for k in [tmp[4],tmp[6],tmp[8],tmp[10]]:
                tmp_fac=10**float(k.split('D')[-1].split('+')[-1])
                tmp_val=float(k.split('D')[0])*tmp_fac
                tmp_list.append(tmp_val)

            tmpP=pd.DataFrame({'it':int(tmp[2]),'ssh_max':tmp_list[0],'u_max':tmp_list[1],'s_min':tmp_list[2],'s_max':tmp_list[3]},index=[0])
            LOGS=pd.concat([LOGS,tmpP],axis=0)
        LOGS=LOGS.reset_index(drop=True)
        
        Title_name='Model stability'
        time1 = LOGS['it'].astype('int')
        time2 = time1
        Label_size = 18
        LOGS.to_csv(self.wpth+'tset.csv')
        plot_subs(time1,time2,LOGS.drop('it',axis=1).round(decimals=8),Title_name,self.wpth,'Model_stability01',fig_bool=1)
            
    def pro_stat2(self):
        NC=Dataset(self.statNC_name)
        Header=NC.variables.keys()
        tmp_pd=pd.DataFrame({},columns=Header)
        for i in Header:
            tmp_pd[i]=NC[i][:]
        Title_name='Model stability (NC)'
        time1 = range(1,len(tmp_pd)+1) 
        time2 = time1
        Label_size = 18
        plot_subs(time1,time2,tmp_pd.round(decimals=8),Title_name,self.wpth,'Model_stability02',fig_bool=1)
        