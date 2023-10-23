import pandas as pd
import pickle

npth='C:/Users\shjo/OneDrive/mySO/mySignals/Ori/Ind_txt_SAM_1957_2020.txt'
wpth='C:/Users\shjo/OneDrive/mySO/mySignals/'

mySig=pd.read_csv(npth,header=None)

mySig_index=pd.date_range('1957-01',periods=len(mySig),freq='1M')

mySig.index=mySig_index

mySig=mySig.drop(columns=[0,1]).rename(columns={2:'SAM'})

with open(wpth+'mySAM.pkl', 'wb') as f:
    pickle.dump(mySig, f)


