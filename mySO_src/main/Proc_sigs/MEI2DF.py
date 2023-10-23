import pandas as pd
import pickle
import matplotlib.pyplot as plt


npth='C:/Users/shjo/OneDrive/mySO/mySignals/Ori/meiv2_1979_2022.csv'
wpth='C:/Users/shjo/OneDrive/mySO/mySignals/'

mySig=pd.read_csv(npth,header=None)
mySig=mySig.values.reshape(44*12)

mySig_index=pd.date_range('1979-01',periods=len(mySig),freq='1M')

mySig=pd.DataFrame({'MEIv2':mySig},index=mySig_index)

with open(wpth+'myMEI.pkl', 'wb') as f:
    pickle.dump(mySig, f)
