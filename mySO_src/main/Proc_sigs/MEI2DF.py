import pandas as pd
import pickle
import matplotlib.pyplot as plt


npth='C:/Users/shjo/OneDrive/mySO/mySignals/Ori/SOI_1951_2022.csv'
wpth='C:/Users/shjo/OneDrive/mySO/mySignals/'

mySig=pd.read_csv(npth,header=None)

mySig=mySig.values.reshape(72*12)

mySig_index=pd.date_range('1951-01',periods=len(mySig),freq='1M')

mySig=pd.DataFrame({'SOI':mySig},index=mySig_index)

with open(wpth+'mySOI.pkl', 'wb') as f:
    pickle.dump(mySig, f)
