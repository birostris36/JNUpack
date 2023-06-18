import  numpy  as np
import matplotlib.pyplot as plt
import xarray as xr
import os

tmp_path='G:/Models/TK0525ED_CLM/'
save_path='D:/tmp/'



ncs=[tmp_path+i for i in os.listdir(tmp_path) if i.endswith(".nc")]



Sample_nc=xr.open_dataset(ncs[0])
ZETA=Sample_nc.zeta

plt.figure(1)
ZETA[0].plot()
plt.savefig(save_path+'tm1')



Sample_nc=xr.open_dataset(ncs[-1])
ZETA=Sample_nc.zeta



plt.figure(2)
ZETA[0].plot()
plt.savefig(save_path+'tm2')
plt.show()