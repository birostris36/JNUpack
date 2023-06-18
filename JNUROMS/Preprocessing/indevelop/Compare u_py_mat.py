# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 17:29:33 2023

@author: shjo9
"""

import matplotlib.pyplot as plt
from scipy import io
from copy import deepcopy



# ubar_1=deepcopy(DATA)[0]
# ubar_2=deepcopy(data_)[0]
# ubar_3=deepcopy(data)[0]

# vbar_1=deepcopy(DATA)[0]
# vbar_2=deepcopy(data_)[0]
# vbar_3=deepcopy(data)[0]

# dv[0][0]
dz
zv[0][0]
tmp_v[0,:,0,:].shape



mat_pth='D:/OneDrive/'

mat01=io.loadmat(mat_pth+'ext_167.mat')
mat_ubar_1=mat01['ubar']
mat_vbar_1=mat01['vbar']
mat_ubar_1[mat_ubar_1<-100]=np.nan
mat_vbar_1[mat_vbar_1<-100]=np.nan

# mat01['dv'][0]
mat01['dz'][0]
mat01['zv'][0]
mat01['v1'][0][0]



mat01['dv'][0]
dv[0][0]

plt.figure(dpi=150)
plt.plot(mat01['dv'][0],linewidth=1,label='mat_ubar_1')
plt.plot(dv[0][0],linewidth=1,label='py_ubar_1')
plt.legend()
plt.show()













plt.figure(dpi=150)
plt.plot(zv[0][0],linewidth=1,label='py_ubar_1')
plt.plot(mat01['zv'][0],linewidth=1,label='mat_ubar_1')
plt.legend()
plt.show()

NN=1

plt.figure(dpi=150)
plt.plot(mat01['v1'][NN,0,:],linewidth=1,label='mat_ubar_1')
plt.plot(tmp_v[0,NN,0,:],linewidth=1,label='py_ubar_1')
plt.legend()
plt.show()

mat01=io.loadmat(mat_pth+'ext_181.mat')
mat_ubar_3=mat01['ubar']
mat_vbar_3=mat01['vbar']

plt.figure(dpi=150)
plt.plot(mat_ubar_1[0],linewidth=1,label='mat_ubar_1')
plt.plot(ubar_1[0],linewidth=1,label='py_ubar_1')
plt.legend()
plt.show()


plt.figure(dpi=150)
plt.plot(mat_ubar_3[0],linewidth=1,label='mat_ubar_3')
plt.plot(ubar_3[0],linewidth=1,label='py_ubar_3')
plt.legend()
plt.show()

plt.figure(dpi=150)
plt.plot(mat_vbar_1[0],linewidth=1,label='mat_vbar_1')
plt.plot(vbar_1[0],linewidth=1,label='py_vbar_1')
plt.legend()
plt.show()



plt.figure(dpi=150)
plt.plot(mat_vbar_3[0],linewidth=1,label='mat_vbar_3')
plt.plot(vbar_3[0],linewidth=1,label='py_vbar_3')
plt.legend()
plt.show()





A.shape










