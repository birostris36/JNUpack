from netCDF4 import MFDataset,Dataset,num2date,date2num
import statsmodels.api as sm
from scipy.stats import  linregress
import numpy as np
import datetime as dt
import xarray as xr
import os



def myfitting(x,y,deg=1):
    '''
    x,y --> 1d numpy.array
    '''
    a, b = np.polyfit(x, y, deg=1)
    n = x.size
    y_est = a * x + b
    y_err = (y-y_est).std() * np.sqrt(1/n + (x - x.mean())**2 / np.sum((x - x.mean())**2))
    print('!!! y_est=ax+b !!!')
    return [a,b],y_err




import seaborn as sns
fig, ax = plt.subplots()

# plot manually calculated interval (std interval) --- the blue one
ax.plot(x, y_est, '-')
ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2,color='r')

# plot seaborn calculated interval (std interval, i.e. when ci=68.27) --- the orange one
# sns.regplot(x=x, y=y, ci=68.27)
plt.show()


NC=sst.assign_coords({'TT':('time',range(len(sst.time)))})
NC=NC.swap_dims({"time":"TT"})
NC_s=NC.polyfit(dim='TT',deg=1,skipna=False)
fit = xr.polyval(NC.TT, NC_s.polyfit_coefficients)
Coef=NC_s.polyfit_coefficients[0]
Coef_var=Coef.values*12 # (m/year)
data_dt=NC-fit
data_dt=data_dt.swap_dims({"TT":"time"})

# calculate interval manually using the formula
my=sst.values
T,A,O=my.shape
y=np.arange(T)
tmp=np.zeros([A,O])
for i in range(A):
    for j in range(O):
        if np.mean(my[:,i,j])!=np.mean(my[:,i,j]):
            tmp[i,j]=np.nan
        else:
            a,_=np.polyfit(y,my[:,i,j], deg=1)
            tmp[i,j]=a
        

