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

def myXfitting2d(Xdata,time):
    '''
    --> 3D xarray.dataarray (myDATA series)
    --> time --> Name of time (str)
    '''
    NC=Xdata.assign_coords({'TT':(time,range(len(Xdata[time])))})
    NC=NC.swap_dims({time:"TT"})
    NC_s=NC.polyfit(dim='TT',deg=1,skipna=False)
    fit = xr.polyval(NC.TT, NC_s.polyfit_coefficients)
    Coef=NC_s.polyfit_coefficients[0]
    Coef_var=Coef.values*12 # (m/year)
    data_dt=NC-fit
    data_dt=data_dt.swap_dims({"TT":time})
    data_dt=data_dt.drop_vars('TT')
    return data_dt, Coef


def myXRegress(Xdata,myIdx):
    '''
    --> 3D xarray.dataarray (myDATA series)
    --> time --> Name of time (str)
    '''
    NC=Xdata.assign_coords({'myIdx':('time',myIdx)})
    NC=NC.swap_dims({'time':"myIdx"})
    NC_s=NC.polyfit(dim='myIdx',deg=1,skipna=False)
    fit = xr.polyval(NC.myIdx, NC_s.polyfit_coefficients)
    Coef=NC_s.polyfit_coefficients[0]
    return Coef



def myfitting2d(data3d):
    '''
    Inputs --> data3d (time, lat,lon)
    '''
    # calculate interval manually using the formula
    T,A,O=data3d.shape
    y=np.arange(T)
    fit=np.zeros([A,O])
    for i in range(A):
        for j in range(O):
            if np.mean(data3d[:,i,j])!=np.mean(data3d[:,i,j]):
                fit[i,j]=np.nan
            else:
                a,_=np.polyfit(y,data3d[:,i,j], deg=1)
                fit[i,j]=a
    return fit

def myfitting2d_sttcs(data3d,threshold=0.05):
    '''
    Inputs --> data3d (time, lat,lon)
    '''
    # calculate interval manually using the formula
    T,A,O=data3d.shape
    y=np.arange(T)
    slope=np.zeros([A,O])
    intercept,r_value,p_value,std_err,Smask=\
        np.zeros_like(slope),np.zeros_like(slope),np.zeros_like(slope),\
            np.zeros_like(slope),np.zeros_like(slope)
    for i in range(A):
        for j in range(O):
            if np.mean(data3d[:,i,j])!=np.mean(data3d[:,i,j]):
                slope[i,j],intercept[i,j],r_value[i,j],p_value[i,j],std_err[i,j],Smask[i,j]=\
                    np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
            else:
                slope[i,j], intercept[i,j], r_value[i,j], p_value[i,j], std_err[i,j]=\
                    linregress(y,data3d[:,i,j])
                if p_value[i,j]<=threshold or p_value[i,j]>=1-threshold:
                    Smask[i,j]=0
                else:
                    Smask[i,j]=np.nan

    return slope, intercept, r_value, p_value, std_err, Smask

def myRegress3d_sttcs(sig,data3d,threshold=0.05):
    '''
    Inputs --> data3d (time, lat,lon)
    '''
    # calculate interval manually using the formula
    T,A,O=data3d.shape
    y=sig
    slope=np.zeros([A,O])
    intercept,r_value,p_value,std_err,Smask=\
        np.zeros_like(slope),np.zeros_like(slope),np.zeros_like(slope),\
            np.zeros_like(slope),np.zeros_like(slope)
    for i in range(A):
        for j in range(O):
            if np.mean(data3d[:,i,j])!=np.mean(data3d[:,i,j]):
                slope[i,j],intercept[i,j],r_value[i,j],p_value[i,j],std_err[i,j],Smask[i,j]=\
                    np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
            else:
                slope[i,j], intercept[i,j], r_value[i,j], p_value[i,j], std_err[i,j]=\
                    linregress(y,data3d[:,i,j])
                if p_value[i,j]<=threshold or p_value[i,j]>=1-threshold:
                    Smask[i,j]=0
                else:
                    Smask[i,j]=np.nan

    return slope, intercept, r_value, p_value, std_err, Smask


def myfitting1d_sttcs(data2d,threshold=0.05):
    '''
    Inputs --> data3d (time, lat,lon)
    '''
    # calculate interval manually using the formula
    T,A=data2d.shape
    y=np.arange(T)
    slope=np.zeros([A])
    intercept,r_value,p_value,std_err,Smask=\
        np.zeros_like(slope),np.zeros_like(slope),np.zeros_like(slope),\
            np.zeros_like(slope),np.zeros_like(slope)
    for i in range(A):
        if np.mean(data2d[:,i])!=np.mean(data2d[:,i]):
            slope[i],intercept[i],r_value[i],p_value[i],std_err[i],Smask[i]=\
                np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
        else:
            slope[i], intercept[i], r_value[i], p_value[i], std_err[i]=\
                linregress(y,data2d[:,i])
            # print(p_value)
            if p_value[i]<=threshold or p_value[i]>=1-threshold:
                Smask[i]=0
            else:
                Smask[i]=np.nan

    return slope, intercept, r_value, p_value, std_err, Smask

'''

import seaborn as sns
fig, ax = plt.subplots()


# plot manually calculated interval (std interval) --- the blue one
ax.plot(x, y_est, '-')
ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2,color='r')

plt.show()

'''
        
