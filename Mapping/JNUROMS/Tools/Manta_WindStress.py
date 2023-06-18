# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 14:27:07 2022

@author: birostris
@email : birostris36@gmail.com

Name : Manta_WindStress.py
Reference :
    This is python version of ra_windstr.m and ra_windstrcurl.m rewritten by birostris.
    Origianl author of ra_windstr.m and ra_windstrcurl.m is Ramkrushn S. Patel (ramkrushn.scrv89@gmail.com)
Description : 
    DESCRIPTION:  Function to compute wind stress from wind field data Based on Gill, 1982 
    Formula and a non-linear Cd based on Large and Pond (1981), modified for low wind 
    speeds (Trenberth et al., 1990)
"""
import numpy as np

def ra_windstr(u,v):
    '''
    INPUTS: 
    u = Zonal wind component [m/s], must be 2D
    v = Meridional wind component [m/s], must be 2D
    
    OUTPUT: 
    Tx = Zonal wind stress [N/m^2]
    Ty = Meridional wind stress [N/m^2]
    '''
    
    # Checking size of U,V
    if u.shape != u.shape:
        raise ValueError
    # Defining Constant
    roh = 1.2 # kg/m^3, air density
    
    # Computation of Wind Stresses
    [lt,ln] = u.shape 
    Tx,Ty = np.zeros([lt,ln]),np.zeros([lt,ln])
    Tx[Tx==Tx] = np.nan; Ty[Ty==Ty] = np.nan
    for ii in range(lt):
        for jj in range(ln):
            U = np.sqrt(u[ii,jj]**2+v[ii,jj]**2) # Wind speed
            if U <= 1:
                Cd = .00218
            elif (U > 1) or (U <=3):
                Cd=(0.62+1.56/U)*0.001
            elif (U > 3) or (U < 10):
                Cd = .00114
            else :
                Cd=(0.49+0.065*U)*0.001
            Tx[ii,jj] = Cd*roh*U*u[ii,jj] # kg/m^3*m/s*m/s= N/m^2
            Ty[ii,jj] = Cd*roh*U*v[ii,jj] # kg/m^3*m/s*m/s= N/m^2
    return Tx, Ty


def ra_windstrcurl(lat,lon,u,v):
    '''
    '''
    # lat,lon = lat_re1.values,lon_re1.values
    # u_,v_ = u_re1.values, v_re1.values
    # u,v= u_[0],v_[0]
    # Degree tp radian
    rad = np.pi/180
    # Wind stresses computation
    Tx,Ty = ra_windstr(u, v)

    # Computation of curl
    lt, ln = u.shape
    a=np.diff(lat)
    aa = np.nan*np.zeros([len(a)-1,1])
    
    for ii in range(len(a)-1):
        if a[ii]==a[ii+1]:
            aa[ii]=a[ii]
        else :
            print('Latitude difference is not consistance')
        dlat = np.mean(aa)
    del ii 
    deltay = dlat*111176
    curlZ,long = np.zeros([lt,ln]), np.zeros([lt,ln])
    curlZ[curlZ==curlZ] = np.nan 
    long[long==long] = np.nan
    for ii in range(lt):
        for jj in range(ln):
            long[ii,jj] = lon[jj]*111176*np.cos(lat[ii]*rad)
            # long(i,j)=lon(j)*6378137*rad*cos(lat(i)*rad)
            # [m] earth radious in meters= 6,378,137.0 m.. from wikipedia.

    del ii,jj
    
    # Centeral difference method in x and y
    for ii in range(2-1,lt-1):
        for jj in range(2-1,ln-1):
            curlZ[ii,jj] = (Ty[ii, jj+1]-Ty[ii, jj-1])/(2*(long[ii, jj+1]-long[ii, jj-1])) \
                - (Tx[ii+1, jj]-Tx[ii-1, jj])/(2*deltay)
    
    del ii,jj
    
    # Forward difference method in x and y
    for jj in range(ln-1):
        curlZ[0,jj] = (Ty[0, jj+1]-Ty[0, jj])/(long[0, jj+1]-long[0, jj]) \
            - (Tx[1, jj]-Tx[0, jj])/deltay 

    for ii in range(lt-1):
        curlZ[ii,0] = (Ty[ii, 1]-Ty[ii, 0])/(long[ii, 1]-long[ii, 0]) \
            - (Tx[ii, 1]-Tx[ii, 0])/deltay

    del ii,jj
    curlZ[0,ln-1] = curlZ[0, ln-2]
    
    # Backward difference method in x and y

    for ii in range(2-1,lt):
        curlZ[ii,ln-1]=(Ty[ii, ln-1]-Ty[ii, ln-2])/(long[ii, ln-1]-long[ii, ln-2])\
            - (Tx[ii, ln-1]-Tx[ii-1, ln-1])/deltay

    for jj in range(2-1,ln-1):
        curlZ[lt-1,jj] = (Ty[lt-1, jj]-Ty[lt-1, jj-1])/(long[lt-1, jj]-long[lt-1, jj-1]) \
            - (Tx[lt-1, jj]-Tx[lt-2, jj])/deltay 

    del ii,jj

    curlZ[lt-1,0] = curlZ[lt-1,lt-2]
    
    return curlZ




'''
% Tau=Rho*Cd*(speed)^2; Tx=Rho*Cd*Speed*u; Ty=Rho*Cd*Speed*v
%===========================================================%
% RA_WINDSTR  $Id: ra_windstr.m, 2014/10/29 $
%          Copyright (C) CORAL-IITKGP, Ramkrushn S. Patel 2014.
%
% AUTHOR: 
% Ramkrushn S. Patel (ramkrushn.scrv89@gmail.com)
% Roll No: 13CL60R05
% Place: IIT Kharagpur.
% This is a part of M. Tech project, under the supervision of DR. ARUN CHAKRABORTY
%===========================================================%
%
% USAGE: [Tx, Ty]=ra_windstr(u,v)
%  
% DESCRIPTION:  Function to compute wind stress from wind field data Based on Gill, 1982 
% Formula and a non-linear Cd based on Large and Pond (1981), modified for low wind 
% speeds (Trenberth et al., 1990)
% 
% INPUTS: 
% u = Zonal wind component [m/s], must be 2D
% v = Meridional wind component [m/s], must be 2D
%
% OUTPUT: 
% Tx = Zonal wind stress [N/m^2]
% Ty = Meridional wind stress [N/m^2]
% 
% DISCLAIMER: 
% Albeit this function is designed only for academic purpose, it can be implemented in 
% research. Nonetheless, author does not guarantee the accuracy.
% 
% REFERENCE:
% A.E. Gill, 1982, 밃tmosphere-Ocean Dynamics�, Academy Press, Vol. 30.
% W. G. Large & S. Pond., 1981,밢pen Ocean Measurements in Moderate to Strong Winds�, 
% J. Physical Oceanography, Vol. 11, pp. 324 - 336.
% K.E. Trenberth, W.G. Large & J.G. Olson, 1990, 밫he Mean Annual Cycle in Global Ocean 
% Wind Stress�, J.Physical Oceanography, Vol. 20, pp. 1742 � 1760.
%
% ACKNOWLEDGMENT:
% Author is eternally grateful to MathWorks for providing in built functions. 
'''







