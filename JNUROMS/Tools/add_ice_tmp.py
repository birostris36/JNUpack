# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 23:40:49 2023

@author: shjo9
"""
from netCDF4 import Dataset

nc=Dataset('G:/MODEL_DATA/SODA_INI/Ini_soda_05d_jhlee_198002_ref70_ice.nc','a')

nc.createVariable('tisrf', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('snow_thick', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('ageice', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('sig11', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('sig22', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('sig12', 'f4', ('ocean_time','eta_rho','xi_rho'))

nc['tisrf'].long_name='???'
nc['tisrf'].units='???'
   
nc['snow_thick'].long_name='???'
nc['snow_thick'].units='???'

nc['ageice'].long_name='???'
nc['ageice'].units='???'

nc['sig11'].long_name='???'
nc['sig11'].units='???'

nc['sig22'].long_name='???'
nc['sig22'].units='???'

nc['sig12'].long_name='???'
nc['sig12'].units='???'


nc['tisrf'][:]=0
nc['snow_thick'][:]=0

nc['ageice'][:]=0
nc['sig11'][:]=0
nc['sig22'][:]=0
nc['sig12'][:]=0



nc.close()
