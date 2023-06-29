# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 23:40:49 2023

@author: shjo9
"""
from netCDF4 import Dataset

nc=Dataset('G:/MODEL_DATA/SODA_INI/Ini_soda_05d_jhlee_198002_ref70_ice.nc','a')

nc.createVariable('hsn', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('ai', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('hi', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('ti', 'f4', ('ocean_time','eta_rho','xi_rho'))

nc.createVariable('tis', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('s0mk', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('t0mk', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('utau_iw', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('chu_iw', 'f4', ('ocean_time','eta_rho','xi_rho'))


nc.createVariable('tisrf', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('snow_thick', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('ageice', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('sig11', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('sig22', 'f4', ('ocean_time','eta_rho','xi_rho'))
nc.createVariable('sig12', 'f4', ('ocean_time','eta_rho','xi_rho'))

nc['tis'].long_name='???'
nc['tis'].units='???'
nc['s0mk'].long_name='???'
nc['s0mk'].units='???'
nc['t0mk'].long_name='???'
nc['t0mk'].units='???'
nc['utau_iw'].long_name='???'
nc['utau_iw'].units='???'
nc['chu_iw'].long_name='???'
nc['chu_iw'].units='???'

nc['hsn'].long_name='???'
nc['hsn'].units='???'
   
nc['ai'].long_name='???'
nc['ai'].units='???'
   
nc['hi'].long_name='???'
nc['hi'].units='???'
   
nc['ti'].long_name='???'
nc['ti'].units='???'
   

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



nc['tis'][:]=-10
nc['s0mk'][:]=0
nc['t0mk'][:]=0
nc['utau_iw'][:]=0.001
nc['chu_iw'][:]=0.001125


nc['hsn'][:]=0
nc['ai'][:]=0
nc['hi'][:]=0
nc['ti'][:]=0


nc['tisrf'][:]=0
nc['snow_thick'][:]=0

nc['ageice'][:]=0
nc['sig11'][:]=0
nc['sig22'][:]=0
nc['sig12'][:]=0



nc.close()
