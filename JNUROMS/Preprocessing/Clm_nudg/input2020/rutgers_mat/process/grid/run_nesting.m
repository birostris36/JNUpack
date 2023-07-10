clear
% Lp=10;
% Mp=10;
% Gname='latter_grd.nc';
% c_grid(Lp,Mp,Gname,varargin)

f_path='/home/jhlee/model/ROMS_latest/project/NWP4/input/';
Ginp=[f_path 'NWP4_grd_3_20m_1.nc'];
% Ginp=[f_path 'NWP4_grd_3_org.nc'];
FGout=[f_path 'NWP4_finer_grd_1.nc'];
Gfactor=3;
Imin=45; Imax=81;
Jmin=67;  Jmax=90;
% Imin=29; Imax=101;
% Jmin=43;  Jmax=125;
S=coarse2fine(Ginp,FGout,Gfactor,Imin,Imax,Jmin,Jmax);

% f_path = '/home/jhlee/model/ROMS_latest/project/NWP4/input/';
% grd_path = [f_path 'NWP4_grd_3_20m.nc'];
% his_path = [f_path 'NWP4_ini_ECCO_20120125.nc'];
% CGout = get_roms_grid(grd_path, his_path, 1);

% F_S=grid_extract(grd_path, FGout,Imin,Imax,Jmin,Jmax);

% break
Gnames={Ginp,FGout};
Cname='/home/jhlee/model/ROMS_latest/project/NWP4/input/NWP4_ngc.nc';

[S,G] =contact(Gnames, Cname);

run('/home/jhlee/Matlab_tools/rutgers_mat/initial/d_roms2roms_NWP4.m');


