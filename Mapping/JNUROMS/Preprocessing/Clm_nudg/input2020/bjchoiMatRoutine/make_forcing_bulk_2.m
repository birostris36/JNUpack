%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS forcing file
%
%  Extrapole and interpole surface data to get surface boundary
%  conditions for ROMS (forcing netcdf file)
%
%  Data input format (netcdf):
%     taux(T, Y, X)
%     T : time [Months]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO Climate Data Library 
%                (Atlas of Surface Marine Data 1994)
%
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.DASILVA/
%
%  Pierrick Penven, IRD, 2002.
%  Byoung J. Choi, Rutgers, 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
%  Title - Grid file name - Forcing file name
%
title=['Forcing (COADS)'];
grdname='D:\Matlab\run\080211\roms_grd_wcsurge_080211.nc';
frcname='roms_wcsurge_forcing_swflux080320.nc';

coadst = (15:30:345);  % time for wind stress [days]
coadsc = 360;          % cycle for wind stress [days]

%
%  Fresh water fluxes (evaporation - precipitation)
%
swf_file='d:/matlab/roms_tool/data/COADS05/emp.cdf';
swf_name='emp';
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp(title)
%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
result=close(nc);
cosa = cos(angle);
sina = sin(angle);
%
% Create the forcing file
%
disp(' ')
disp(' Create the forcing file...')
create_forcing_bulk(frcname,grdname,title,coadst,...
               coadsc)
%
% Loop on time
%
nc=netcdf(frcname,'write');

for tindex=1:length(coadst)
  time=nc{'coads_time'}(tindex);
%
% coeff = mm/(3hour) -> centimeter day-1 (!!!!!)
%
  nc{'swflux'}(tindex,:,:)=0.8*ext_data(swf_file,swf_name,tindex,...
                                        lon,lat,time);

end

close(nc)
%