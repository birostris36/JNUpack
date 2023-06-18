%
% Pair Bulk Forcing
%
title=['Forcing (COADS)'];

grdname= 'D:\Matlab\run\080211\roms_grd_wcsurge_080211.nc';

%frcname='d:/matlab/roms_tool/make_forcing/bulk_heat/roms_eas_frc_tide_sin.nc';
frcname='roms_wcsurge_forcing_Pair080211.nc';
%
%  Wind stress
%
U_file='d:/matlab/roms_tool/data/COADS05/u3.cdf';
U_name='u3';
V_file='d:/matlab/roms_tool/data/COADS05/v3.cdf';
V_name='v3';
coadst = (15:30:345);  % time for wind stress [days]
coadsc = 360;          % cycle for wind stress [days]
%
%  Air temperature
%
sat_file='d:/matlab/roms_tool/data/COADS05/sat.cdf';
sat_name='sat';
%
%  Relative humidity
%
rh_file='d:/matlab/roms_tool/data/COADS05/rh.cdf';
rh_name='rh';
%
%  Short wave radiation
%
srf_file='d:/matlab/roms_tool/data/COADS05/shortrad.cdf';
srf_name='shortrad';
%
%
%  sea level pressure
%
slp_file='d:/matlab/roms_tool/data/COADS05/slp.cdf';
slp_name='slp';
%
%  Fresh water fluxes (evaporation - precipitation)
%
% swf_file='d:/matlab/roms_tool/data/COADS05/emp.cdf';
% swf_name='emp';
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
create_forcing_bulk_Pair(frcname,grdname,title,coadst,...
               coadsc);
%
% Loop on time
%
nc=netcdf(frcname,'write');
for tindex=1:length(coadst)
  time=nc{'coads_time'}(tindex);
  u=ext_data(U_file,U_name,tindex,...
             lon,lat,time);
  v=ext_data(V_file,V_name,tindex,...
             lon,lat,time);
%
%  Rotation (if not rectangular lon/lat grid)
%
  nc{'Uwind'}(tindex,:,:)=(u.*cosa + v.*sina);
  nc{'Vwind'}(tindex,:,:)=(v.*cosa - u.*sina);
end

for tindex=1:length(coadst)
  time=nc{'coads_time'}(tindex);
%
coeff = mm/(3hour) -> centimeter day-1 (!!!!!)

  nc{'swflux'}(tindex,:,:)=0.8*ext_data(swf_file,swf_name,tindex,...
                                        lon,lat,time);

end


for tindex=1:length(coadst)
  time=nc{'coads_time'}(tindex);
  nc{'Tair'}(tindex,:,:)=ext_data(sat_file,sat_name,tindex,...
                                  lon,lat,time);
  nc{'Qair'}(tindex,:,:)=ext_data(rh_file,rh_name,tindex,...
                                  lon,lat,time);
  nc{'Pair'}(tindex,:,:)=ext_data(slp_file,slp_name,tindex,...
                                  lon,lat,time);
  nc{'swrad'}(tindex,:,:)=ext_data(srf_file,srf_name,tindex,...
                                  lon,lat,time);
end

close(nc)