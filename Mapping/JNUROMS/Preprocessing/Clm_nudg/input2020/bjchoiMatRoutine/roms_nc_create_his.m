% function roms_nc_create_his(out_file,roms)
% Create a ROMS netcdf history/initial file
% 
% The file format is consistent with that assumed by a ROMS subsurface 
% temperature assimilation/nudging file, so the data must already be on a
% ROMS 3-D grid
%
% Data required in the workspace are:
%
%   Structure grd (see roms_get_grid)
%        
%   out_file = the name of the netcdf file generated here
%   grd_file = name of grid file for the record

%
% This script could be called from e.g. init_assim_from_sbt.m
%
% John Wilkin

disp(' ')
disp(['The output netcdf file (out_file) will be ' out_file])

switch noclobber
  case 1
    nc = netcdf(out_file,'noclobber');
  case 0
    nc = netcdf(out_file,'clobber');
  otherwise
end
if isempty(nc)
  error(['Failed to open ' out_file])
end

% global attributes

try
  nc.type = ncchar(roms.type);
catch
  nc.type = ncchar('ROMS history/initial file ');
end

nc.out_file = ncchar(out_file);
nc.grd_file = ncchar(grd_file);
if exist('sourcestr')
  nc.source = ncchar(sourcestr);
end
if exist('details')
  nc.details = ncchar(details);
end
if exist('reference_string')
  nc.reference = ncchar(reference_string); 
end
nc.history = ncchar(['Created by ' which(mfilename) ' - ' datestr(now)]);

% dimensions

eta_rho = size(grd.lon_rho,1);
xi_rho = size(grd.lon_rho,2);
eta_u = size(grd.lon_u,1);
xi_u = size(grd.lon_u,2);
eta_v = size(grd.lon_v,1);
xi_v = size(grd.lon_v,2);

s_rho = length(grd.sc_r);
s_w = length(grd.sc_w);

nc('xi_rho') = xi_rho;
nc('xi_u') = xi_u;
nc('xi_v') = xi_v;
nc('eta_rho') = eta_rho;
nc('eta_u') = eta_u;
nc('eta_v') = eta_v;
nc('s_rho') = s_rho;
nc('s_w') = s_w;
nc('one') = 1;
nc('N') = s_rho;
time_dimension = 'time';
nc(time_dimension) = 0; % UNLIMITED

% coordinates

% THE TIME

if ~exist('time_variable')
  time_variable = 'ocean_time';
  time_variable_units = 'seconds';
end
theVarname = time_variable;
nc{theVarname} = ncfloat(time_dimension);
nc{theVarname}.long_name = ncchar('subsurface temp/salt observations time');
nc{theVarname}.units = ncchar(time_variable_units);
nc{theVarname}.field = ncchar([theVarname ', scalar, series']);
  
% THE VERTICAL COORDINATE VARIABLES

theVarname = 'spherical';
%nc{theVarname} = ncchar('one');
nc{theVarname} = ncchar;
nc{theVarname}.long_name = ncchar('grid type logical switch');

theVarname = 'theta_s';
%nc{theVarname} = ncfloat('one');
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate surface control parameter');
nc{theVarname}.units = ncchar('nondimensional');

theVarname = 'theta_b';
%nc{theVarname} = ncfloat('one');
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate bottom control parameter');
nc{theVarname}.units = ncchar('nondimensional');

theVarname = 'Tcline';
%nc{theVarname} = ncfloat('one');
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{theVarname}.units = ncchar('meter');

theVarname = 'hc';
%nc{theVarname} = ncfloat('one');
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{theVarname}.units = ncchar('meter');

theVarname = 'sc_r';
nc{theVarname} = ncfloat('s_rho');
nc{theVarname}.long_name = ncchar('S-coordinate at RHO-points');
nc{theVarname}.units = ncchar('nondimensional');
nc{theVarname}.valid_min = -1;
nc{theVarname}.valid_max = 0;
nc{theVarname}.field = ncchar('sc_r, scalar');

theVarname = 'sc_w';
nc{theVarname} = ncfloat('s_w');
nc{theVarname}.long_name = ncchar('S-coordinate at W-points');
nc{theVarname}.units = ncchar('nondimensional');
nc{theVarname}.valid_min = -1;
nc{theVarname}.valid_max = 0;
nc{theVarname}.field = ncchar('sc_w, scalar');

theVarname = 'Cs_r';
nc{theVarname} = ncfloat('s_rho');
nc{theVarname}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{theVarname}.units = ncchar('nondimensional');
nc{theVarname}.valid_min = -1;
nc{theVarname}.valid_max = 0;
nc{theVarname}.field = ncchar('Cs_r, scalar');

theVarname = 'Cs_w';
nc{theVarname} = ncfloat('s_w');
nc{theVarname}.long_name = ncchar('S-coordinate stretching curves at W-points');
nc{theVarname}.units = ncchar('nondimensional');
nc{theVarname}.valid_min = -1;
nc{theVarname}.valid_max = 0;
nc{theVarname}.field = ncchar('Cs_r, scalar');

theVarname = 'Lev';
nc{theVarname} = ncint('s_rho');
nc{theVarname}.long_name = ncchar('output levels');
nc{theVarname}.units = ncchar('nondimensional');
nc{theVarname}.valid_min = 1;
nc{theVarname}.valid_max = s_rho;
nc{theVarname}.field = ncchar('Cs_w, scalar');

% enter the values

nc{'theta_s'}(:) = grd.theta_s;
nc{'theta_b'}(:) = grd.theta_b;
nc{'Tcline'}(:)  = grd.Tcline;
nc{'hc'}(:)      = grd.hc;
nc{'sc_r'}(:)    = grd.sc_r;
nc{'sc_w'}(:)    = grd.sc_w;
nc{'Cs_r'}(:)    = grd.Cs_r;
nc{'Cs_w'}(:)    = grd.Cs_w;
nc{'Lev'}(:)     = 1:s_rho;
nc{'spherical'}(:) = 'T';

% The horizontal coordinates

theVarname = 'h';
nc{theVarname} = ncfloat('eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('bathymetry at RHO-points');
nc{theVarname}.units = ncchar('meter');

theVarname = 'lon_rho';
nc{theVarname} = ncfloat('eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('longitude at RHO-points');
nc{theVarname}.units = ncchar('degree_east');

theVarname = 'lat_rho';
nc{theVarname} = ncfloat('eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('latitude at RHO-points');
nc{theVarname}.units = ncchar('degree_north');

theVarname = 'mask_rho';
nc{theVarname} = ncfloat('eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('mask on RHO-points');
nc{theVarname}.option_0 = ncchar('land');
nc{theVarname}.option_1 = ncchar('water');

% enter the values
nc{'h'}(:) = grd.h;
nc{'lon_rho'}(:) = grd.lon_rho;
nc{'lat_rho'}(:) = grd.lat_rho;
nc{'mask_rho'}(:) = grd.mask_rho;

% THE CLIMATOLOGY VARIABLES

% NOTE: These variable names (and the associated time variable names) need
% to correspond the the entires in roms:init_scalars.F

theVarname = 'temp';
nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('potential temperature');
nc{theVarname}.units = ncchar('Celcius');
nc{theVarname}.field = ncchar('temperature, scalar, series');
nc{theVarname}.time = ncchar(time_variable);

if dosalt
    theVarname = 'salt';
    nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
    nc{theVarname}.long_name = ncchar('salinity');
    nc{theVarname}.units = ncchar('PSU');
    nc{theVarname}.field = ncchar('salinity, scalar, series');
    nc{theVarname}.time = ncchar(time_variable);
end  

if doerrs
  theVarname = 'temp_err';
  nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
  nc{theVarname}.long_name = ncchar('subsurface temperature normalized error variance');
  nc{theVarname}.units = ncchar('nondimensional');
end

if dosalt
  if doerrs
    theVarname = 'salt_err';
    nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
    nc{theVarname}.long_name = ncchar('subsurface salinity normalized error variance');
    nc{theVarname}.units = ncchar('nondimensional');
  end
end

if dou
  theVarname = 'u';
  nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_u', 'xi_u');
  nc{theVarname}.long_name = ncchar('u-momentum component');
  nc{theVarname}.units = ncchar('meter second-1');
  nc{theVarname}.field = ncchar('u-velocity, scalar, series');
  nc{theVarname}.time = ncchar(time_variable);
  
  theVarname = 'v';
  nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_v', 'xi_v');
  nc{theVarname}.long_name = ncchar('v-momentum component');
  nc{theVarname}.units = ncchar('meter second-1');
  nc{theVarname}.field = ncchar('v-velocity, scalar, series');
  nc{theVarname}.time = ncchar(time_variable);
  
  theVarname = 'ubar';
  nc{theVarname} = ncfloat(time_dimension, 'eta_u', 'xi_u');
  nc{theVarname}.long_name = ncchar('vertically averaged u-momentum component');
  nc{theVarname}.units = ncchar('meter second-1');
  nc{theVarname}.field = ncchar('ubar-velocity, scalar, series');
  nc{theVarname}.time = ncchar(time_variable);
  
  theVarname = 'vbar';
  nc{theVarname} = ncfloat(time_dimension, 'eta_v', 'xi_v');
  nc{theVarname}.long_name = ncchar('vertically averaged v-momentum component');
  nc{theVarname}.units = ncchar('meter second-1');
  nc{theVarname}.field = ncchar('vbar-velocity, scalar, series');
  nc{theVarname}.time = ncchar(time_variable);

  theVarname = 'zeta';
  nc{theVarname} = ncfloat(time_dimension, 'eta_rho', 'xi_rho');
  nc{theVarname}.long_name = ncchar('free-surface');
  nc{theVarname}.units = ncchar('meter');
  nc{theVarname}.field = ncchar('free-surface, scalar, series');
  nc{theVarname}.time = ncchar(time_variable);

end  

if dou
  if douerrs
    theVarname = 'Evel';
    nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
    nc{theVarname}.long_name = ncchar('subsurface velocity normalized error variance');
    nc{theVarname}.units = ncchar('nondimensional');
  end
end

if donuts
  tnum = 1;
  tracer(tnum).name      = 'NO3';
  tracer(tnum).long_name = 'nitrate concentration';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
  tnum = tnum+1;
  tracer(tnum).name      = 'phytoplankton';
  tracer(tnum).long_name = 'phytoplankton concentration';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
  tnum = tnum+1;
  tracer(tnum).name      = 'zooplankton';
  tracer(tnum).long_name = 'zooplankton concentration';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
  tnum = tnum+1;
  tracer(tnum).name      = 'NH4';
  tracer(tnum).long_name = 'ammonium concentration';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
  tnum = tnum+1;
  tracer(tnum).name      = 'Ldetritus';
  tracer(tnum).long_name = 'large fraction detritus concentration';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
  tnum = tnum+1;
  tracer(tnum).name      = 'Sdetritus';
  tracer(tnum).long_name = 'small fraction detritus concentration';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
  tnum = tnum+1;
  tracer(tnum).name      = 'CPratio';
  tracer(tnum).long_name = 'chlorophyll concentration';
  tracer(tnum).units     = 'milligram meter-3';
  % add new tracers here
  % tnum = tnum+1;
  % tracer(tnum).name      = ' ';
  % tracer(tnum).long_name = ' ';
  % tracer(tnum).units     = ' ';
  for nuts=1:size(tracer,2)
    theVarname = [tracer(nuts).name];
    nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
    nc{theVarname}.long_name = ncchar(tracer(nuts).long_name);
    nc{theVarname}.units = ncchar(tracer(nuts).units);
    nc{theVarname}.field = ncchar([theVarname ', scalar, series']);
    nc{theVarname}.FillValue_ = 0.0;
    nc{theVarname}.time = ncchar(time_variable);    
  end  
end

% force write-to-disk of the configuration and coordinates
% will need to re-open the file to write the data

close(nc)
