% Create a netcdf file of initial condition data for ROMS
%
% It is assumed the data to be written are on a ROMS 3-D grid, but only the
% boundary segments will be used. (OK, so it uses a lot of space but it
% makes it easy to recycle code).
%
% The ROMS 3-D data must be a structure named 'roms'
%   roms.time = the times of the data
%   roms.base_date = COARDS format string describing base date
%   roms.temp      = temperature observations (such as from Ridgway's SBT
%                     method)
%   roms.salt
%   roms.{u,v,ubar,vbar} on their respective Arakawa-C grid locations
%
% Additional data required in the workspace are:
%
%   out_file = the name of the netcdf file generated here
%   grd_file (for the record)
%
%
% John Wilkin

disp(' ')
disp(['The output netcdf file (out_file) will be ' out_file])

% /matlab/toolbox/netcdf/@netcdf/netcdf.m
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
% ./matlab/toolbox/netcdf/ncchar.m

try
  nc.type = ncchar(roms.type);
catch
  nc.type = ncchar('ROMS INITIAL file');
end

if exist('titlestr')
  nc.title = ncchar(titlestr);
end
if exist('out_file')
  nc.out_file = ncchar(out_file);
end
if exist('grd_file')
  nc.grd_file = ncchar(grd_file);
end
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

eta_rho = size(roms.grd.lon_rho,1);
xi_rho = size(roms.grd.lon_rho,2);
eta_u = size(roms.grd.lon_u,1);
xi_u = size(roms.grd.lon_u,2);
eta_v = size(roms.grd.lon_v,1);
xi_v = size(roms.grd.lon_v,2);

s_rho = length(roms.grd.sc_r);
s_w = length(roms.grd.sc_w);

nc('eta_rho') = eta_rho;
nc('xi_rho') = xi_rho;
nc('eta_u') = eta_u;
nc('xi_u') = xi_u;
nc('eta_v') = eta_v;
nc('xi_v') = xi_v;
nc('s_rho') = s_rho;
nc('s_w') = s_w;
nc('one') = 1;

% The unlimited dimension
time_dimension = 'ocean_time';
nc(time_dimension) = 0; % UNLIMITED

if exist('generic_tracer')
  if generic_tracer
    other_time_dimension = 'ones_time'; % for the generic unscaled tracer
    nc(other_time_dimension) = 2; 
  end
end

% coordinates ------------------------------------------------------------

% THE TIME

theVarname = time_variable;
nc{theVarname} = ncfloat(time_dimension);
nc{theVarname}.long_name = ncchar('subsurface temp/salt observations time');
nc{theVarname}.units = ncchar(time_variable_units);
nc{theVarname}.field = ncchar([theVarname ', scalar, series']);

if exist('generic_tracer')
  if generic_tracer
    theVarname = 'ones_time';
    nc{theVarname} = ncfloat(other_time_dimension);
    nc{theVarname}.long_name = ncchar('generic unit tracer time');
    nc{theVarname}.units = ncchar(time_variable_units);
    nc{theVarname}.field = ncchar([theVarname ', scalar, series']);
  end  
end

% THE VERTICAL COORDINATE VARIABLES

theVarname = 'theta_s';
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate surface control parameter');
nc{theVarname}.units = ncchar('nondimensional');

theVarname = 'theta_b';
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate bottom control parameter');
nc{theVarname}.units = ncchar('nondimensional');

theVarname = 'Tcline';
nc{theVarname} = ncfloat;
nc{theVarname}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{theVarname}.units = ncchar('meter');

theVarname = 'hc';
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

% enter the values

nc{'theta_s'}(:) = roms.grd.theta_s;
nc{'theta_b'}(:) = roms.grd.theta_b;
nc{'Tcline'}(:)  = roms.grd.Tcline;
nc{'hc'}(:)      = roms.grd.hc;
nc{'sc_r'}(:)    = roms.grd.sc_r;
nc{'sc_w'}(:)    = roms.grd.sc_w;
nc{'Cs_r'}(:)    = roms.grd.Cs_r;
nc{'Cs_w'}(:)    = roms.grd.Cs_w;

% THE CLIMATOLOGY VARIABLES

theVarname = 'temp';
nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('potential temperature');
nc{theVarname}.units = ncchar('Celcius');
nc{theVarname}.field = ncchar('temperature, scalar, series');
nc{theVarname}.time = ncchar(time_variable);

theVarname = 'salt';
nc{theVarname} = ncfloat(time_dimension, 's_rho', 'eta_rho', 'xi_rho');
nc{theVarname}.long_name = ncchar('salinity');
nc{theVarname}.units = ncchar('PSU');
nc{theVarname}.field = ncchar('salinity, scalar, series');
nc{theVarname}.time = ncchar(time_variable);

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
  tracer(tnum).long_name = 'chlorophyll-phytoplankton ratio';
  tracer(tnum).units     = 'millimole nitrogen meter-3';
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
