% Create a netcdf file of forcing condition data for ROMS
%
% It is assumed the data to be written are on a ROMS 2-D grid
%
% The ROMS 2-D data must be a structure named 'roms'
%   roms.time      = the times of the data
%   roms.base_date = COARDS format string describing base date

%   roms.data      = the data
%   roms.frcvname  = string with the variable name (This is required because 
%                    this script is being used such that  


% Additional data required in the workspace are:
%
%   time_variable = name of the time variable (string)
%   out_file = the name of the netcdf file generated here
%   grd_file (for the record)
%
% Optional data in the workspace
%   titlestr
%   sourcestr = something aobut where the data came from
%   details = something about computations relevant to these data
%
% John Wilkin
% Updated for EAC 15-Aug-2003

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
  nc.type = ncchar('ROMS forcing file');
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

nc('eta_rho') = eta_rho;
nc('xi_rho') = xi_rho;
nc('eta_u') = eta_u;
nc('xi_u') = xi_u;
nc('eta_v') = eta_v;
nc('xi_v') = xi_v;

% The unlimited dimension
time_dimension = 'time';
nc(time_dimension) = 0; % UNLIMITED

% time coordinates --------------------------------------------------------

theVarname = time_variable;
nc{theVarname} = ncfloat(time_dimension);
nc{theVarname}.long_name = ncchar('forcing observations time');
nc{theVarname}.units = ncchar(time_variable_units);
nc{theVarname}.field = ncchar([theVarname ', scalar, series']);

% forcing data ------------------------------------------------------------

theVarname = frcvarname;
switch theVarname
  case { 'sustr'}
    eta_dim = 'eta_u';
    xi_dim = 'xi_u';    
  case { 'svstr'}
    eta_dim = 'eta_v';
    xi_dim = 'xi_v';
  otherwise
    eta_dim = 'eta_rho';
    xi_dim = 'xi_rho';    
end
nc{theVarname} = ncfloat(time_dimension,eta_dim,xi_dim);
nc{theVarname}.long_name = ncchar(frcvarlongname);
nc{theVarname}.units = ncchar(frcvarunits);
nc{theVarname}.positive = ncchar(frcvarconvention);
nc{theVarname}.field = ncchar([theVarname ', scalar, series'])
nc{theVarname}.time = ncchar(time_variable);

close(nc)
