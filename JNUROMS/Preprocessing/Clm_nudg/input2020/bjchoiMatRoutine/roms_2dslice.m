function [data,x,y,t,grd] = roms_2dslice(file,var,time,grd)
% Get a horiztonal slice of a 2-D variable out of a ROMS history, averages or restart file
% [data,x,y,t] = roms_2dslice(file,var,time,grd)
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file
%    grd (optional) is the structure of grid coordinates from get_roms_grd 
%
% Outputs
%    
%    data = the 2d slice at requested depth 
%    x,y = horizontal coordinates
%    t = time in days for the data
%
% John Wilkin

% open the history or averages file
nc = netcdf(file);

if isempty(nc{var})
  error([ 'Variable ' var ' is not present in file ' file])
end

% get the time
time_variable = nc{var}.time(:);
if isempty(time_variable)
  time_variable = 'scrum_time';
end
ot = nc{time_variable}(:);
time_units = nc{time_variable}.units(:);
if strcmp(lower(time_units(1:3)),'sec');
  ot = ot/86400;
end

% get the time
switch var
  case { 'h','f','pm','pn'}
    % dont need time
  otherwise    
    if time > length(ot)
      disp(['Requested time index ' int2str(time) ' not available'])
      disp(['There are ' int2str(length(ot)) ' time records in ' file])
      disp(['from day ' num2str(ot(1)) ' to day ' num2str(ot(length(ot)))])
      error(' ')
    end    
    t = ot(time);
end

% get the data to slice up
switch var
  case { 'h','f','pm','pn'}
    data = nc{var}(:);
  otherwise
    data = nc{var}(time,:,:);
    if ~size(data)
      error(['Variable ' var ' not found in ' file])      
    end
end

% check the grid information
if nargin<4 | (nargin==4 & isempty(grd))
  % no grd input given so try to get grd_file name from the history file
  grd_file = nc.grd_file(:);
  grd = get_roms_grid(grd_file,file);
else
  if isstr(grd)
    grd = get_roms_grid(grd,file);
  else
    % input was a grd structure but check that it includes the z values    
    % if ~isfield(grd,'z_r')
    %  error('grd does not contain z values');
    % end
  end
end

% close the nc file
close(nc)

switch var
  case {'ubar','sustr'}   
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
    
  case {'vbar','svstr'}
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v; 
    
  otherwise    
    % for zeta, Hsbl etc
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho; 
    
end

% Apply mask
dry = find(mask==0);
mask(dry) = NaN;
data = data.*mask;
