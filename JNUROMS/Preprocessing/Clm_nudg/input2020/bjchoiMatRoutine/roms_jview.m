function [Data,han] = roms_jview(file,var,time,jindex,grd)
% [data,han] = roms_jview(file,var,time,jindex,grd)
% 
% file   = roms his/avg/rst etc nc file
% var    = variable to plot
% time   = time index into nc file
% jindex = jindex for slice
% grd can be 
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
%
% John Wilkin

if nargin < 5
  grd = [];
end

[data,z,lon,lat,t] = roms_jslice(file,var,time,jindex,grd);

latstr = [' - Lat ' num2str(mean(lat(:)),4)];

hant = pcolorjw(lon,z,data);

% pcolor plot of the variable

titlestr = ...
    {['file: ' strrep_(file) ],...
    [upper(var) ' - Day ' num2str(t) latstr]};

title(titlestr)

if nargout > 0
  Data.var = data;
  Data.lon = lon;
  Data.lat = lat;
  Data.z = z;
  Data.t = t;
end

if nargout > 1
  han = hant;
end
