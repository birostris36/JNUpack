function [Data,han] = roms_iview(file,var,time,iindex,grd)
% [data,han] = roms_iview(file,var,time,iindex,grd)
%
% file   = roms his/avg/rst etc nc file
% var    = variable to plot
% time   = time index into nc file
% iindex = iindex for slice
% grd can be 
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
%
% John Wilkin

if nargin < 5
  grd = [];
end

[data,z,lon,lat,t] = roms_islice(file,var,time,iindex,grd);

lonstr = [' - Lon ' num2str(mean(lon(:)),4)];

hant = pcolorjw(lat,z,data);

% pcolor plot of the variable

titlestr = ...
    {['file: ' strrep_(file) ],...
    [upper(var) ' - Day ' num2str(t) lonstr]};

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
