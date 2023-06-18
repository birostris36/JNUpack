function [data,t] = roms_tslice(file,var,time,method)
% Get a 3-d field (time slice) out of a ROMS history, averages or restart file
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file
%
% Outputs
%    
%    data = the 3d slice at requested time
%    t = time in days for the data
%
% John Wilkin
%
% BUGS: To use ACOM snapshots to fill where there are no SBT data you eed to
% be using the 'nearest' option. There is little checking to see that the
% time handling is right. 

warning(['This has not been checked carefully for time handling'])

% open the history or averages file
nc = netcdf(file);

% find the time
time_variable = nc{var}.time(:);
ot = nc{time_variable}(:);
time_units = nc{time_variable}.units;
if strcmp(lower(time_units(1:3)),'sec');
  ot = ot/86400;
end

disp(['There are ' int2str(length(ot)) ' time records in ' file])
disp([' from day ' num2str(ot(1)) ' to day ' num2str(ot(length(ot)))])

% check for case of only one time slice and return it (not need for any
% interpolation 
if length(size(nc{var}))<4
  warning(['There is only one time slice in ' file])
  t = ot;
  return
end

% if the data are time cyclic we need to handle the possibility of
% finding a requested time that falls in the wrap-around interval
if length(nc{time_variable}.cycle_length)
  cycle_length = nc{time_variable}.cycle_length(:);
  disp([' and the data have attribute cycle_length of ' num2str(cycle_length)])
  lot = length(ot);
  ot = [ot(lot)-cycle_length; ot; cycle_length+ot(1)];
  data = nc{var}(:);
  data = [data(lot,:,:,:); data; data(1,:,:,:)];
else
  cycle_length = 0;
end

% default is to find the nearest time slice
% The other option allowed is 'linear' interpolation between time slices
if nargin < 4
  method = 'nearest';
end

tdiff = abs(ot-time);
[tdiff_sorted,itime] = sort(tdiff);
if(min(tdiff_sorted)==0)
  % no need to interpolate
  method = 'nearest';
end

switch method
  
  case 'nearest'
    
    % find nearest time to analysis date
    itime = itime(1); % it was already sorted
    
    % get the data for this time
    if cycle_length
      data = data(itime,:,:,:);
    else
      data = nc{var}(itime,:,:,:);
    end
    t = ot(itime);

  case 'linear'
    
    % linearly interpolate in time between nearest times
    itime = sort(itime([1 2])); % sort again to get indices in ascending order
    
    % interpolate
    data2 = data(itime,:,:,:);
    a = (time-ot(itime(1)))/diff(ot(itime));
    data = (1-a)*data2(1,:,:,:)+a*data2(2,:,:,:);    
    t = ot(itime);

    if abs(a)>1
      % then time is not between the bracketing values, which could be 
      % a requested date outside a cycle_length, or from the wrong history file
      error(['Time was not between bracketing values - extrapolation not allowed'])
    end
    
end

data = squeeze(data);

% close the nc file
close(nc)
