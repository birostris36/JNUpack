function [data,x,y,t,grd] = roms_zslice(file,var,time,depth,grd)
% Get a constant-z slice out of a ROMS history, averages or restart file
% [data,x,y] = roms_zslice(file,var,time,depth,grd)
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file
%    depth = depth in metres of the required slice
%    grd (optional) is the structure of grid coordinates from get_roms_grd 
%
% Outputs
%    
%    data = the 2d slice at requested depth 
%    x,y = horizontal coordinates
%    t = time in days for the data
%
% John Wilkin

depth = -abs(depth);

% open the history or averages file
nc = netcdf(file);

if isempty(nc{var})
  error([ 'Variable ' var ' is not present in file ' file])
end

% get the time
time_variable = nc{var}.time(:);
if isempty(time_variable)
  time_variable = 'ocean_time'; % default used to be scrum_time and this 
                 % might be required to plot some really old output file   
else
  ot = nc{time_variable}(:);
  time_units = nc{time_variable}.units(:);
  if strcmp(lower(time_units(1:3)),'sec');
    ot = ot/86400;
  end
  if time > length(ot)
    disp(['Requested time index ' int2str(time) ' not available'])
    disp(['There are ' int2str(length(ot)) ' time records in ' file])
    disp(['from day ' num2str(ot(1)) ' to day ' num2str(ot(length(ot)))])
    error(' ')
  end
  t = ot(time);
end

% get the data to slice up

data = nc{var}(time,:,:,:);
if ~size(data)
  error(['Variable ' var ' not found in ' file])
end

% check the grid information
if nargin<5 | (nargin==5 & isempty(grd))
  % no grd input given so try to get grd_file name from the history file
  grd_file = nc.grd_file(:);
  grd = get_roms_grid(grd_file,file);
else
    if isstr(grd)
        %     grd = get_roms_grid(grd,file);
        grd = roms_get_grid(grd,file);
    else
        % input was a grd structure but check that it includes the z values
    if ~isfield(grd,'z_r')
      error('grd does not contain z values');
    end
  end
end

% close the nc file
close(nc)

% interpolate to requested depth

% (the following code could be replaced by a call to roms_zslice_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make 2d 'ribbon' out of the data
[N L M] = size(data);
data = data(:,1:L*M);

% make 2d 'ribbon' out of the depth coordinates
z_r = grd.z_r;
switch var
  case 'u'    
    % average z_r to Arakawa-C u points
    zM = size(z_r,2);
    zMm = zM-1;
    zL = size(z_r,3);
    zLm = zL-1;
    z = 0.5*(z_r(:,:,1:zLm)+z_r(:,:,2:zL));
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
    
  case 'v'
    % average z_r to Arakawa-C v points
    zM = size(z_r,2);
    zMm = zM-1;
    zL = size(z_r,3);
    zLm = zL-1;
    z = 0.5*(z_r(:,1:zMm,:)+z_r(:,2:zM,:));
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v; 
    
  otherwise    
    % for temp, salt, rho, w
    z = z_r;
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho; 
    if size(data,1) ~= size(z,1)
      % trap the var='omega' case
      % but omega can be N or N+1 depending on whether a rst or his file
      z = grd.z_w;
    end
    
end

z = reshape(z,[size(z,1) size(z,2)*size(z,3)]);

% code lifted from omviz/scrum_zslice:

% pad out bottom and surface z values with -Inf and 0 respectively
z = [-Inf*ones([1 L*M]); z; zeros([1 L*M])];

% pad out bottom and surface data values
data = [NaN*ones([1 L*M]); data; data(N,:)];

z = flipud(z);
data = flipud(data);

% Find the indices of data values that have just greater depth than
% depth

zg_ind = find(diff(z<depth)~=0);
zg_ind = zg_ind + [0:1:length(zg_ind)-1]';
data_greater_z = data(zg_ind);
depth_greater_z = z(zg_ind);
        
% Find the indices of the data values that have just lesser depth
% than depth
zl_ind = find(diff(z>depth)~=0);
zl_ind = zl_ind + [1:1:length(zg_ind)]';
data_lesser_z = data(zl_ind);
depth_lesser_z = z(zl_ind);
        
% Interpolate between the data values.
alpha = (depth-depth_greater_z)./(depth_lesser_z-depth_greater_z);
data_at_depth = (data_lesser_z.*alpha)+(data_greater_z.*(1-alpha));
data = reshape(data_at_depth,[L M]);

% Apply mask to catch shallow water values where the z interpolation does
% not create NaNs in the data
dry = find(mask==0);
mask(dry) = NaN;
data = data.*mask;
