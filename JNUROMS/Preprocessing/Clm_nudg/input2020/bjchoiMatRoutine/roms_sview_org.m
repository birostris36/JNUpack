function [thedata,thegrid,han] = roms_sview(file,var,time,k,grd,vec_d,uscale,varargin)
% [thedata,thegrid,han] = roms_zview(file,var,time,k,grd,vec_d,uscale,varargin)
% 
% Inputs:
%
% file  = roms his/avg/rst etc nc file
% var   = variable to plot
% time  = time index into nc file
% k     = index of vertical (s-coordinate) level of horizontal slice 
%       if 0 any vector plot will be for ubar,vbar
% grd can be 
%       grd structure (from roms_get_grid)
%       grd_file name
%       [] (will attempt to get grid from roms file)
% vec_d = density (decimation factor) of velocity vectors to plot over 
%       if 0 no vectors are plotted
% varargin are quiver arguments passed on to roms_quiver
%
% Outputs:
% 
% thedata = structure of pcolored data and velocities
% thegrid = roms grid structure
% han = structure of handles for pcolor, quiver and title objects
%
% John Wilkin

if nargin == 0
  error('Usage: roms_sview(file,var,time,k,grd,vec_d,uscale,varargin)');
end
  
if nargin < 5
  grd = [];
end

% a nasty little trick to allow me to send a caxis range through the input
if isstruct(var)
  cax = var.cax;  
  var = var.name;
else 
  cax = 'auto';
end

varlabel = var;

% check the grid information

% pcolor plot of the variable

nc = netcdf(file);
switch var
  case { 'ubar','vbar','zeta','Hsbl','h','f','pm','pn',...
	'swrad','SST','dQdSST','shflux','swflux','SSS',...
	'sustr','svstr','Uwind','Vwind','Tair','Pair','Qair',...
	'lwrad','lwrad_down','lwrad_net','sensible','latent'}
    data = nc{var}(time,:);
    depstr = [];
  case 'stress'
    error('option not debugged yet')
    datau = nc{'sustr'}(time,:);
    datau = datau(:,[1 1:end end]);
    datau = av2(datau')';
    datav = nc{'svstr'}(time,:);
    datav = datav([1 1:end end],:);
    datav = av2(datav);
    data = abs(datau+sqrt(-1)*datav);
    depstr = [ ' - at surface '];
  case 'wind'
    datau = nc{'Uwind'}(time,:);
    datav = nc{'Vwind'}(time,:);
    data = abs(datau+sqrt(-1)*datav);
    depstr = [ ' - 10 m above surface '];
    var = 'Uwind'; % for time handling
  otherwise    
    data = nc{var}(time,k,:);
    depstr = [ ' - Level ' int2str(k) ' '];
end
close(nc)

switch var
  case { 'u','ubar','sustr'}
    mask = change(grd.mask_u,'==',0,NaN);
    x = grd.lon_u;
    y = grd.lat_u;
  case { 'v','vbar','svstr'}
    mask = change(grd.mask_v,'==',0,NaN);
    x = grd.lon_v;
    y = grd.lat_v;
  otherwise
    mask = grd.mask_rho_nan;
    x = grd.lon_rho;
    y = grd.lat_rho;
end

nc = netcdf(file);
t = nc{nc{var}.time(:)}(time);
time_units = nc{nc{var}.time(:)}.units(:);
if strcmp(time_units(1:3),'sec')
  t = t/86400;
end
got_date = 0;
try
  got_date = 0;
  basedate = parsetnc(time_units);
  if basedate(1)>0
    tdate = julian(basedate)+t;
    got_date = 1;
  end
catch
  disp([ ' Could not parse base date from ' file ' for ' var])
  got_date = 0;
end
close(nc)

hanpc = pcolorjw(x,y,data.*mask);

caxis(cax);

if nargin > 5
  if vec_d

    nc = netcdf(file);
    % add vectors
    % ! sorry, this doesn't allow for {u,v}bar vectors on a 3d variable
    if k>0 
      u = nc{'u'}(time,k,:);
      v = nc{'v'}(time,k,:);
      depstr = [depstr ' - Vectors at level ' int2str(k) ' '];
    else  
      u = nc{'ubar'}(time,:);
      v = nc{'vbar'}(time,:);
      % a forcing file won't have u,v ...
      if isempty(u)
        u = nc{'sustr'}(time,:);
        v = nc{'svstr'}(time,:);
	depstr = [depstr ' - Wind stress vectors '];
      else
	depstr = [depstr ' - Depth average velocity vectors '];
      end
    end
    close(nc)
    if nargin < 7
      uscale = 1;
    end    
    hanquiver = roms_quiver(grd.lon_rho,grd.lat_rho,u,v,grd.angle,...
	vec_d,uscale,varargin{:});
  end
end

amerc

try
  if findstr('leeuwin',grd.grd_file)
    gebco_eez(0,'k')
  elseif findstr('eauc',grd.grd_file)
    plotnzb
  elseif findstr('nena',grd.grd_file)
    plotnenacoast(2,'k')
  end
catch
end

if nargout > 0
  thedata.data = data;
  thedata.t = t;
  try
     nc = netcdf(file);
     thedata.base_date = nc{nc{var}.time(:)}.units(:);
     close(nc)
  catch
  end   
  if nargin > 5
    if vec_d
      thedata.u = u;
      thedata.v = v;
    end
  end
end
if nargout > 1
  thegrid = grd;
end

titlestr{1} = ['file: ' strrep_(file)];
titlestr{2} = [upper(varlabel) ' - Day ' num2str(t) depstr];
if got_date
  titlestr{2} = [titlestr{2} ' - Date ' greg2str(tdate,0)];
end

hantitle = title(titlestr);

if nargout > 2
  han.title = hantitle;
  han.pcolor = hanpc;
  if exist('hanquiver')
    han.quiver = hanquiver;
  end
end

