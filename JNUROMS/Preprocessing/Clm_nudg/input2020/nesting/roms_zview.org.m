function [thedata,thegrid,han] = roms_zview(file,var,time,depth,grd,vec_d,uscale,varargin)
% [thedata,thegrid,han] = roms_zview(file,var,time,depth,grd,vec_d,uscale,varargin)
% 
% Inputs:
%
% file  = roms his/avg/rst etc nc file
% var   = variable to plot
% time  = time index into nc file
% depth = z depth of horizontal slice (m)
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
  error('Usage: roms_zview(file,var,time,depth,grd,vec_d,uscale,varargin)');
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

% check the grid information

% pcolor plot of the variable
switch var
  case { 'ubar','vbar','zeta','Hsbl','h','f','pm','pn',...
	'swrad','SST','dQdSST','shflux','swflux','SSS',...
	'sustr','svstr','Uwind','Vwind','Tair','Pair',...
	'sensible','latent'}
    [data,x,y,t,grd] = roms_2dslice(file,var,time,grd);
    depstr = [];
  otherwise    
    [data,x,y,t,grd] = roms_zslice(file,var,time,depth,grd);
    depstr = [' - Depth ' num2str(abs(depth)) ' m '];
end

switch var
  case { 'u','ubar','sustr'}
    mask = change(grd.mask_u,'==',0,NaN);
  case { 'v','vbar','svstr'}
    mask = change(grd.mask_v,'==',0,NaN);
  otherwise
    mask = grd.mask_rho_nan;
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
  try
    got_date = 0;
    basedate = parsetnc(nc{'dstart'}(:));
    tdate = julian(basedate)+t;
    got_date = 1;
  catch    
    disp([ ' Could not parse base date from ' file ' for ' var])
    got_date = 0;
  end
end
close(nc)

hanpc = pcolorjw(x,y,data.*mask);

caxis(cax);

if nargin > 5
  if vec_d
    % add vectors
    % ! sorry, this doesn't allow for {u,v}bar vectors on a 3d variable
    if depth 
      u = roms_zslice(file,'u',time,depth,grd);
      v = roms_zslice(file,'v',time,depth,grd);
      if isempty(depstr)
	depstr = [' - Vectors at depth ' num2str(abs(depth)) ' m '];
      end
    else  
      try 
        % a forcing file won't have u,v ...
        u = roms_2dslice(file,'ubar',time,grd);
        v = roms_2dslice(file,'vbar',time,grd);
      catch	
        % ... failing that look for wind stress
	% (should make this more general to allow for plotting
	% wind at 10m, but that would be on rho-points no u,v-points
	% so some extra checking is required -- maybe try to use
	% roms_addvect or roms_sview instead ... something for later)
	u = roms_2dslice(file,'sustr',time,grd);
	v = roms_2dslice(file,'svstr',time,grd);
	depstr = [ ' - Wind stress vectors '];
      end
      if isempty(depstr)
	depstr = [' - Vectors for depth-average velocity '];
      end
    end
    if nargin < 7
      uscale = 1;
    end    
    % hanquiver = roms_quiver(grd.lon_rho,grd.lat_rho,u,v,grd.angle,...
    %	vec_d,uscale,varargin{:});
    hanquiver = roms_quivergrd(u,v,grd,vec_d,uscale,varargin{:});
  end
end

amerc

try
  if findstr('leeuwin',grd.grd_file)
    gebco_eez(0,'k')
  elseif findstr('eauc',grd.grd_file)
    plotnzb
  elseif findstr('nena',grd.grd_file)
    plotnenacoast(3,'k')
  end
catch
end

if nargout > 0
  thedata.data = data;
  thedata.t = t;
  try
     nc = netcdf(file);
     thedata.base_date = nc{nc{'temp'}.time(:)}.units(:);
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
titlestr{2} = [upper(var) ' - Day ' num2str(t) depstr];
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

