function rslice_vertical(ncfile,varname,g)
% Inputs:
%
% ncfile   = roms his/avg/rst etc nc file
% varname  = variable to plot% grd can be 
% g        = grd structure (from roms_get_grid)
%
% Outputs:
% temperature and salinity vertical sections along the transect
% ROMS u and v along the transect ( u and v are not physical u and v)
%



% grid information

if exist('g')~=1
%   g = grd('hudson');
    grd_file = '/home/bchoi/roms/hudson/grid/roms_latte_grid_3b.nc';
    scoord = [5 0.4 50 30];
    g = roms_get_grid(grd_file,scoord);
end

switch varname
 case { 'v' }
  grid.x=squeeze(g.lon_v);
  grid.y=squeeze(g.lat_v);
  grid.N=g.N;
  grid.z_r=g.z_r;
  grid.mask=g.mask_v;
  grid.mask_nan=grid.mask;
  Inan=find( grid.mask_nan < 1 );
  grid.mask_nan(Inan)=g.mask_v(Inan)*NaN;
 case { 'u' }
  grid.x=squeeze(g.lon_u);
  grid.y=squeeze(g.lat_u);
  grid.N=g.N;
  grid.z_r=g.z_r;
  grid.mask=g.mask_u;
  grid.mask_nan=grid.mask;
  Inan=find( grid.mask_nan < 1 );
  grid.mask_nan(Inan)=g.mask_u(Inan)*NaN;
 case { 'salt','temp'}
  grid.x=squeeze(g.lon_rho);
  grid.y=squeeze(g.lat_rho);
  grid.N=g.N;
  grid.z_r=g.z_r;
  grid.mask=g.mask_rho;
  grid.mask_nan=g.mask_rho_nan;
end
  
  var=varname;

% read data and store them 
 
  nc=netcdf(ncfile,'read');
  gdata=nc{varname}(:);
  close(nc)

%  graphical input or text input

  qq = input('Do you want to use the keyboard or the mouse? ','s');
  if isempty(qq)
    qq = 'm';
  end
if strcmpi(qq(1),'m');
  disp(' ');
  f = input('Which Figure? (Return to create new figure) ');
  if isempty(f)
     f = figure;
     figure(f)
     set(f,'Position',[300 300 900 600]);
     pcolorjw(grid.x,grid.y,squeeze(gdata(grid.N,:,:)).*grid.mask_nan);
     colorbar
     title([varname ' at the surface level '])
  else;
     figure(f)
  end;
  for n=1:2;
    figure(f)
    [endpt_lon(n),endpt_lat(n)]=ginput(1);
  end;
else;
    endpt_lon(1) = input('Starting Lon = ');
    endpt_lat(1) = input('Starting Lat = ');
    endpt_lon(2) = input('  Ending Lon = ');
    endpt_lat(2) = input('  Ending Lat = ');
end;
   figure(f)
   hold on 
   line(endpt_lon,endpt_lat)
   disp([' You choose a transect from ' ...
           num2str(endpt_lon(1)) ' ' num2str(endpt_lat(1)) ])
   disp(['                       to   ' ... 
           num2str(endpt_lon(2)) ' ' num2str(endpt_lat(2)) ])


% transection information

  pt1.x = endpt_lon(1);
  pt1.y = endpt_lat(1);
  pt2.x = endpt_lon(2);
  pt2.y = endpt_lat(2);


% interpolation method: 'nearest_neighbor'
% [x,z,z_yx_data,dist]=rslice_vertical_interp0(pt1,pt2,grid,gdata);
% interpolation method: 'bilinear'
  [x,z,z_yx_data,dist]=rslice_vertical_interp1(pt1,pt2,grid,gdata);

  save slice_data.mat x z z_yx_data varname;

  figure
  pcolor ( x/1000, z, z_yx_data ); shading flat;
  colorbar;
  xlabel ( 'kilometers' );
  ylabel ( 'meters' );



