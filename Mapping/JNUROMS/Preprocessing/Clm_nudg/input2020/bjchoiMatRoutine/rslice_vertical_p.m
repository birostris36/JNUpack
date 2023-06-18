function rslice_vertical(ncfile,varname,g)
% Inputs:
%
% ncfile   = roms his/avg/rst etc nc file
%          = 'his_playlow_05s_0121.nc';
% varname  = variable to plot
%          = 'salt' or 'temp' or 'u' or 'v' or 'speed'
% g        = grd structure (from roms_get_grid)
%
% Outputs:
% temperature and salinity vertical sections along the transect
% u and v are physical u and v.( u and v are not ROMS u and v)
% normalvel is the normal velocity across the transect.



% grid information

if exist('g')~=1
%   g = grd('hudson');
    grd_file = '/home/bchoi/roms/hudson/grid/roms_latte_grid_3b.nc';
    scoord = [5 0.4 50 30];
    g = roms_get_grid(grd_file,scoord);
end

  grid.x=squeeze(g.lon_rho);
  grid.y=squeeze(g.lat_rho);
  grid.N=g.N;
  grid.z_r=g.z_r;
  grid.mask=g.mask_rho;
  grid.mask_nan=g.mask_rho_nan;

  var=varname;

% read data and store them 
 
switch varname
 case{ 'temp', 'salt'}
  nc=netcdf(ncfile,'read');
  gdata=nc{varname}(:);
  close(nc)
 case{ 'u', 'v', 'speed'}
  nc=netcdf(ncfile,'read');
  uroms=nc{'u'}(:);
  vroms=nc{'v'}(:);
  close(nc)
 otherwise
  error(' error in selecting variable name')
end

% coordinate rotation from ROMS u v to physical space u v

switch varname
 case{ 'u', 'v', 'speed' }

  [N M L]=size(uroms);
  uroms_e = ones(N,M,L+2)*NaN;
  uroms_e(:,:,1)   = uroms(:,:,1);
  uroms_e(:,:,end) = uroms(:,:,end);
  uroms_e(:,:,2:end-1) = uroms;

  [N M L]=size(vroms);
  vroms_e = ones(N,M+2,L)*NaN;
  vroms_e(:,1,:)   = vroms(:,1,:);
  vroms_e(:,end,:) = vroms(:,end,:) ;
  vroms_e(:,2:end-1,:) = vroms;

  uroms = ( uroms_e(:,:,1:end-1) + uroms_e(:,:,2:end) )*0.5;
  vroms = ( vroms_e(:,1:end-1,:) + vroms_e(:,2:end,:) )*0.5;           
  cosmat=cos(g.angle);
  sinmat=sin(g.angle);
  cosz=repmat(cosmat,[1 1 g.N]);
  sinz=repmat(sinmat,[1 1 g.N]);
  cosz=permute(cosz,[3 1 2]);
  sinz=permute(sinz,[3 1 2]);
  uphysical = uroms.*cosz - vroms.*sinz;
  vphysical = uroms.*sinz + vroms.*cosz;

end

switch varname
case{ 'u' }
  gdata=uphysical;
case{ 'v' }
  gdata=vphysical;
case{ 'speed' }
  speed=sqrt( uphysical.*uphysical + vphysical.*vphysical );
  gdata=speed;
end


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
    f = figure;
    figure(f)
    set(f,'Position',[300 300 900 600]);
    pcolorjw(grid.x,grid.y,squeeze(gdata(grid.N,:,:)).*grid.mask_nan);
    colorbar
    title([varname ' at the surface level '])
    endpt_lon(1) = input('Starting Lon = ');
    endpt_lat(1) = input('Starting Lat = ');
    endpt_lon(2) = input('  Ending Lon = ');
    endpt_lat(2) = input('  Ending Lat = ');
end;
   figure(f)
   hold on 
   line(endpt_lon,endpt_lat)
   disp([' You chose a transect from ' ...
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
% pcolor ( x/1000, z, z_yx_data ); shading interp;
  colorbar;
  xlabel ( 'kilometers' );
  ylabel ( 'meters' );
  title( varname );

switch varname 
  case{ 'speed' }

     tangential_vector=[pt2.x-pt1.x pt2.y-pt2.y];
     tan_unit_vec=tangential_vector/norm(tangential_vector);
%    tangential_vector=[m_lldist([pt1.x pt2.x],[pt1.y pt1.y]), ...
%         	       m_lldist([pt1.x pt1.x],[pt1.y pt2.y])]
%    tan_unit_vec=tangential_vector/m_lldist([pt1.x pt2.x],[pt1.y pt2.y])
    nor_unit_vec=[-tan_unit_vec(2) tan_unit_vec(1)];
    gdata=uphysical.*nor_unit_vec(1) + vphysical.*nor_unit_vec(2);

    [x,z,z_yx_data,dist]=rslice_vertical_interp1(pt1,pt2,grid,gdata);

    figure
    pcolor ( x/1000, z, z_yx_data ); shading flat;
    colorbar;
    xlabel ( 'kilometers' );
    ylabel ( 'meters' );
    title( 'normal speed: into the page +, out of page -' );

end
