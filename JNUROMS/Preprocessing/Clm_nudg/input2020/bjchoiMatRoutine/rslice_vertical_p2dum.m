function  [x,z,z_yx_data,izeta,dist]=rslice_vertical_p2dum(ncfile,varname,itime,g,point1,point2)
%
% rslice_vertical_p2(ncfile,varname,g)
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
%
%
% transect information
% point1 = [lon1 lat1]  starting point
% point2 = [lon2 lat2]  ending   point
% =============================================================

%
% rslice_vertical_p2dum.m ==> dum use x-axis from dist=x_transect.
%



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
  grid.sc_r=g.sc_r;
  grid.Cs_r=g.Cs_r;
  grid.hc=g.hc;
  grid.h=g.h;

  var=varname;
  if ( length(varname) == 5 )
  var='u';
  end

% read data and store them 
 
  nc=netcdf(ncfile,'read');
  t = nc{nc{var}.time(:)}(itime);
  time_units = nc{nc{var}.time(:)}.units(:);
  if strcmp(time_units(1:3),'sec')
    t = t/86400;
  end

  time=nc{'ocean_time'}(:);
  zeta=nc{'zeta'}(:);
  ntime=length(time);
  if( ntime > 1 ) 
     zeta=squeeze( zeta(itime,:,:,:) );
     disp([' Display the a time frame '])
     disp([' at ' num2str(t) ' day']) 
  end 

  switch varname
   case{ 'temp', 'salt'}
    gdata=nc{varname}(:);
    if( ntime > 1 ) 
     gdata=squeeze( gdata(itime,:,:,:) );
    end
   case{ 'u', 'v', 'speed'}
    uroms=nc{'u'}(:);
    vroms=nc{'v'}(:);
    if( ntime > 1 ) 
     uroms=squeeze( uroms(itime,:,:,:) );
     vroms=squeeze( vroms(itime,:,:,:) );
    end
   otherwise
    error(' error in selecting variable name')
  end
  close(nc)

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



if ( nargin < 5) 

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
     titlestr{1} = ['file: ' strrep_(ncfile)];
     titlestr{2} = [upper(varname) ' - Day ' num2str(t) ' at Nth level'];
     title(titlestr)
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


else

% transection information

  pt1.x = point1(1);
  pt1.y = point1(2);
  pt2.x = point2(1);
  pt2.y = point2(2);

end


% interpolation method: 'nearest_neighbor'
% [x,z,z_yx_data,izeta,dist]=rslice_vertical_interp0z(pt1,pt2,grid,gdata,zeta);
% interpolation method: 'bilinear'
  [x,z,z_yx_data,izeta,dist]=rslice_vertical_interp1z(pt1,pt2,grid,gdata,zeta);

if (nargin < 5 | length(varname) < 5)
  if(nargin < 5)
  figure
  end
% pcolor ( x/1000, z, z_yx_data ); shading flat;
% pcolor ( x/1000, z, z_yx_data ); shading interp;
  pcolor ( dist, z, z_yx_data ); shading interp;
  hold on
% plot(x/1000,izeta','k')
  plot(dist,izeta','k')
  colorbar;
  xlabel ( 'kilometers' );
  ylabel ( 'meters' );
    titlestr{1} = ['file: ' strrep_(ncfile)];
    titlestr{2} = [upper(varname) ' - Day ' num2str(t)];
    title(titlestr)
end


switch varname 
  case{ 'speed' }

     tangential_vector=[pt2.x-pt1.x pt2.y-pt1.y];
     tan_unit_vec=tangential_vector/norm(tangential_vector);
%    tangential_vector=[m_lldist([pt1.x pt2.x],[pt1.y pt1.y]), ...
%         	       m_lldist([pt1.x pt1.x],[pt1.y pt2.y])]
%    tan_unit_vec=tangential_vector/m_lldist([pt1.x pt2.x],[pt1.y pt2.y])
    nor_unit_vec=[-tan_unit_vec(2) tan_unit_vec(1)];
    gdata=uphysical.*nor_unit_vec(1) + vphysical.*nor_unit_vec(2);

% interpolation method: 'nearest_neighbor'
% [x,z,z_yx_data,izeta,dist]=rslice_vertical_interp0z(pt1,pt2,grid,gdata,zeta);
% interpolation method: 'bilinear'
  [x,z,z_yx_data,izeta,dist]=rslice_vertical_interp1z(pt1,pt2,grid,gdata,zeta);

    if (nargin < 5)
    figure
    end
%   pcolor ( x/1000, z, z_yx_data ); shading flat;
    pcolor ( x/1000, z, z_yx_data ); shading interp;
    hold on
    plot(x/1000,izeta','k')
    colorbar;
    xlabel ( 'kilometers' );
    ylabel ( 'meters' );
    titlestr{1} = ['file: ' strrep_(ncfile)];
    titlestr{2} = [upper(varname) ' - Day ' num2str(t)];
    title(titlestr)
%   title( 'normal speed: into the page +, out of page -' );
end

  save slice_data.mat x z z_yx_data izeta varname;
