function  [x,z,z_yx_data,izeta,dist]=rslice_vertical_p2dataonly(ncfile,varname,itime,g,point1,point2);
%
% [x,z,z_yx_data,izeta,dist]=rslice_vertical_p2(ncfile,varname,g)
% Inputs:
%
% ncfile   = roms his/avg/rst etc nc file
%          = 'his_playlow_05s_0121.nc';
% varname  = variable to plot
%          = 'salt' or 'temp' or 'u' or 'v' or 'speedn', or 'speedt'
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

% grid information

if exist('g')~=1
%  g = grd('hudson');
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
  if ( length(varname) == 6 )
  var='u';
  end

% read data and store them 
 
  nc=netcdf(ncfile,'read');
  t = nc{nc{var}.time(:)}(1);
  time_units = nc{nc{var}.time(:)}.units(:);
  if strcmp(time_units(1:3),'sec')
    t = t/86400;
  end

  time=nc{'ocean_time'}(:);
  zeta=nc{'zeta'}(:);
  ntime=length(time);
  if( ntime > 1 ) 
     zeta=squeeze( zeta(itime,:,:,:) );
  end 

  switch varname
   case{ 'temp', 'salt'}
    gdata=nc{varname}(:);
    if( ntime > 1 ) 
     gdata=squeeze( gdata(itime,:,:,:) );
    end
   case{ 'u', 'v', 'speedn', 'speedt'}
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
 case{ 'u', 'v', 'speedn', 'speedt' }
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

  uroms_rho = ( uroms_e(:,:,1:end-1) + uroms_e(:,:,2:end) )*0.5;
  vroms_rho = ( vroms_e(:,1:end-1,:) + vroms_e(:,2:end,:) )*0.5;           
  cosmat=cos(g.angle);
  sinmat=sin(g.angle);
  cosz=repmat(cosmat,[1 1 g.N]);
  sinz=repmat(sinmat,[1 1 g.N]);
  cosz=permute(cosz,[3 1 2]);
  sinz=permute(sinz,[3 1 2]);
  uphysical = uroms_rho.*cosz - vroms_rho.*sinz;
  vphysical = uroms_rho.*sinz + vroms_rho.*cosz;
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


  pt1.x = point1(1);
  pt1.y = point1(2);
  pt2.x = point2(1);
  pt2.y = point2(2);

switch varname 
  case{ 'speedn', 'speedt'  }
    tangential_vector=[pt2.x-pt1.x pt2.y-pt1.y];
    tan_unit_vec=tangential_vector/norm(tangential_vector);
    nor_unit_vec=[-tan_unit_vec(2) tan_unit_vec(1)];
    if (varname=='speedn')
     gdata=uphysical.*nor_unit_vec(1) + vphysical.*nor_unit_vec(2);  % normal velocity
    elseif (varname=='speedt') 
     gdata=uphysical.*tan_unit_vec(1) + vphysical.*tan_unit_vec(2);  % tangential velocity
    end
end

  % interpolation method: 'bilinear'
  [x,z,z_yx_data,izeta,dist]=rslice_vertical_interp1z(pt1,pt2,grid,gdata,zeta);

  save slice_data.mat x z z_yx_data izeta varname;
