function [volume,heat,salinity] = ...
          ij_transport(outputfile,grid_nick_name);

% calculate transports of volume (in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
% 
% across a line segment (slice along a constant I or J)
%
% 'outputfile' data in file name history or average file
% 'grid_nick_name' = grid information such as 'eas' 'hudson' 'latte'
%
% vertical coordinate changes in time = h + zeta(t)
%

%   grid information

    g = grd(grid_nick_name);
   %g = grd('hudson');
   %g = grd('latte');
   %g = grd('eas');

    N = g.N;  
    hz=g.z_w(N+1:2,:,:)-g.z_w(N:1,:,:); % z_w: [31x142x254]

%   choose two grid point (i,j)

    disp(' Enter two grid points (i,j): ')

    endpt_lon(1) = input('Starting i = ');
    endpt_lat(1) = input('Starting j = ');
    endpt_lon(2) = input('  Ending i = ');
    endpt_lat(2) = input('  Ending j = ');

%   load model output 

    nc=netcdf(outputfile,'read');
    zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
    u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
    v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
    temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
    salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
    close(nc)

%

    h_total = g.h + squeeze(zeta(1,:,:)); % total water thickness
    for level=1:N
     Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./h);
    end

%

  if ( endpt_lon(1) == endpt_lon(2) )

    % average grd.z_r to Arakawa-C u points
    % grd.z_r = (level, j, i)
    zM = size(grd.z_r,2);
    zMm = zM-1;
    zL = size(grd.z_r,3);
    zLm = zL-1;
    
    zh = 0.5*(grd.z_w(1,:,1:zLm)+grd.z_w(1,:,2:zL));
    z = 0.5*(grd.z_r(:,:,1:zLm)+grd.z_r(:,:,2:zL));
    z0 = 0.5*(grd.z_w(Np,:,1:zLm)+grd.z_w(Np,:,2:zL));
    
    zw = 0.5*(grd.z_w(:,:,1:zLm)+grd.z_w(:,:,2:zL));
    pm = 0.5*(grd.pm(:,1:(end-1))+grd.pm(:,2:end));
    pn = 0.5*(grd.pn(:,1:(end-1))+grd.pn(:,2:end));
    h = 0.5*(grd.h(:,1:(end-1))+grd.h(:,2:end));
    
    x = grd.lon_u;
    y = grd.lat_u;
    mask = grd.mask_u;
    
    dx = 1./pm;
    dy = 1./pn;
    dz = diff(zw,1);

    u = getnc(file,'u');
    u = squeeze(u(:,:,index));
    
    temp = getnc(file,'temp');
    temp = 0.5*(squeeze(temp(:,:,index)) + squeeze(temp(:,:,index+1)));
    
    salt = getnc(file,'salt');
    salt = 0.5*(squeeze(salt(:,:,index)) + squeeze(salt(:,:,index+1)));;
    salt = salt./1e3;
    
    h = h(:,index);
    x = x(:,index);
    y = y(:,index);
    m = mask(:,index);
    z = z(:,:,index);
    dx = dx(:,index);
    dy = dy(:,index);
    dz = dz(:,:,index);
    
     trans = u;
    strans = u.*salt;
    ttrans = u.*temp;

  else if  ( endpt_lat(1) == endpt_lat(2) )

    % average grd.z_r to Arakawa-C v points
    zM = size(grd.z_r,2);
    zMm = zM-1;
    zL = size(grd.z_r,3);
    zLm = zL-1;
    zh = 0.5*(grd.z_w(1,1:zMm,:)+grd.z_w(1,2:zM,:));
    z = 0.5*(grd.z_r(:,1:zMm,:)+grd.z_r(:,2:zM,:));
    z0 = 0.5*(grd.z_w(Np,1:zMm,:)+grd.z_w(Np,2:zM,:));
    
    zw = 0.5*(grd.z_w(:,1:zMm,:)+grd.z_w(:,2:zM,:));
    pm = 0.5*(grd.pm(1:(end-1),:)+grd.pm(2:end,:));
    pn = 0.5*(grd.pn(1:(end-1),:)+grd.pn(2:end,:));
    h = 0.5*(grd.h(1:(end-1),:)+grd.h(2:end,:));
    
    x = grd.lon_v;
    y = grd.lat_v;
    mask = grd.mask_v; 
    
    dx = 1./pm;
    dy = 1./pn;
    dz = diff(zw,1);

    v = getnc(file,'v');
    v = squeeze(v(:,index,:));
    
    temp = getnc(file,'temp');
    temp = 0.5*(squeeze(temp(:,index,:)) + squeeze(temp(:,index+1,:)));
    
    salt = getnc(file,'salt');
    salt = 0.5*(squeeze(salt(:,index,:)) + squeeze(salt(:,index+1,:)));;
    salt = salt./1e3;
    
    fw = (sref-salt)./sref;
        
    h = h(index,:);
    x = x(index,:);
    y = y(index,:);
    m = mask(index,:);
    z = z(:,index,:);
    dx = dx(index,:);
    dy = dy(index,:);
    dz = dz(:,index,:);
    
     trans = v;
    strans = v.*salt;
    ttrans = v.*temp;

  else
  
      error(' either i or j must be same ! ')
   
  end

           
z = squeeze(z);
dz = squeeze(dz);

lon = x(:);
lat = y(:);
s = cumsum([0; sw_dist(lat,lon,'km')]);

x = repmat(x(:)',[grd.N 1]);
y = repmat(y(:)',[grd.N 1]);
m = repmat(m(:)',[grd.N 1]);
dx = repmat(dx(:)',[grd.N 1]);
dy = repmat(dy(:)',[grd.N 1]);

switch direction
  case {'xi','i','u'}
       volume = sumnan(trans.*dy.*dz)./1e6;
         heat = sumnan(ttrans.*dy.*dz).*4.2e-9;
     salinity = sumnan(strans.*dy.*dz).*1e-6;

   disp(['At index ' num2str(index) ' in the ' direction ' direction : '])
   disp(['     Volume transport is ' num2str(volume) ' Sv']);
   disp(['       Heat transport is ' num2str(heat) ' PW']);
   disp(['       Salt transport is ' num2str(salinity) ' x 10^9 Kg']);
   disp([' Freshwater transport is ' num2str(freshwater) ' mSv']);
   disp(' ');
   
  case {'eta','j','v'}
       volume = sumnan(trans.*dx.*dz)./1e6;
         heat = sumnan(ttrans.*dx.*dz).*4.2e-9;
     salinity = sumnan(strans.*dx.*dz).*1e-6;
   freshwater = sumnan(fwtrans.*dx.*dz)./1e3;

   disp(['At index ' num2str(index) ' in the ' direction ' direction : '])
   disp(['     Volume transport is ' num2str(volume) ' Sv']);
   disp(['       Heat transport is ' num2str(heat) ' PW']);
   disp(['       Salt transport is ' num2str(salinity) ' x 10^9 Kg']);
   disp([' Freshwater transport is ' num2str(freshwater) ' mSv']);
   disp(' ');
   
end;
        
