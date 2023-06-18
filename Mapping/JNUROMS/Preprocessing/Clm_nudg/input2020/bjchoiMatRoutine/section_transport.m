function [volume,heat,salinity,freshwater] = section_transport(file,direction,index,grd,sref);

%function [volume,heat,salinity,freshwater] = section_transport(file,direction,index,grd,sref);
% calculate transports of volume (in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
%                    fresh water (in mSv, 10^3 m^3/sec, relative to sref) 
% 
% across the section in a given direction
%
% [volume,heat,salinity,freshwater] =  cross_section_transports(file,direction,grd,sref);
%
% data in file name 'file'
% direction = (slice along a constant XI or ETA or I or J)
% grd = grid information 
% sref = reference salinity
%

status = 0;
   disp(' ')

if nargin < 5;
   sref = 0.0348;
   disp('Fresh Water transport calculated relative to 34.8 psu')
   disp(' ')
end;

if sref > 1;
   sref = sref / 1000;
end;

if nargin < 4;
    disp('Getting grid info')
    grd = roms_get_grid(file,file);
    disp(' ')
end;
Np = grd.N+1;

direction = lower(direction);

switch direction
  case {'xi','i','u'}    
    % average grd.z_r to Arakawa-C u points
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
    
    fw = (sref-salt)./sref;
        
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
   fwtrans = u.*fw;

  case {'eta','j','v'}
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
   fwtrans = v.*fw;

  otherwise   
  
  error('NOT CONFIGURED FOR Z DIRECTION YET')
   
    % for temp, salt, rho, w
    zh = grd.z_w(1,:,:);
    z = grd.z_r;
    z0 = grd.z_w(Np,:,:);
    x = grd.lon_rho;
    y = grd.lat_rho;
    mask = grd.mask_rho_nan;

    zw = grd.z_w;
    pm = grd.pm;
    pn = grd.pn;
    h = grd.h;
    
    % might have to do something special to handle omega
    
end
   
%keyboard

           
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
   freshwater = sumnan(fwtrans.*dy.*dz)./1e3;

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
        
