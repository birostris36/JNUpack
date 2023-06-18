function [transport,los,angle0]=roms_section(outputfile,gridfile);

%   [transport,los,angle0]=roms_transect(outputfile,gridfile);
%
%     This will calculate the net transport across a user-defined 
%       trasnsect (entered via mouse or keyboard) from a 
%       ROMS output (avg or his) file.
%
%     file is the complete path to the output file
%      grd is the grid information (optional if NENA or LATTE run)
%
%        los is the length of the transect in KM
%     angle0 is the direction of transect from TRUE EAST
%  transport is the net flow across the transect in Sv 
%
%    transect  _ transport
%          \   /|
%           \ /
%            \
%             \ 

addpath  /home/bchoi/matlab/netcdf_csiro

% degree to radian

d2r = pi/180;
r2d = 180/pi;

%get grid info 

      nc = netcdf(outputfile);  
      theta_s = nc{'theta_s'}(:);
      theta_b = nc{'theta_b'}(:);
      Tcline = nc{'Tcline'}(:);
      N = length(nc{'sc_r'}(:));
      Np = N+1;
      time = nc{'ocean_time'}(:);
      units = nc{'ocean_time'}.units(:);
      if strcmp(units(1:3),'sec')
         time = time/86400;
      end
      jtime = time+2448623; 
      result=close(nc);
     
      grd = roms_get_grid(gridfile,[theta_s theta_b Tcline N]);

% get spatial coordinates

     zhr = grd.z_w(1,:,:);
      zr = grd.z_r;
     z0r = grd.z_w(Np,:,:);
      xr = grd.lon_rho;
      yr = grd.lat_rho;
    angr = grd.angle;
   maskr = grd.mask_rho_nan;

     zwr = grd.z_w;
     pmr = grd.pm;
     pnr = grd.pn;
      hr = grd.h;

      zM = size(grd.z_r,2);
     zMm = zM-1;
      zL = size(grd.z_r,3);
     zLm = zL-1;

% get free surface info
%
    fsr = getnc(file,'zeta');

% get velocity data
%
%   first U
%
     zhu = 0.5*(grd.z_w(1,:,1:zLm)+grd.z_w(1,:,2:zL));
      zu = 0.5*(grd.z_r(:,:,1:zLm)+grd.z_r(:,:,2:zL));
     z0u = 0.5*(grd.z_w(Np,:,1:zLm)+grd.z_w(Np,:,2:zL));
    
     zwu = 0.5*(grd.z_w(:,:,1:zLm)+grd.z_w(:,:,2:zL));
     pmu = 0.5*(grd.pm(:,1:(end-1))+grd.pm(:,2:end));
     pnu = 0.5*(grd.pn(:,1:(end-1))+grd.pn(:,2:end));
      hu = 0.5*(grd.h(:,1:(end-1))+grd.h(:,2:end));
    angu = 0.5*(grd.angle(:,1:(end-1))+grd.angle(:,2:end));
    
      xu = grd.lon_u;
      yu = grd.lat_u;
 warning off;
   masku = grd.mask_u./grd.mask_u;
 warning on;
      
     dxu = 1./pmu;
     dyu = 1./pnu;
     dzu = diff(zwu,1);

       u = getnc(file,'u');
     fsu = 0.5*(fsr(:,1:(end-1))+fsr(:,2:end));

%
%   then V
%

     zhv = 0.5*(grd.z_w(1,1:zMm,:)+grd.z_w(1,2:zM,:));
      zv = 0.5*(grd.z_r(:,1:zMm,:)+grd.z_r(:,2:zM,:));
     z0v = 0.5*(grd.z_w(Np,1:zMm,:)+grd.z_w(Np,2:zM,:));
    
     zwv = 0.5*(grd.z_w(:,1:zMm,:)+grd.z_w(:,2:zM,:));
     pmv = 0.5*(grd.pm(1:(end-1),:)+grd.pm(2:end,:));
     pnv = 0.5*(grd.pn(1:(end-1),:)+grd.pn(2:end,:));
      hv = 0.5*(grd.h(1:(end-1),:)+grd.h(2:end,:));
    angv = 0.5*(grd.angle(1:(end-1),:)+grd.angle(2:end,:));
    
      xv = grd.lon_v;
      yv = grd.lat_v;
 warning off;
   maskv = grd.mask_v./grd.mask_v;
 warning on;
    
     dxv = 1./pmv;
     dyv = 1./pnv;
     dzv = diff(zwv,1);

       v = getnc(file,'v');
     fsv = 0.5*(fsr(1:(end-1),:)+fsr(2:end,:));

% Integrate U and V with depth

% velocity
   udz = u.*dzu;
   vdz = v.*dzv;
   
   udz(20,:,:) = squeeze(u(20,:,:)).*(fsu + squeeze(dzu(20,:,:)));
   vdz(20,:,:) = squeeze(v(20,:,:)).*(fsv + squeeze(dzv(20,:,:)));

    uz = squeeze(sum(udz,1));
   uzr = av2(uz);
    
    vz = squeeze(sum(vdz,1));
   vzr = av2(vz')';
    
%
%   UZ and VZ in m^2/sec
%

  xuv = av2(xu);
  yuv = av2(yu);
anguv = av2(av2(angr)')';

if 0;
cax=[-1 1];
f = figure;
figure(f)
set(f,'Position',[300 300 900 600]);
pcolorjw(xr,yr,fsr.*maskr);caxis(cax);colorbar
hold on 
roms_quiver(xuv,yuv,uz.*masku,vz.*maskv,anguv,3,0.001,'k')
hold off
pnc;
title(['Free Surface Height on ' greg2str(jtime)],'fontsize',14);
end;

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
     pcolorjw(grd.lon_rho,grd.lat_rho,grd.mask_rho_nan);
%    pnc;
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

hold on;
h=plot(endpt_lon,endpt_lat,'k');
set(h,'LineWidth',2);
hold off;
drawnow

los_in_m=m_lldist(endpt_lon,endpt_lat);
los_in_km = los_in_m/1e3;
myangle    = angle(diff(endpt_lon)*d2r + sqrt(-1)*diff(endpt_lat)*d2r);
myangledeg = myangle*r2d;

if (myangle >= pi/2)
     myangle = myangle-pi;
  myangledeg = myangledeg-180;
   endpt_lon = fliplr(endpt_lon);
   endpt_lat = fliplr(endpt_lat);
elseif (myangle < -pi/2)
     myangle = myangle+pi;
  myangledeg = myangledeg+180;
   endpt_lon = fliplr(endpt_lon);
   endpt_lat = fliplr(endpt_lat);
end

disp(['Length of Section is ' num2str(los_in_km) ' km']);
disp(['   Rotation Angle is ' num2str(myangledeg) ' degrees from EAST']);

% DEFINE NEW COORDINATES
%
%   break the transect into ~km sized sections
%

%  check for true E-W or N-S section
%
check = diff(endpt_lon)*diff(endpt_lat);

if (check ~= 0)

lonq = [endpt_lon(1):(diff(endpt_lon)/floor(los_in_km)):endpt_lon(2)];
latq = [endpt_lat(1):(diff(endpt_lat)/floor(los_in_km)):endpt_lat(2)];

elseif (diff(endpt_lon) == 0)

latq = [endpt_lat(1):(diff(endpt_lat)/floor(los_in_km)):endpt_lat(2)];
lonq = ones(size(latq))*endpt_lon(1);

else

lonq = [endpt_lon(1):(diff(endpt_lon)/floor(los_in_km)):endpt_lon(2)];
latq = ones(size(lonq))*endpt_lat(1);

end;

lons = meshgrid(lonq',latq);
lats = meshgrid(latq,lonq')';

%  
% 
eitheta = exp(-sqrt(-1)*(myangle+mean(mean(anguv))));

disp('Regridding U')
uz2 = griddata(xuv,yuv,uzr,lonq,latq);
disp('Regridding V')
vz2 = griddata(xuv,yuv,vzr,lonq,latq);

disp('Calculating Along and Across Transect Transports')
U = (uz2+sqrt(-1)*vz2).*eitheta;

along = imag(U);
across = real(U);

totalacross = sumnan(av2(across).*m_lldist(lonq,latq))/1e6;

disp (['Total transport ACROSS the Section is ' num2str(round(totalacross*100)/100) ' Sv'])
disp (' ')

transport = totalacross;
los = los_in_km;
angle0 = myangledeg;

