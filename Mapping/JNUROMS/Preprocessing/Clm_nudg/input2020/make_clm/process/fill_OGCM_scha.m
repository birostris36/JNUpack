         function fill_OGCM_scha(OGCM_dir,OGCM_prefix,year,month,day,Roa,interp_method,...
                     lonU,latU,lonV,latV,lonT,latT,Z,tin,...
		     nc_clm,nc_bry,lon,lat,angle,h,tout,obc,vtransform,vstretching)
         
%          year = 2017;
%          month = 01;
%          day = 01;
%          tin = 1;
%          nc_bry = [];
%          tout = 1;
%         interp_OGCM_scha(OGCM_dir,OGCM_prefix,Time(1),Time(2),Time(3),Roa,interp_method,...
%         lonU,latU,lonV,latV,lonT,latT,Z,1,...
%         [],nc_clm,lon,lat,angle,h,ii,obc,vtransform,vstretching)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Read the local OGCM files and perform the interpolations
% Adapted for reducing computationnal time for bry file
% by S. Illig, IRD-LEGOS
%
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr 
%
%  Updated    6-Sep-2006 by Pierrick Penven : Nothing special for the bry file 
%  Updated    5-Nov-2006 by Pierrick Penven : A bit of cleaning...
%  Update    13-Sep-2009 by Gildas Cambon :   Begin treatments case  for the bry
%  Update    March-2011 by Serena Illig : optimisation bry
%  Update    Nov-2011 by Pierrick Penven : cleaning
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
conserv=1; % same barotropic velocities as the OGCM
%
disp(['  Horizontal interpolation: ',...
      OGCM_prefix,num2str(year),...
      num2str(month,'%2.2i'),num2str(day,'%2.2i'),'.cdf'])
%
%
% ROMS grid angle
%
cosa=cos(angle);
sina=sin(angle);
%
% Open the OGCM file
%
nc=netcdf([OGCM_dir,OGCM_prefix,num2str(year),...
    num2str(month,'%2.2i'),num2str(day,'%2.2i'),'.cdf']);
%
% Interpole data on the OGCM Z grid and ROMS horizontal grid
%
% Get zeta because it is needed to compute vertical levels of ROMS grid
zeta=ext_data_OGCM(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method);

%
if ~isempty(nc_clm)
%
% Read and extrapole the 2D variables
%
%   u2d=ext_data_OGCM(nc,lonU,latU,'ubar',tin,lon,lat,1,Roa,interp_method);
%   v2d=ext_data_OGCM(nc,lonV,latV,'vbar',tin,lon,lat,1,Roa,interp_method);
%   ubar=rho2u_2d(u2d.*cosa+v2d.*sina);
%   vbar=rho2v_2d(v2d.*cosa-u2d.*sina);
%
% Read and extrapole the 3D variables
% 
  NZ=length(Z);
  [M,L]=size(lon);
  dz=gradient(Z);
  temp=zeros(NZ,M,L);
  salt=zeros(NZ,M,L);
%   u=zeros(NZ,M,L-1);
%   v=zeros(NZ,M-1,L);
  for k=1:NZ
    if rem(k,10)==0
      disp(['  Level ',num2str(k),' of ',num2str(NZ)])
    end
%     u2d=ext_data_OGCM(nc,lonU,latU,'u',tin,lon,lat,k,Roa,interp_method);
%     v2d=ext_data_OGCM(nc,lonV,latV,'v',tin,lon,lat,k,Roa,interp_method);
%     u(k,:,:)=rho2u_2d(u2d.*cosa+v2d.*sina);
%     v(k,:,:)=rho2v_2d(v2d.*cosa-u2d.*sina);
    temp(k,:,:)=ext_data_OGCM(nc,lonT,latT,'temp',tin,lon,lat,k,Roa,interp_method);
    salt(k,:,:)=ext_data_OGCM(nc,lonT,latT,'salt',tin,lon,lat,k,Roa,interp_method);
  end
end


%
% Close the OGCM file
%
close(nc)
%
%
% Get the ROMS vertical grid
%
disp('  Vertical interpolations')
if ~isempty(nc_clm)
  theta_s=nc_clm{'theta_s'}(:);
  theta_b=nc_clm{'theta_b'}(:);
  hc=nc_clm{'hc'}(:);
  N=length(nc_clm('s_rho'));
end
% if ~isempty(nc_bry)
%   theta_s=nc_bry{'theta_s'}(:);
%   theta_b=nc_bry{'theta_b'}(:);
%   hc=nc_bry{'hc'}(:);
%   N=length(nc_bry('s_rho'));
% end
%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
Z=[100;Z;-100000];
%
% ROMS vertical grid
%
zr=zlevs(vtransform,vstretching,theta_s,theta_b,hc,N,1,h,zeta,1);
zu=rho2u_3d(zr);
zv=rho2v_3d(zr);
zw=zlevs(vtransform,vstretching,theta_s,theta_b,hc,N,5,h,zeta,1);
dzr=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=rho2u_3d(dzr);
dzv=rho2v_3d(dzr);
%
%
% Vertical interpolation in case of clim file
%
%

if ~isempty(nc_clm)
%
% Add a level on top and bottom with no-gradient
%

%   u=cat(1,u(1,:,:),u);
%   u=cat(1,u,u(end,:,:));
%   v=cat(1,v(1,:,:),v);
%   v=cat(1,v,v(end,:,:));
  temp=cat(1,temp(1,:,:),temp);
  temp=cat(1,temp,temp(end,:,:));
  salt=cat(1,salt,salt(end,:,:));
  salt=cat(1,salt(1,:,:),salt);
% 
% Perform the vertical interpolations 
%

%   u=ztosigma(flipdim(u,1),zu,flipud(Z));
%   v=ztosigma(flipdim(v,1),zv,flipud(Z));
  temp=ztosigma(flipdim(temp,1),zr,flipud(Z));
  salt=ztosigma(flipdim(salt,1),zr,flipud(Z));  
%
% Correct the horizontal transport 
% i.e. remove the interpolated tranport and add 
%      the OGCM transport
%

%   if conserv==1
%     u=u-tridim(squeeze(sum(u.*dzu)./sum(dzu)),N);
%     v=v-tridim(squeeze(sum(v.*dzv)./sum(dzv)),N);
%     u=u+tridim(ubar,N);
%     v=v+tridim(vbar,N);
%   end

%
% Barotropic velocities
%
% %   ubar=squeeze(sum(u.*dzu)./sum(dzu));
% %   vbar=squeeze(sum(v.*dzv)./sum(dzv));
%
end   %~isempty(nc_clm)

%
%
% Vertical interpolation in case of bry files
%
%
% if ~isempty(nc_bry)
%   if obc(1)==1
%     [u_south,v_south,ubar_south,vbar_south,...
%      temp_south,salt_south]=vinterp_OGCM_bry(zr(:,1,:),zu(:,1,:),zv(:,1,:),...
%                                              dzr(:,1,:),dzu(:,1,:),dzv(:,1,:),...
%                                              u_south,v_south,ubar_south,vbar_south,...
%                                              temp_south,salt_south,...
%                                              N,Z,conserv);
%   end
%   if obc(2)==1 
%     [u_east,v_east,ubar_east,vbar_east,...
%      temp_east,salt_east]=vinterp_OGCM_bry(zr(:,:,end),zu(:,:,end),zv(:,:,end),...
%                                            dzr(:,:,end),dzu(:,:,end),dzv(:,:,end),...
%                                            u_east,v_east,ubar_east,vbar_east,...
%                                            temp_east,salt_east,...
%                                            N,Z,conserv);
%   end
%   if obc(3)==1
%     [u_north,v_north,ubar_north,vbar_north,...
%      temp_north,salt_north]=vinterp_OGCM_bry(zr(:,end,:),zu(:,end,:),zv(:,end,:),...
%                                              dzr(:,end,:),dzu(:,end,:),dzv(:,end,:),...
%                                              u_north,v_north,ubar_north,vbar_north,...
%                                              temp_north,salt_north,...
%                                              N,Z,conserv);
%   end
%   if obc(4)==1  
%     [u_west,v_west,ubar_west,vbar_west,...
%      temp_west,salt_west]=vinterp_OGCM_bry(zr(:,:,1),zu(:,:,1),zv(:,:,1),...
%                                            dzr(:,:,1),dzu(:,:,1),dzv(:,:,1),...
%                                            u_west,v_west,ubar_west,vbar_west,...
%                                            temp_west,salt_west,...
%                                            N,Z,conserv);
%   end
% end   %~isempty(nc_bry)
%--------------------------------------------------------------

%
%  fill the files
%
% Climatology file
%
if ~isempty(nc_clm)
%   nc_clm{'zeta'}(tout,:,:)=zeta;
%   nc_clm{'SSH'}(tout,:,:)=zeta;
  nc_clm{'temp'}(tout,:,:,:)=temp;
  nc_clm{'salt'}(tout,:,:,:)=salt;
%   nc_clm{'u'}(tout,:,:,:)=u;
%   nc_clm{'v'}(tout,:,:,:)=v;
%   nc_clm{'ubar'}(tout,:,:,:)=ubar;
%   nc_clm{'vbar'}(tout,:,:,:)=vbar;
end
%
% Boundary file
%
% if ~isempty(nc_bry)
%   if obc(1)==1
%     nc_bry{'zeta_south'}(tout,:)=zeta_south;
%     nc_bry{'temp_south'}(tout,:,:)=temp_south;
%     nc_bry{'salt_south'}(tout,:,:)=salt_south;
%     nc_bry{'u_south'}(tout,:,:)=u_south;
%     nc_bry{'v_south'}(tout,:,:)=v_south;
%     nc_bry{'ubar_south'}(tout,:,:)=ubar_south;
%     nc_bry{'vbar_south'}(tout,:,:)=vbar_south;
%   end  
%   if obc(2)==1
%     nc_bry{'zeta_east'}(tout,:)=zeta_east;
%     nc_bry{'temp_east'}(tout,:,:)=temp_east;
%     nc_bry{'salt_east'}(tout,:,:)=salt_east;
%     nc_bry{'u_east'}(tout,:,:)=u_east;
%     nc_bry{'v_east'}(tout,:,:)=v_east;
%     nc_bry{'ubar_east'}(tout,:,:)=ubar_east;
%     nc_bry{'vbar_east'}(tout,:,:)=vbar_east;  
%   end  
%   if obc(3)==1
%     nc_bry{'zeta_north'}(tout,:)=zeta_north;
%     nc_bry{'temp_north'}(tout,:,:)=temp_north;
%     nc_bry{'salt_north'}(tout,:,:)=salt_north;
%     nc_bry{'u_north'}(tout,:,:)=u_north;
%     nc_bry{'v_north'}(tout,:,:)=v_north;
%     nc_bry{'ubar_north'}(tout,:,:)=ubar_north;
%     nc_bry{'vbar_north'}(tout,:,:)=vbar_north;    
%   end  
%   if obc(4)==1
%     nc_bry{'zeta_west'}(tout,:)=zeta_west;
%     nc_bry{'temp_west'}(tout,:,:)=temp_west;
%     nc_bry{'salt_west'}(tout,:,:)=salt_west;
%     nc_bry{'u_west'}(tout,:,:)=u_west;
%     nc_bry{'v_west'}(tout,:,:)=v_west;
%     nc_bry{'ubar_west'}(tout,:,:)=ubar_west;
%     nc_bry{'vbar_west'}(tout,:,:)=vbar_west;  
%   end
% end
