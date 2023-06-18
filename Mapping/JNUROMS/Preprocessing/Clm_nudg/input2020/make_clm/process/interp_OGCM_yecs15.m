function interp_OGCM_yecs15(OGCM_dir,OGCM_prefix,year,month,day,Roa,interp_method,...
                     lonU,latU,lonV,latV,lonT,latT,Z,tin,...
		     nc_clm,nc_bry,lon,lat,angle,h,tout,obc,vtransform,vstretching)
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
      OGCM_prefix,'Y',num2str(year),'M',num2str(month)])
  disp(num2str(datenum([year,month,day]) - datenum([1993,01,01])))
%
%
% ROMS grid angle
% %
cosa=cos(angle);
sina=sin(angle);

% 
% OGCM_dir = Model_dir;
% % OGCM_prefix
% TimeS = datenum([2017,3,31]);
% Time = datevec(TimeS);
% year = Time(1);
% month = Time(2);
% day = Time(3);
% Roa = 0;
% interp_method = 'cubic';
% tin = 1;
% % 
% nc_clm = nc_ini;
% % nc_clm = [];
% % nc_bry=netcdf(bryname,'write');
% nc_bry = [];
% tout = TimeS +1;
% % 
% obc = [1 1 1 1];
% vtransform = 2;


%
% Open the OGCM file
%
% nc=netcdf([OGCM_dir,OGCM_prefix,'Y',num2str(year),'M',num2str(month),'.cdf']);
dnum = datenum([year,month,day])-datenum([2017,01,01]) + 1;

nc=netcdf([OGCM_dir,OGCM_prefix,num2str(dnum,'%4.4i'),'.nc']);


%
% Interpole data on the OGCM Z grid and ROMS horizontal grid
%
% Get zeta because it is needed to compute vertical levels of ROMS grid
% zeta=ext_data_OGCM_yecs(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method);
zeta=ext_data_OGCM_yecs15(nc,lonT,latT,'zeta',tin,lon,lat,1,Roa,interp_method);



%
if ~isempty(nc_clm)
%
% Read and extrapole the 2D variables
%
  u2d=ext_data_OGCM_yecs15(nc,lonU,latU,'ubar',tin,lon,lat,1,Roa,interp_method);
  v2d=ext_data_OGCM_yecs15(nc,lonV,latV,'vbar',tin,lon,lat,1,Roa,interp_method);
  ubar=rho2u_2d(u2d.*cosa+v2d.*sina);
  vbar=rho2v_2d(v2d.*cosa-u2d.*sina);
  
  
  
  
%
% Read and extrapole the 3D variables
% 
%   NZ=length(Z);
  NZ=Z;
  [M,L]=size(lon);
  dz=gradient(Z);
  temp=zeros(NZ,M,L);
  salt=zeros(NZ,M,L);
  u=zeros(NZ,M,L-1);
  v=zeros(NZ,M-1,L);
  for k=1:NZ
       k1 = NZ - k + 1;
%     if rem(k,10)==0
      disp(['  Level ',num2str(k),' of ',num2str(NZ)])
%     end
    u2d=ext_data_OGCM_yecs15(nc,lonU,latU,'u_eastward',tin,lon,lat,k,Roa,interp_method);
    v2d=ext_data_OGCM_yecs15(nc,lonV,latV,'v_northward',tin,lon,lat,k,Roa,interp_method);
%     u(k1,:,:)=rho2u_2d(u2d.*cosa+v2d.*sina); %level change
%     v(k1,:,:)=rho2v_2d(v2d.*cosa-u2d.*sina); %level change
    u(k1,:,:)=rho2u_2d(u2d); %level change
    v(k1,:,:)=rho2v_2d(v2d); %level change
    temp(k1,:,:)=ext_data_OGCM_yecs15(nc,lonT,latT,'temp',tin,lon,lat,k,Roa,interp_method); %level change
    salt(k1,:,:)=ext_data_OGCM_yecs15(nc,lonT,latT,'salt',tin,lon,lat,k,Roa,interp_method); %level change
      
  end
end


%
%Initialisation in case of bry files
%
if ~isempty(nc_bry)
%   NZ=length(Z);
  NZ=Z;
  [M,L]=size(lon); 
% Read and extrapole the 2D variables
  if obc(1)==1
    zeta_south=squeeze(zeta(1,:));
    ubar1d=ext_data_OGCM_yecs15(nc,lonU,latU,'ubar',tin,...
	   squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),1,Roa,interp_method);
    vbar1d=ext_data_OGCM_yecs15(nc,lonV,latV,'vbar',tin,...
	   squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),1,Roa,interp_method);
    ubar_south=squeeze(rho2u_2d(ubar1d(1,:).*cosa(1,:)+vbar1d(1,:).*sina(1,:)));
    vbar_south=squeeze(rho2v_2d(vbar1d.*cosa(1:2,:)-ubar1d.*sina(1:2,:)));
  end
  if obc(2)==1
    zeta_east=squeeze(zeta(:,end));
    ubar1d=ext_data_OGCM_yecs15(nc,lonU,latU,'ubar',tin,...
	   squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),1,Roa,interp_method);
    vbar1d=ext_data_OGCM_yecs15(nc,lonV,latV,'vbar',tin,...
	   squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),1,Roa,interp_method);
    ubar_east=squeeze(rho2u_2d(ubar1d.*cosa(:,end-1:end)+vbar1d.*sina(:,end-1:end)));
    vbar_east=squeeze(rho2v_2d(vbar1d(:,end).*cosa(:,end)-ubar1d(:,end).*sina(:,end)));
  end
  if obc(3)==1
    zeta_north=squeeze(zeta(end,:));
    ubar1d=ext_data_OGCM_yecs15(nc,lonU,latU,'ubar',tin,...
	   squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),1,Roa,interp_method);
    vbar1d=ext_data_OGCM_yecs15(nc,lonV,latV,'vbar',tin,...
	   squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),1,Roa,interp_method);
    ubar_north=squeeze(rho2u_2d(ubar1d(end,:).*cosa(end,:)+vbar1d(end,:).*sina(end,:)));
    vbar_north=squeeze(rho2v_2d(vbar1d.*cosa(end-1:end,:)-ubar1d.*sina(end-1:end,:)));
  end
  if obc(4)==1
    zeta_west=squeeze(zeta(:,1));
    ubar1d=ext_data_OGCM_yecs15(nc,lonU,latU,'ubar',tin,...
	   squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),1,Roa,interp_method);
    vbar1d=ext_data_OGCM_yecs15(nc,lonV,latV,'vbar',tin,...
	   squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),1,Roa,interp_method);
    ubar_west=squeeze(rho2u_2d(ubar1d.*cosa(:,1:2)+vbar1d.*sina(:,1:2)));
    vbar_west=squeeze(rho2v_2d(vbar1d(:,1).*cosa(:,1)-ubar1d(:,1).*sina(:,1)));
  end   
%
% Read and extrapole the 3D variables
%
  if obc(1)==1
    temp_south=zeros(NZ,L);
    salt_south=zeros(NZ,L);
    u_south=zeros(NZ,L-1);
    v_south=zeros(NZ,L);
  end
  if obc(2)==1
    temp_east=zeros(NZ,M);
    salt_east=zeros(NZ,M);
    u_east=zeros(NZ,M);
    v_east=zeros(NZ,M-1);
  end
  if obc(3)==1
    temp_north=zeros(NZ,L);
    salt_north=zeros(NZ,L);
    u_north=zeros(NZ,L-1);
    v_north=zeros(NZ,L);
  end
  if obc(4)==1
    temp_west=zeros(NZ,M);
    salt_west=zeros(NZ,M);
    u_west=zeros(NZ,M);
    v_west=zeros(NZ,M-1);
  end
%
  for k=1:NZ
     k1 = NZ - k + 1;
%     if rem(k,10)==0
      disp(['  Level bry ',num2str(k),' of ',num2str(NZ)])
%     end
    if obc(1)==1 % Southern boundary
      t1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k1,Roa,interp_method));
      s1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k1,Roa,interp_method));
      temp_south(k1,:)=squeeze(t1d(1,:));
      salt_south(k1,:)=squeeze(s1d(1,:));
      u1d=ext_data_OGCM_yecs15(nc,lonU,latU,'u_eastward',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k1,Roa,interp_method);
      v1d=ext_data_OGCM_yecs15(nc,lonV,latV,'v_northward',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k1,Roa,interp_method);
%       u_south(k1,:)=squeeze(rho2u_2d(u1d(1,:).*cosa(1,:)+...
%                                     v1d(1,:).*sina(1,:)));
%       v_south(k1,:)=squeeze(rho2v_2d(v1d.*cosa(1:2,:)-...
%                                     u1d.*sina(1:2,:)));         
      u_south(k1,:)=squeeze(rho2u_2d(u1d(1,:)));
      v_south(k1,:)=squeeze(rho2v_2d(v1d)); 
                                
    end
    if obc(2)==1  % Eastern boundary 
      t1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k1,Roa,interp_method));
      s1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k1,Roa,interp_method));
      temp_east(k1,:)=squeeze(t1d(:,end)');
      salt_east(k1,:)=squeeze(s1d(:,end)');
      u1d=ext_data_OGCM_yecs15(nc,lonU,latU,'u_eastward',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k1,Roa,interp_method);
      v1d=ext_data_OGCM_yecs15(nc,lonV,latV,'v_northward',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k1,Roa,interp_method);
%       cff=squeeze(rho2u_2d(u1d.*cosa(:,end-1:end)+...
%                            v1d.*sina(:,end-1:end)));
%       u_east(k1,:)=cff';			   
%       cff=squeeze(rho2v_2d(v1d(:,end).*cosa(:,end)-...
%     	                   u1d(:,end).*sina(:,end)));
%       v_east(k1,:)=cff';
      cff=squeeze(rho2u_2d(u1d));
      u_east(k1,:)=cff';			   
      cff=squeeze(rho2v_2d(v1d(:,end)));
      v_east(k1,:)=cff';
    end
    if obc(3)==1  % Northern boundary
      t1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k1,Roa,interp_method));
      s1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k1,Roa,interp_method));
      temp_north(k1,:)=squeeze(t1d(end,:));
      salt_north(k1,:)=squeeze(s1d(end,:));
      u1d=ext_data_OGCM_yecs15(nc,lonU,latU,'u_eastward',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k1,Roa,interp_method);
      v1d=ext_data_OGCM_yecs15(nc,lonV,latV,'v_northward',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k1,Roa,interp_method);
%       u_north(k1,:)=squeeze(rho2u_2d(u1d(end,:).*cosa(end,:)+...
%                                     v1d(end,:).*sina(end,:)));
%       v_north(k1,:)=squeeze(rho2v_2d(v1d.*cosa(end-1:end,:)-...
%                                     u1d.*sina(end-1:end,:)));
      u_north(k1,:)=squeeze(rho2u_2d(u1d(end,:)));
      v_north(k1,:)=squeeze(rho2v_2d(v1d));
    end
    if obc(4)==1  % Western boundary
      t1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k1,Roa,interp_method));
      s1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k1,Roa,interp_method));
      temp_west(k1,:)=squeeze(t1d(:,1)');
      salt_west(k1,:)=squeeze(s1d(:,1)');
      u1d=ext_data_OGCM_yecs15(nc,lonU,latU,'u_eastward',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k1,Roa,interp_method);
      v1d=ext_data_OGCM_yecs15(nc,lonV,latV,'v_northward',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k1,Roa,interp_method);
%       cff=squeeze(rho2u_2d(u1d.*cosa(:,1:2)+...
%                            v1d.*sina(:,1:2)));
%       u_west(k1,:)=cff';
%       cff=squeeze(rho2v_2d(v1d(:,1).*cosa(:,1)-...
%                            u1d(:,1).*sina(:,1)));
%       v_west(k1,:)=cff';
      cff=squeeze(rho2u_2d(u1d));
      u_west(k1,:)=cff';
      cff=squeeze(rho2v_2d(v1d(:,1)));
      v_west(k1,:)=cff';
    end
  end
end

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
if ~isempty(nc_bry)
  theta_s=nc_bry{'theta_s'}(:);
  theta_b=nc_bry{'theta_b'}(:);
  hc=nc_bry{'hc'}(:);
  N=length(nc_bry('s_rho'));
end
%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
% Z=[100;Z;-100000];
Z = 20;
%
% ROMS vertical grid
%

                   
% zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
% % zr=zlevs(vtransform,vstretching,theta_s,theta_b,hc,N, ,h,zeta,0);
% 
% zu=rho2u_3d(zr);
% zv=rho2v_3d(zr);
% zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
% % zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
% dzr=zw(2:end,:,:)-zw(1:end-1,:,:);
% dzu=rho2u_3d(dzr);
% dzv=rho2v_3d(dzr);

zr=zlevs(vtransform,vstretching,theta_s,theta_b,hc,N,1,h,zeta,1);
zu=rho2u_3d(zr);
zv=rho2v_3d(zr);
zw=zlevs(vtransform,vstretching,theta_s,theta_b,hc,N,5,h,zeta,1);
dzr=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=rho2u_3d(dzr);
dzv=rho2v_3d(dzr);


% data=ext_data_OGCM_yecs(nc,X,Y,vname,tndx,lon,lat,k,Roa,interp_method);
% t1d=squeeze(ext_data_OGCM_yecs15(nc,lonT,latT,'temp',tin,...
%                   squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
%                   k1,Roa,interp_method));
              
disp('  YECS15 data')

ytheta_s = nc{'theta_s'}(:);
ytheta_b = nc{'theta_b'}(:);
yhc = nc{'hc'}(:);
yN = length(nc{'s_rho'}(:));

yecs_h=ext_data_OGCM_yecs15(nc,lonT,latT,'h',tin,lon,lat,1,Roa,interp_method);

% zeta=ext_data_OGCM_yecs15(nc,lonT,latT,'zeta',tin,lon,lat,1,Roa,interp_method);

yzr=zlevs(vtransform,vstretching,ytheta_s,ytheta_b,yhc,yN,1,yecs_h,zeta,1);
yzu=rho2u_3d(yzr);
yzv=rho2v_3d(yzr);
yzw=zlevs(vtransform,vstretching,ytheta_s,ytheta_b,yhc,yN,5,yecs_h,zeta,1);
ydzr=yzw(2:end,:,:)-yzw(1:end-1,:,:);
ydzu=rho2u_3d(ydzr);
ydzv=rho2v_3d(ydzr);

%
% Close the OGCM file
%
close(nc)


%
%
%
% Vertical interpolation in case of clim file
%
%
if ~isempty(nc_clm)
%
% Add a level on top and bottom with no-gradient
%
% u2 = u;
% v2 = v;
% temp2 = temp;
% salt2 = salt;


  u=cat(1,u(1,:,:),u);
  u=cat(1,u,u(end,:,:));
  v=cat(1,v(1,:,:),v);
  v=cat(1,v,v(end,:,:));
  temp=cat(1,temp(1,:,:),temp);
  temp=cat(1,temp,temp(end,:,:));
  salt=cat(1,salt,salt(end,:,:));
  salt=cat(1,salt(1,:,:),salt);
% 
% Perform the vertical interpolations 
%


yzr=cat(1,yzr(1,:,:),yzr);
yzr=cat(1,yzr,yzr(end,:,:));
yzr(1,:,:) =  yzr(1,:,:) - 50;
yzr(end,:,:) =  yzr(end,:,:) + 50;

yzu=cat(1,yzu(1,:,:),yzu);
yzu=cat(1,yzu,yzu(end,:,:));
yzu(1,:,:) =  yzu(1,:,:) - 50;
yzu(end,:,:) =  yzu(end,:,:) + 50;

yzv=cat(1,yzv(1,:,:),yzv);
yzv=cat(1,yzv,yzv(end,:,:));
yzv(1,:,:) =  yzv(1,:,:) - 50;
yzv(end,:,:) =  yzv(end,:,:) + 50;
%   u=ztosigma(flipdim(u,1),zu,flipud(Z));
%   v=ztosigma(flipdim(v,1),zv,flipud(Z));
%   temp=ztosigma(flipdim(temp,1),zr,flipud(Z));
%   salt=ztosigma(flipdim(salt,1),zr,flipud(Z));

for cy1 = 1 : min(size(lon))
    for cy2 = 1 : max(size(lon))
temp1(1:N,cy1,cy2) = interp1(yzr(:,cy1,cy2),temp(:,cy1,cy2),zr(1:N,cy1,cy2),'pchip');
salt1(1:N,cy1,cy2) = interp1(yzr(:,cy1,cy2),salt(:,cy1,cy2),zr(1:N,cy1,cy2),'pchip');
    end
end

for cy1 = 1 : min(size(lon))
    for cy2 = 1 : max(size(lon))-1
u1(1:N,cy1,cy2) = interp1(yzu(:,cy1,cy2),u(:,cy1,cy2),zu(1:N,cy1,cy2),'pchip');
    end
end

for cy1 = 1 : min(size(lon))-1
    for cy2 = 1 : max(size(lon))
v1(1:N,cy1,cy2) = interp1(yzv(:,cy1,cy2),v(:,cy1,cy2),zv(1:N,cy1,cy2),'pchip');
    end
end
  temp=flip(temp1,1);
  salt=flip(salt1,1);
  u=flip(u1,1);
  v=flip(v1,1);
  
  
  
%
% Correct the horizontal transport 
% i.e. remove the interpolated tranport and add 
%      the OGCM transport
%
  if conserv==1
    u=u-tridim(squeeze(sum(u.*dzu)./sum(dzu)),N);
    v=v-tridim(squeeze(sum(v.*dzv)./sum(dzv)),N);
    u=u+tridim(ubar,N);
    v=v+tridim(vbar,N);
  end
%
% Barotropic velocities
%
  ubar=squeeze(sum(u.*dzu)./sum(dzu));
  vbar=squeeze(sum(v.*dzv)./sum(dzv));
%
end   %~isempty(nc_clm)
%
%
% Vertical interpolation in case of bry files
%
%
% % if ~isempty(nc_bry)
% %   if obc(1)==1
% %     [u_south,v_south,ubar_south,vbar_south,...
% %      temp_south,salt_south]=vinterp_OGCM_bry(zr(:,1,:),zu(:,1,:),zv(:,1,:),...
% %                                              dzr(:,1,:),dzu(:,1,:),dzv(:,1,:),...
% %                                              u_south,v_south,ubar_south,vbar_south,...
% %                                              temp_south,salt_south,...
% %                                              N,Z,conserv);
% %   end
% %   if obc(2)==1 
% %     [u_east,v_east,ubar_east,vbar_east,...
% %      temp_east,salt_east]=vinterp_OGCM_bry(zr(:,:,end),zu(:,:,end),zv(:,:,end),...
% %                                            dzr(:,:,end),dzu(:,:,end),dzv(:,:,end),...
% %                                            u_east,v_east,ubar_east,vbar_east,...
% %                                            temp_east,salt_east,...
% %                                            N,Z,conserv);
% %   end
% %   if obc(3)==1
% %     [u_north,v_north,ubar_north,vbar_north,...
% %      temp_north,salt_north]=vinterp_OGCM_bry(zr(:,end,:),zu(:,end,:),zv(:,end,:),...
% %                                              dzr(:,end,:),dzu(:,end,:),dzv(:,end,:),...
% %                                              u_north,v_north,ubar_north,vbar_north,...
% %                                              temp_north,salt_north,...
% %                                              N,Z,conserv);
% %   end
% %   if obc(4)==1  
% %     [u_west,v_west,ubar_west,vbar_west,...
% %      temp_west,salt_west]=vinterp_OGCM_bry(zr(:,:,1),zu(:,:,1),zv(:,:,1),...
% %                                            dzr(:,:,1),dzu(:,:,1),dzv(:,:,1),...
% %                                            u_west,v_west,ubar_west,vbar_west,...
% %                                            temp_west,salt_west,...
% %                                            N,Z,conserv);
% %   end
% % end   %~isempty(nc_bry)
%--------------------------------------------------------------

%
%
% Vertical interpolation in case of bry files by scha
%
%
if ~isempty(nc_bry)
  if obc(1)==1
    [u_south,v_south,ubar_south,vbar_south,...
     temp_south,salt_south]=vinterp_OGCM_bry_yecs15(zr(:,1,:),zu(:,1,:),zv(:,1,:),...
                                             dzr(:,1,:),dzu(:,1,:),dzv(:,1,:),...
                                             yzr(:,1,:),yzu(:,1,:),yzv(:,1,:),...
                                             ydzr(:,1,:),ydzu(:,1,:),ydzv(:,1,:),...
                                             u_south,v_south,ubar_south,vbar_south,...
                                             temp_south,salt_south,...
                                             N,Z,conserv);
  end
  if obc(2)==1 
    [u_east,v_east,ubar_east,vbar_east,...
     temp_east,salt_east]=vinterp_OGCM_bry_yecs15(zr(:,:,end),zu(:,:,end),zv(:,:,end),...
                                           dzr(:,:,end),dzu(:,:,end),dzv(:,:,end),...
                                           yzr(:,:,end),yzu(:,:,end),yzv(:,:,end),...
                                           ydzr(:,:,end),ydzu(:,:,end),ydzv(:,:,end),...
                                           u_east,v_east,ubar_east,vbar_east,...
                                           temp_east,salt_east,...
                                           N,Z,conserv);
  end
  if obc(3)==1
    [u_north,v_north,ubar_north,vbar_north,...
     temp_north,salt_north]=vinterp_OGCM_bry_yecs15(zr(:,end,:),zu(:,end,:),zv(:,end,:),...
                                             dzr(:,end,:),dzu(:,end,:),dzv(:,end,:),...
                                             yzr(:,end,:),yzu(:,end,:),yzv(:,end,:),...
                                             ydzr(:,end,:),ydzu(:,end,:),ydzv(:,end,:),...
                                             u_north,v_north,ubar_north,vbar_north,...
                                             temp_north,salt_north,...
                                             N,Z,conserv);
  end
  if obc(4)==1  
    [u_west,v_west,ubar_west,vbar_west,...
     temp_west,salt_west]=vinterp_OGCM_bry_yecs15(zr(:,:,1),zu(:,:,1),zv(:,:,1),...
                                           dzr(:,:,1),dzu(:,:,1),dzv(:,:,1),...
                                           yzr(:,:,1),yzu(:,:,1),yzv(:,:,1),...
                                           ydzr(:,:,1),ydzu(:,:,1),ydzv(:,:,1),...
                                           u_west,v_west,ubar_west,vbar_west,...
                                           temp_west,salt_west,...
                                           N,Z,conserv);
  end
end   %~isempty(nc_bry)

%     [u_south,v_south,ubar_south,vbar_south,...
%      temp_south,salt_south]=vinterp_OGCM_bry(zr(:,1,:),zu(:,1,:),zv(:,1,:),...
%                                              dzr(:,1,:),dzu(:,1,:),dzv(:,1,:),...
%                                              yzr(:,1,:),yzu(:,1,:),yzv(:,1,:),...
%                                              ydzr(:,1,:),ydzu(:,1,:),ydzv(:,1,:),...
%                                              u_south,v_south,ubar_south,vbar_south,...
%                                              temp_south,salt_south,...
%                                              N,Z,conserv,lon1d,lat1d);




%
%  fill the files
%
% Climatology file
%
if ~isempty(nc_clm)
  nc_clm{'zeta'}(tout,:,:)=zeta;
  nc_clm{'SSH'}(tout,:,:)=zeta;
  nc_clm{'temp'}(tout,:,:,:)=temp;
  nc_clm{'salt'}(tout,:,:,:)=salt;
  nc_clm{'u'}(tout,:,:,:)=u;
  nc_clm{'v'}(tout,:,:,:)=v;
  nc_clm{'ubar'}(tout,:,:,:)=ubar;
  nc_clm{'vbar'}(tout,:,:,:)=vbar;
end
%
% Boundary file
%
if ~isempty(nc_bry)
  if obc(1)==1
    nc_bry{'zeta_south'}(tout,:)=zeta_south;
    nc_bry{'temp_south'}(tout,:,:)=temp_south;
    nc_bry{'salt_south'}(tout,:,:)=salt_south;
    nc_bry{'u_south'}(tout,:,:)=u_south;
    nc_bry{'v_south'}(tout,:,:)=v_south;
    nc_bry{'ubar_south'}(tout,:,:)=ubar_south;
    nc_bry{'vbar_south'}(tout,:,:)=vbar_south;
  end  
  if obc(2)==1
    nc_bry{'zeta_east'}(tout,:)=zeta_east;
    nc_bry{'temp_east'}(tout,:,:)=temp_east;
    nc_bry{'salt_east'}(tout,:,:)=salt_east;
    nc_bry{'u_east'}(tout,:,:)=u_east;
    nc_bry{'v_east'}(tout,:,:)=v_east;
    nc_bry{'ubar_east'}(tout,:,:)=ubar_east;
    nc_bry{'vbar_east'}(tout,:,:)=vbar_east;  
  end  
  if obc(3)==1
    nc_bry{'zeta_north'}(tout,:)=zeta_north;
    nc_bry{'temp_north'}(tout,:,:)=temp_north;
    nc_bry{'salt_north'}(tout,:,:)=salt_north;
    nc_bry{'u_north'}(tout,:,:)=u_north;
    nc_bry{'v_north'}(tout,:,:)=v_north;
    nc_bry{'ubar_north'}(tout,:,:)=ubar_north;
    nc_bry{'vbar_north'}(tout,:,:)=vbar_north;    
  end  
  if obc(4)==1
    nc_bry{'zeta_west'}(tout,:)=zeta_west;
    nc_bry{'temp_west'}(tout,:,:)=temp_west;
    nc_bry{'salt_west'}(tout,:,:)=salt_west;
    nc_bry{'u_west'}(tout,:,:)=u_west;
    nc_bry{'v_west'}(tout,:,:)=v_west;
    nc_bry{'ubar_west'}(tout,:,:)=ubar_west;
    nc_bry{'vbar_west'}(tout,:,:)=vbar_west;  
  end
end
