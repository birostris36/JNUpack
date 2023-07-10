clear; clc


nc = netcdf('JEJU2_ini_by_yecs15_M04D01.nc');
% ncdisp('JEJU2_ini_by_yecs15_M04D01.nc')
time = nc{'ocean_time'}(:);
temp = nc{'temp'}(:);
salt = nc{'salt'}(:);
u = nc{'u'}(:);
v = nc{'v'}(:);
clear nc
mapdir = '/data/Lab/SCHA/DATA_SCHA/MAP/';
mapname = [mapdir,'navy_coastline_u.mat'];

yecsdir = '/data/JHMoon/model/ROMS839/project/yecs15_auto/output/';


nw = netcdf([yecsdir,'yecs15_avg_0091.nc']);
yt = nw{'temp'}(:);
ys = nw{'salt'}(:);
yu = nw{'u'}(:);
yv = nw{'v'}(:);
clear nw

nw = netcdf('JEJU2_ini_by_yecs15_M04D01_v1.nc');
yt = nw{'temp'}(:);
ys = nw{'salt'}(:);
yu = nw{'u'}(:);
yv = nw{'v'}(:);
clear nw

nc = netcdf('/data/Lab/SCHA/MODEL/ROMS839/project/JEJU_MODEL/INPUT/JEJU2_grd_201000.nc');
lonn = nc{'lon_rho'}(:);
latn = nc{'lat_rho'}(:);
lonnu = nc{'lon_u'}(:);
latnu = nc{'lat_u'}(:);
clear nc

nw = netcdf('/data/JHMoon/model/ROMS839/project/yecs15_auto/input/yecs15_grd2.nc');
lonl = nw{'lon_rho'}(:);
latl = nw{'lat_rho'}(:);
lonlu = nw{'lon_u'}(:);
latlu = nw{'lat_u'}(:);
clear nw



% aa = squeeze(temp(16,:,:));
% bb = squeeze(yt(20,:,:));
% aa = squeeze(temp(1,:,:));
% bb = squeeze(yt(1,:,:));
% idn = find(lonn == 125);
% idl = find(lonn == 125);
aa = squeeze(temp(:,1,:)); % vertical
bb = squeeze(yt(:,1,:));
% aa = squeeze(salt(16,:,:));
% bb = squeeze(ys(20,:,:));
% aa = squeeze(salt(1,:,:));
% bb = squeeze(ys(1,:,:));

% aa = squeeze(u(16,:,:));
% bb = squeeze(yu(20,:,:));

close all;clf

m_proj('mercator','lon',[124.5 128],'lat',[32.5 34])
m_contourf(lonn,latn,aa,8:0.5:20); % temp
caxis([8 20]);
% m_contourf(lonn,latn,aa,30:0.5:35); % salt
% caxis([30 35]);
% m_contourf(lonnu,latnu,aa,-0.06:0.005:0.06); % uv bar
% caxis([-0.06 0.06]);

colorbar;colormap(jet)
m_usercoast(mapname,'patch',[.7 .7 .7])
m_grid('tickdir','out','box','on','fontname','helvetica','fontsize',20)

figure
m_proj('mercator','lon',[124.5 128],'lat',[32.5 34])
m_contourf(lonl,latl,bb,8:0.5:20);%temp
caxis([8 20]);
% m_contourf(lonl,latl,bb,30:0.5:35); % salt
% caxis([30 35]);
% m_contourf(lonlu,latlu,bb,-0.06:0.005:0.06); % uv bar
% caxis([-0.06 0.06]);
colorbar;colormap(jet)
m_usercoast(mapname,'patch',[.7 .7 .7])
m_grid('tickdir','out','box','on','fontname','helvetica','fontsize',20)



clf
contourf(aa);
colorbar;colormap(jet);caxis([8 20]);
figure
contourf(bb);
colorbar;colormap(jet);caxis([8 20]);


%% check bry

clear; clc

nc = netcdf('JEJU2_bry_by_yecs15_M04.nc');
ncdisp('JEJU2_bry_by_yecs15_M04.nc')
te = nc{'temp_east'}(:);
tw = nc{'temp_west'}(:);
ts = nc{'temp_south'}(:);
tn = nc{'temp_north'}(:);

se = nc{'salt_east'}(:);
sw = nc{'salt_west'}(:);
ss = nc{'salt_south'}(:);
sn = nc{'salt_north'}(:);

ue = nc{'u_east'}(:);
uw = nc{'u_west'}(:);
us = nc{'u_south'}(:);
un = nc{'u_north'}(:);

contourf(squeeze(te(:,16,:)));colorbar;colormap(jet);caxis([8 20])








