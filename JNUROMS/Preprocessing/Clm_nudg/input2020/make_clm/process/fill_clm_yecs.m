clear; clc

% par_dir = '/data/Lab/SCHA/MODEL/ROMS839/project/JEJU_MODEL/Make_yecs/make_clm/M07_natural_vlinear/';
% par_dir = '/data/Lab/SCHA/MODEL/ROMS839/project/JEJU_MODEL/Make_yecs/make_clm/M01_natural_vlinear/';
par_dir = '/data/Lab/SCHA/MODEL/ROMS839/project/Jelly_dr/input2020/make_clm/process/clmfiles_yecs15/';


% ncdisp(par_file)

clm_dir = './';
% clmfile = ([clm_dir,'JEJU2_clm_yecs15_Y2017M7-inter1-vlinear.nc']);
% clmfile = ([clm_dir,'JEJU2_clm_yecs15_Y2017M7-inter1-vspline.nc']);

% clmfile = ([clm_dir,'JEJU2_clm_yecs15_Y2017M1-inter1-vlinear.nc']);
clmfile = ([clm_dir,'Jelly2021_clm_yecs15_Y2021M3-inter1-vlinear.nc']);

ncdisp(clmfile)

nc_clm=netcdf(clmfile,'write');

tout = 0;
% for iday = 181 : 212 % start 6/30 end 7/31
% for iday = 0 : 31 % start 12/31 end 1/31
for iday = 10286 : 10290 % start 12/31 end 1/31    
    tout = tout + 1 

% parfile = ([par_dir,'JEJU2_clm_partition_yecs_v2_17_0',num2str(iday,'%3.3i'),'_natural_vspline.nc']);
% parfile = ([par_dir,'JEJU2_clm_partition_yecs_v2_17_0',num2str(iday,'%3.3i'),'_natural_vlinear.nc']);
parfile = ([par_dir,'Jelly2021_clm_partition_yecs15_',num2str(iday,'%5.5i'),'_natural_vlinear.nc']);

nc = netcdf(parfile);
temp = nc{'temp'}(:);
salt = nc{'salt'}(:);
time = nc{'ocean_time'}(:);
% theta_s = nc{'theta_s'}(:);
close(nc);clear nc

  nc_clm{'temp'}(tout,:,:,:)=temp;
  nc_clm{'salt'}(tout,:,:,:)=salt;
  nc_clm{'tclm_time'}(tout,:,:,:)=time./86400;
  nc_clm{'temp_time'}(tout,:,:,:)=time./86400;
  nc_clm{'sclm_time'}(tout,:,:,:)=time./86400;
  nc_clm{'salt_time'}(tout,:,:,:)=time./86400;
  
end

close(nc_clm); clear nc_clm

%% check
% 
% nc_clm=netcdf(clmfile);
% temp1 = nc_clm{'temp'}(:);
% tt = nc_clm{'tclm_time'}(:);
% aa = squeeze(temp1(5,:,:,:));
% 
% nc = netcdf([par_dir,'JEJU2_clm_partition_yecs_v2_17_0185_natural_vspline.nc']);
% bb = nc{'temp'}(:);
% 
% cc = aa - bb;

