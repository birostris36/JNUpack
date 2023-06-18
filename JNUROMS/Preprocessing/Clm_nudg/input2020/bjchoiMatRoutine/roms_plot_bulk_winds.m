function han = roms_plot_winds(frc_file,var,day,grd_file)
% roms_plot_winds(frc_file,var,day,grd_file)
%
% plot the surface fluxes (wind,heat etc) from a ROMS forcing nc file, using
% the data nearest to the requested day
%
% Input 'var' can be any valid surface flux variable in the netcdf file
% If var = 'wind' a vector plot is made of the wind stress components
%
% John Wilkin

nc = netcdf(frc_file);

if isempty(day)
  sms_time = nc{'sms_time'}(:);
  disp(['sms_time values are in the range ' mat2str(range(sms_time))])
  disp(['apparent time interval is ' num2str(diff(sms_time(1:2))) ' day'])
  return
end

sms_time = nc{'wind_time'}(:);
tmp = abs(day-sms_time(:));
tindex = min(find(tmp==min(tmp)));

grd = roms_get_grid(grd_file);
lon = grd.lon_rho;
lat = grd.lat_rho;
mask = grd.mask_rho_nan;

j = tindex;
clf
uscale=.1;
switch var
    case 'wind'
        sustr = nc{'Uwind'}(tindex,:,:);
        svstr = nc{'Vwind'}(tindex,:,:);
%         u = sustr.*grd.mask_rho;
%         v = svstr.*grd.mask_rho;
        u = sustr;
        v = svstr;
        %Peter-------------------------------------
        %     roms_quiver(lon,lat,u,v,grd.angle,2,0)
    roms_quiver(lon,lat,u,v,grd.angle,4,uscale,'linewidth',1,'color','k')
    roms_addvect_scale([128.5 36],[5 0],uscale,'m/s','k');
    plot_korea_coast(2,'k')
    otherwise
end
close(nc)
day_str = [' day ' num2str(sms_time(j))];
title([var ': ' day_str])
%Peter-------------------------------------
% set(gcf,'color','w','position',[461 23 1100 950])

% plotnzb(0,'k');axiseauc %peter------------------


