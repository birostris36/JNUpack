clear all


gn=grd('ecsy12');


endpt_lon(1) = 117;
endpt_lat(1) = 33;
endpt_lon(2) = 126.5;
endpt_lat(2) = 41;


lon_index=find(gn.lon_rho >endpt_lon(1) & gn.lon_rho <endpt_lon(2));
r_lon1=gn.lon_rho(lon_index);
r_lat1=gn.lat_rho(lon_index);
r_mask1=gn.mask_rho(lon_index);

lat_index=find(r_lat1 >endpt_lat(1) & r_lat1 <endpt_lat(2));
r_lon2=r_lon1(lat_index);
r_lat2=r_lat1(lat_index);
r_mask2=r_mask1(lat_index);
clear r_lon1 r_lat1 r_mask1
mask_index=find(r_mask2==1);
r_lon3=r_lon2(mask_index);
r_lat3=r_lat2(mask_index);
r_mask3=r_mask2(mask_index);
% clear r_lon2 r_lat2

clf
plot(r_lon2, r_lat2,'ko')
hold on
plot(r_lon3, r_lat3,'ro')

f_path=['/home/jhlee/model/ROMS35/project/ECSY12_2010_ghrsst_kkl/output/'];

for day=1:365
    filename=[f_path 'avg_ecsy12_td_2010_kkl_v2_a4_wt1_' num2str(day,'%4.4i') '.nc'];
    nc=netcdf(filename,'read');
    disp([' Day : ', num2str(day)])
    salt=squeeze(nc{'salt'}(1,20,:,:)); % salt(time, s_rho, eta_rho, xi_rho)
    ssflux=squeeze(nc{'ssflux'}(1,:,:)); % salt(time, s_rho, eta_rho, xi_rho)
    evaporation=squeeze(nc{'evaporation'}(1,:,:)); % salt(time, s_rho, eta_rho, xi_rho)
    rain=squeeze(nc{'rain'}(1,:,:)); % salt(time, s_rho, eta_rho, xi_rho)
    close(nc)
    % for salinity
    r_salt1=salt(lon_index);
    r_salt2=r_salt1(lat_index);
    r_salt3=r_salt2(mask_index);
    % for ssflux
    r_ssflux1=ssflux(lon_index);
    r_ssflux2=r_ssflux1(lat_index);
    r_ssflux3=r_ssflux2(mask_index);
    % for evaporation
    r_evap1=evaporation(lon_index);
    r_evap2=r_evap1(lat_index);
    r_evap3=r_evap2(mask_index);
    % for ssflux
    r_rain1=rain(lon_index);
    r_rain2=r_rain1(lat_index);
    r_rain3=r_rain2(mask_index);
    clear r_salt1 r_salt2 ssflux1 ssflux2 r_evap1 r_evap2 r_rain1 r_rain2
    YS_salt(day)=mean(r_salt3);
    YS_ssflux(day)=mean(r_ssflux3);
    YS_evap(day)=mean(r_evap3);
    YS_rain(day)=mean(r_rain3);
end

figure(2);clf
subplot(411)
plot(YS_salt)
ylabel('YS salt')
subplot(412)
plot(YS_ssflux)
ylabel('YS ssflux')
subplot(413)
plot(YS_evap)
ylabel('YS evap')
subplot(414)
plot(YS_rain)
ylabel('YS rain')

figure(3);clf
title('2010')
subplot(311)
plot(152:243,YS_salt(152:243))
ylabel('YS salt')
subplot(312)
plot(152:243,YS_ssflux(152:243))
ylabel('YS ssflux')
subplot(313)
plot(152:243,YS_evap(152:243)-YS_rain(152:243))
ylabel('YS evap-precip')

save YS_ssflux_2010.mat YS_salt YS_ssflux YS_evap YS_rain