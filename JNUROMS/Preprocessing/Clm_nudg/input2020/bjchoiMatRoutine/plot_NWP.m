% f='D:\NWP6_his_0132.nc';
% clear close all
g=grd('NWP6');
% xy=[ 130 25; 135 38];
% plot_roms_zprofile(f,g,'salt',xy);
roms_plot_bathy(g,'cmap',[100 250 500 1000:1000:4000],'h');