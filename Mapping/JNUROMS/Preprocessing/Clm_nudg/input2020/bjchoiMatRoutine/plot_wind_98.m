% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color
figure(1);clf
f_path='D:\old_Roms_tools\Aforc_bulk\split_bulk_forcing\';
fn=[f_path 'ecsy10_fr_wd_weekly_1998.nc'];
gn=grd('ecsy10');

plot_path='D:\figure/';
for day=5:7:264
    clf
% day=1;
roms_plot_bulk_winds(fn,'wind',day,gn)
eval([plot_path 'print -depsc2 -painters ana_wind' num2str(day) '.eps'])
end
% close