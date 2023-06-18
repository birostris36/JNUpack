% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color
figure(1);clf
f_path='D:\MODEL\DATA\ecsy10\';
fn=[f_path 'ecsy10_frc_wd_1996.nc'];
gn=grd('ecsy10');

plot_path='D:\figure/';
for day=1:270
    clf
% day=1;
roms_plot_winds(fn,'wind',day,gn)
eval(['print -depsc2 -painters wind96_day_' num2str(day) '.eps'])
end
% close