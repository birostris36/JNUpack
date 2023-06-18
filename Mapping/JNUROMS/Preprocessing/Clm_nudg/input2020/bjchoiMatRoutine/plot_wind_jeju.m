% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color
figure(1);clf
f_path='D:\MODEL\DATA\jeju_kor\jeju_new_10lev\';
fn=[f_path 'jeju90_frc_td07.nc'];
gn=grd('jeju');

plot_path='D:\';
for day=15:30:345
% for day=15:15
    clf
% day=1;
roms_plot_winds(fn,'wind',day,gn)
eval(['print -depsc2 -painters wind_jeju_day_' num2str(day) '.eps'])
end
% close