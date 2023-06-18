% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color
% f='his_Navy_mod3_tpxo7_1-30_0014.nc';
% 18
% gn=grd('navy');
% clf
% time=20 ;
% roms_zview(f,'zeta',time,-1,gn,0,3,.02,'k')
% fn='F:\case\jeju\new_jeju90_lev10\jeju90_frc_Y2008.nc';
% gn=grd('jeju');
% d_path='./figure/';
for day=171:273
    clf
% day=1;
roms_plot_winds(fn,'wind',day,gn)
eval(['print -depsc2 -painters ana_wind' num2str(day) '.eps'])
end
close