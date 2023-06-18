% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color
% figure(1);clf
% f_path='F:\case\ecsy10\bulk\';
% fn=[f_path 'ecsy10_bulk_wd_96_6hour.nc'];
clear
gn=grd('ecsy10');
% 
% plot_path='D:\figure/';
% for day=210:.25:233
%     clf
% % day=1;
% roms_plot_bulk_winds(fn,'wind',day,gn)
% eval(['print -depsc2 -painters ana_wind' num2str(day) '.eps'])
% end
% close
plot_path='D:\figure/';
file='F:\case\ecsy10\1996\avg_ecsy10_blk_96_lmd_0004.nc';
for ii=1:60
    figure(1);    clf
    roms_sview(file,'salt',ii,20,gn)%,0,2,1,'b')
    roms_addvect(file,'sustr',ii,0,gn,3,10,'k')
    
%     roms_addvect_scale(ginput(1),[1],10)
    
    print ('-dpng',[plot_path '96_CDW_wd_str' num2str(ii) '.png'])
end

