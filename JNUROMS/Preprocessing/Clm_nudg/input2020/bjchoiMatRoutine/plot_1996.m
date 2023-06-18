% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color

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
% plot_path='F:\case\ecsy10\figure_96\pcolor_20_35_wstr\';

plot_path='F:\case\ecsy10\test\';
f_pth='F:\case\ecsy10\1996\mod_td\';
file=[f_pth 'avg_ecsy10_blk_96_lmd_96td_0004.nc'];

warning off
var1='salt';
% var2='sustr';
var2='u';
vec_d=2;  % plot vector with every X space

for ii=1:1
    figure(1);    clf
    switch var2
        case 'u'
            uscale=1;
            roms_sview(file,var1,ii,20,gn,0,vec_d,uscale,'k')
            roms_addvect_scale([120.3 29],[1 0],uscale,'m/s','k');
            roms_addvect_scale([120.6 28.5],[.5 0],uscale,'m/s','k');
        case 'sustr'
            uscale=10;
            roms_sview(file,var1,ii,20,gn)
            roms_addvect(file,'sustr',ii,20,gn,vec_d,uscale,'k')
            roms_addvect_scale([121.2 29.5],[.05 0],uscale,'N/m^2','k');
        otherwise
            title('Check out variable name')
            return
    end
    hold on
    axis([118.8 131 23.5 37])
    hold on
    text(120,32.3,'C R','fontsize',12,'fontname','times',...
        'fontangle','italic','fontweight','bold')
    text(127,36,'Korea','fontsize',13,'fontname','times',...
        'fontangle','italic','fontweight','bold')
    text(119,33.5,'China','fontsize',13,'fontname','times',...
        'fontangle','italic','fontweight','bold')
    hold on
    colorscale([0 64],[20 35],2,'horiz','position',...
        [0.57 0.15 0.15 0.02],'fontsize',8,'fontweight','bold');
    hold off
% % % %     ii=ii+60;
    if  ii<10
        print ('-dpng',[plot_path '96_CDW_str0' num2str(ii) '.png'])
    else
        print ('-dpng',[plot_path '96_CDW_str' num2str(ii) '.png'])
    end
end
