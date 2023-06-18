% file,var,time,depth,grd, contr,vec_d,uscale,varargin
% eg). roms_zview(f,'temp',1,-10,3,1,'k');
%                          every 3 vector, uscale, color

clear
gn=grd('ecsy10');

plot_path='G:\case\ecsy10\1996\figure_96\';
f_pth='G:\case\ecsy10\1996\mod_tide\';
file=[f_pth 'avg_ecsy10_blk_96_lmd_96td_0004.nc'];

% warning off
var1='salt';
ii=1;
clf
point1=[124 126];
point2=[32 33];


rslice_vertical_p2(file,var1,1,gn,point1,point2)


% roms_iview(file,var1,ii,132,gn);

for ii=1:1
    figure(1);    clf
            roms_iview(file,var1,ii,32,gn);
%             roms_addvect_scale([120.3 29],[1 0],uscale,'m/s','k');
%             roms_addvect_scale([120.6 28.5],[.5 0],uscale,'m/s','k');
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
