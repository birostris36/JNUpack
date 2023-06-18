clear
g=grd('ecsy10');
varname='temp';
itime=6;
point1=[124 126];
point2=[33 35];
f_pth='G:\case\ecsy10\1996\mod_tide\';
ncfile=[f_pth 'avg_ecsy10_blk_96_lmd_96td_0004.nc'];
rslice_vertical_p2(ncfile,varname,itime,g,point1,point2)