function plot_vertical_profile=plot_roms_zprofile(fname,grids,var,xy)
%  plot_roms_zprofile.m uses roms_zprofile.m
%  and plot a number of z-profiles.
%
%  You may want to change output plot format.
%
%  oct 2004, BJ Choi
%
%=====================================================
%  fname = history or average file name: fname='his_eas_xxx.nc'
%  grids = grid structure use grd.m    : grids=grd('eas');
%  var   = 'temp' 'salt' 'u' 'v' 'w'
%  xy    = a matrix of points where you want to plot the vertical profile.
%            | lon1 lat1 |
%            | lon2 lat3 |
%            | lon3 lat4 |
%            |  ... ...  |
%            | lonN latN |
%===================================================== 
% EXAMPLE on matlab window
% >> f='his_eas_008_0043.nc';
% >> g=grd('eas');
% >> xy=[ 130 25; 135 38; 125 35];
% >> plot_roms_zprofile(f,g,'salt',xy)
% >> plot_roms_zprofile(f,g,'temp',xy)
% >> plot_roms_zprofile(f,g,'u',xy)
% >> plot_roms_zprofile(f,g,'v',xy)
% >> plot_roms_zprofile(f,g,'w',xy)
%=====================================================

 nc=netcdf(fname,'read');
 data_ori=squeeze( nc{var}(:) );
 dimen=size(data_ori);
 if (length(dimen) == 3)
   data=data_ori;
 elseif (length(dimen) == 4)
   data=squeeze(data_ori(1,:,:,:));
 else
   disp(' Error in dimension of the variable ')
   stop
 end
 close(nc)
 [npoint dum]=size(xy);   %find the number of points
 [out,zout]=roms_zprofile(data,grids,xy);

 for i=1:npoint

  figure
  depth=squeeze( zout(:,i) );
  zprof=squeeze( out(:,i)  );
  plot(zprof,depth,'o-')

  ylabel('depth (m)')
  if (var == 'temp')
   xlabel('temperature (^{o}C)')
  elseif (var == 'salt')
   xlabel('Salinity (psu)')
  elseif (var == 'u')
   xlabel('north-south velocity (m/s)')
  elseif (var == 'v')
   xlabel('east-west velocity (m/s)')
  elseif (var == 'w')
   xlabel('vertical velocity (m/s)')
  end
  title([ var ' at ' num2str(xy(i,1)) '^{o}E ' num2str(xy(i,2)) '^{o}N' ])

 end

