
addpath('D:\OneDrive\Sources\Tools\Roms_tools_ori\Preprocessing_tools')
Grd_name='G:/MODEL_DATA/Grd/Grd_SO_05d_sponge.nc';
Data_name='D:/OneDrive/base142/Factory/MantaROMS/Test_nc/Ini_soda_05d_jhlee_198002.nc';

data=netcdf(Data_name);
ncG=netcdf(Grd_name);

LON=ncG{'lon_rho'}(:);
LAT=ncG{'lat_rho'}(:);

Coord1=find( LAT(:,1)>=-60.1 & LAT(:,1)<=-59.9) ;
zeta=data{'zeta'}(:);

z_rho=squeeze(zlevs(1, 1, 5, 0.1, 300, 50, ...
                       1, ncG{'h'}(Coord1,:),zeta(Coord1,:) , 0));
z_w=squeeze(zlevs(1, 1, 5, 0.1, 300, 50, ...
                       5, ncG{'h'}(Coord1,:), zeta(Coord1,:), 0));


z_rho = zlevs(ncG{'h'}(Coord1,:),zeta(Coord1,:),5,0.1,300,50,'r',1);
z_w = zlevs(ncG{'h'}(Coord1,:),zeta(Coord1,:),5,0.1,300,50,'w',1);




% [z]=zlevs(Vtransform, Vstretching,theta_s, theta_b, hc, N, ...
%                        igrid, h, zeta, report);
%    igrid         Staggered grid C-type (integer):
%                    igrid=1  => density points
%                    igrid=2  => streamfunction points
%                    igrid=3  => u-velocity points
%                    igrid=4  => v-velocity points
%                    igrid=5  => w-velocity points

figure(1)
hold on
for i=1:50
    plot(LON(1,:),z_rho(i,:),color='k')
    plot(LON(1,:),z_w(i,:),color='b')
end
 plot(LON(1,:),zeta(Coord1,:),color='r')














