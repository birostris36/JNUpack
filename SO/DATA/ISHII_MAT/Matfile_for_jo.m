clear;
run('C:/Users/shjo/OneDrive/Sources/Tools/nctoolbox-master/setup_nctoolbox.m')

% addpath('D:\matlab\ÇÔ¼ö\ocean\')

Syear = 1955;
Eyear = 2019;
Smon = 1;
Emon = 12;
Mrho = 28.106331; % sigma(0,0,35) = 28.106331 kg/m^3

Filedir_salt = 'D:/JMA/salt/';
Filedir_temp = 'D:/JMA/temp/';

grb = ncgeodataset([Filedir_temp,'temp.1955.grb2']);
lon= grb{'lon'}(:);
lat= grb{'lat'}(:);
time= grb{'time'}(:);
depth= grb{'depth_below_sea'}(:);
clear grb

% sal1 = grb.data('VAR10-4-192_FROM_34-0-0_depth_below_sea');
% sal2 = grb.data('VAR10-4-201_FROM_34-0-0_depth_below_sea');
% sal3 = grb.data('VAR10-4-202_FROM_34-0-0_depth_below_sea');

[ lonr , latr ] = meshgrid(lon,lat);
[ M , L ] = size(lonr);

cp = 4000;

temp = zeros((Eyear-Syear+1).*12,M,L);
salt = zeros((Eyear-Syear+1).*12,M,L);
nm = 0;
for Year = Syear : Eyear
    disp(['========= Processing ',num2str(Year),' Year ========='])

    Filename1 = ([Filedir_temp,'temp.',num2str(Year),'.grb2']);
    grb1 = ncgeodataset(Filename1);
    temp = grb1.data('VAR10-4-192_FROM_34-0-0_depth_below_sea');
    
    Filename2 = ([Filedir_salt,'sal.',num2str(Year),'.grb2']);
    grb2 = ncgeodataset(Filename2);
    salt = grb2.data('VAR10-4-193_FROM_34-0-0_depth_below_sea');
    save(['D:/JMA/myMAT/','Ishii_TS_',num2str(Year),'.mat'],'temp','salt','lon','lat','depth')
    clear grb1 grb2
    
end


