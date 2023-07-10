function grd = gg(location)
  
if nargin == 0
  location = 'latte'; % default
end

switch location
  case 'd_PAC'
    grd_file = 'down_PAC_grd.nc';
    scoord = [5 0.4 50 21];
%   case 'ecsy20'
%     grd_file = 'D:\MODEL\DATA\ecsy20\ECSY_ke5_grd.nc';
%     scoord = [5 0.4 50 20];
%   case 'ecs10'
%     grd_file = 'D:\MODEL\DATA\ecs10\NWP8_grd.nc';
%     scoord = [5 0.4 50 30];
%   case 'west_sea'
%     grd_file = '../J5/westsea/westsea_grd.nc';
%     scoord = [5 0.4 50 30];
%   case 'PAC'
%     grd_file = 'D:\MODEL\DATA\PAC\PWS_grid_2-1.nc';
%     scoord = [5 0.4 50 21];
%   case 'nw025'
%     grd_file = 'D:\MODEL\DATA\nwp025_n\nw025_grd9.nc';
%     scoord = [5 0.4 10 20];
%   case 'ecsk_1'
%     grd_file = 'D:\MODEL\DATA\ecsk_1\ecsk1_grd.nc';
%     scoord = [5 0.4 50 30];
%   case 'ecsk_2'
%     grd_file = 'D:\MODEL\DATA\ecsk_2\ecsk2_grd.nc';
%     scoord = [5 0.4 50 30];
%   case 'ecsk_3'
%     grd_file = 'D:\MODEL\DATA\ecsk_2\ecsk2_grd_2.nc';
%     scoord = [5 0.4 50 30];
%   case 'ecsy12_jao'
%     grd_file = 'D:\MODEL\DATA\ecsy12\ecsy12_jao_grd.nc';
%     scoord = [5 0.4 50 20];
%   case 'spirit'
%     grd_file = 'D:\MODEL\DATA\spirit\spirit_grd_mask.nc';
%     scoord = [5 0.4 50 20];
%   case 'spirit2'
%     grd_file = 'D:\MODEL\DATA\spirit\spirit_grd_24.nc';
%     scoord = [5 0.4 50 20];
%   case 'ecsy16'
%     grd_file = 'D:\MODEL\DATA\ecsy16\ecsy16_grd.nc';
%     scoord = [5 0.4 50 20];
%   case 'ecsy10'
%     grd_file = 'D:\MODEL\DATA\ecsy10\ecsy10_grd.nc';
%     scoord = [5 0.4 50 20];
%   case 'ecsy10_30'
%     grd_file = 'D:\MODEL\DATA\ecsy10\ecsy10_grd.nc';
%     scoord = [5 0.4 50 30];
%   case 'ecsy10_flat'
%     grd_file = 'D:\MODEL\DATA\ecsy10\ecsy10_flat_40_grd.nc';
%     scoord = [5 0.4 40 20];
%   case 'ecsy12'
%     grd_file = 'D:\MODEL\DATA\ecsy12\ecsy12_grd.nc';
%     scoord = [5 0.4 5 20];
%   case 'ecss'
%     grd_file = 'D:\MODEL\DATA\ecss\ecss_grd.nc';
%     scoord = [7 0.1 30 20 2 2];
%   case 'ecsy12_1'
%     grd_file = 'D:\MODEL\DATA\ecsy12\ecsy12_grd_geb.nc';
%     scoord = [5 0.4 15 20];
%   case 'yecs22_50'
%     grd_file = 'D:\MODEL\DATA\yecs22\yecs22_grd.nc';
%     scoord = [7 0.1 50 20 2 2];
%   case 'yecs22'
%     grd_file = 'D:\MODEL\DATA\yecs22\yecs22_grd.nc';
%     scoord = [5 0.4 15 20 2 2];
%   case 'yecs22_geb'
%     grd_file = 'D:\MODEL\DATA\yecs22\yecs22_grd_geb.nc';
%     scoord = [5 0.4 15 20 2 2];
%   case 'yecs22_etp5'
%     grd_file = 'D:\MODEL\DATA\yecs22\yecs22_grd_etp5.nc';
%     scoord = [5 0.4 15 20 2 2];
%   case 'semi_NWP12'
%     grd_file = 'D:\MODEL\DATA\Semi_NWP12\Semi_NWP12_grd_R0p25_smth_50.nc';
%     scoord = [7 0.1 50 30 2 2];
%   case 'ecsy6'
%     grd_file = 'D:\MODEL\DATA\ecsy6\ecsy6_grd.nc';
%     scoord = [5 0.4 50 20];
%   case 'ecsy6_v2'
%     grd_file = 'D:\MODEL\DATA\ecsy6\ecsy6_grd.nc';
%     scoord = [7 0.1 50 20 2 2];
%   case 'NWP4'
%     grd_file = 'D:\MODEL\DATA\NWP4\NWP4_grd.nc';
%     scoord = [5 0.4 50 20];
%   case 'Navy'
%     grd_file = 'D:\MODEL\DATA\kangjung\navy4\Navy_mod4_grd.nc';
%     scoord = [7 0.1 50 10];
%   case 'Navy_mod2'
%     grd_file = 'D:\MODEL\DATA\kangjung\navy5\Navy_mod5_grd.nc';
%     scoord = [7 0.1 50 10];
%   case 'Navy_mod10'
%     grd_file = 'D:\MODEL\DATA\kangjung\navy5\Navy_mod5_grd.nc';
%     scoord = [5 0.4 50 5];
%   case 'Navy_mod7'
%     grd_file = 'D:\MODEL\DATA\kangjung\navy7\Navy_mod7_mask_grd.nc';
%     scoord = [5 0.4 50 5];
%   case 'jeju'
%     grd_file = 'D:\MODEL\DATA\jeju_kor\jeju_new_10lev\jeju90_grd.nc';
%     scoord = [7 0.1 50 10];
%   case 'pacific'
%     grd_file = 'D:\MODEL\DATA\pacific\pac015_grd_adopt_high_skk_r15_depth.nc';
%     scoord = [7 0.1 150 20 2 2];
%   case 'PAC50'
%     grd_file = 'D:\MODEL\DATA\PACIFIC_50\PAC50_grd.nc';
%     scoord = [5 0.4 50 20];
%     grd_file = 'D:\MODEL\DATA\PACIFIC_50\PAC50_grd_hc20.nc';
%     scoord = [5 0.4 20 20];
%   case 'GNWP6'
%     grd_file = 'D:\MODEL\DATA\GNWP6\GNWP6_grd.nc';
%     scoord = [5 0.4 20 20];
%   case 'NWP6'
%     grd_file = 'D:\MODEL\DATA\NWP6\NWP6_grd.nc';
%     scoord = [5 0.4 10 20];
%   case 'NWP6_PAC'
%     grd_file = 'D:\MODEL\DATA\NWP6_PAC\NWP6_grd_pactopo.nc';
%     scoord = [5 0.4 10 20];
%   case 'NWP6_PAC2'
%     grd_file = 'D:\MODEL\DATA\NWP6_PAC\NWP6_grd_pactopo_hc10.nc';
%     scoord = [5 0.4 10 20];
%   case 'NWP6_ETP'
%     grd_file = 'D:\MODEL\DATA\NWP6_PAC\NWP6_grd_etp5_30sec.nc';
%     scoord = [5 0.4 10 20];
%   case 'NWP6_ETP_KT_KURHO'
%     grd_file = 'D:\MODEL\DATA\NWP6_PAC\NWP6_grd_etp5_30sec_kt_kurho.nc';
%     scoord = [5 0.4 10 20];
  case 'jeju2'
    grd_file = '/data/Lab/SCHA/MODEL/ROMS839/project/JEJU_MODEL/INPUT/JEJU2_grd_for_yecs15.nc';
%     scoord = [7 0.1 10 16];
    scoord = [7 0.1 10 16 2 2];
%     theta_s = scoord(1);
%     theta_b = scoord(2);
%     Tcline  = scoord(3);
%     N       = scoord(4);
%     Vtransform = scoord(5);
%     Vstretching = scoord(6);
  case 'yecs15'
    grd_file = '/data/JHMoon/model/ROMS839/project/yecs15_auto/input/yecs15_grd2.nc';
    scoord = [7 0.1 30 20 2 2];
end

disp(' ')
disp([ 'Loading ROMS grd for application: ' location])
disp([ 'using grid file ' grd_file])
disp(' ')
grd = roms_get_grid(grd_file,scoord);
