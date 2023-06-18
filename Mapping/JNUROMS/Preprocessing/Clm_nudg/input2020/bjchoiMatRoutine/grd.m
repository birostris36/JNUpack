function grd = gg(location)
  
if nargin == 0
  location = 'latte'; % default
end

switch location
    case 'ecs10'
        grd_file = '\Roms_tools/ecs/ecs10_moon/ecs10_grd.nc';
        scoord = [5 0.4 50 30];
    case 'ecsy10'
        grd_file = 'D:\MODEL\DATA\ecsy10\ecsy10_grd.nc';
        scoord = [5 0.4 50 20];
    case 'ecsy12'
        grd_file = 'D:\MODEL\DATA\ecsy12\ecsy12_grd.nc';
        scoord = [5 0.4 50 20];
    case 'power_s'
        grd_file = '\Roms_tools\ecs\project\power\power_s\power_s_grd.nc';
        scoord = [5 0.4 50 10];
    case 'NWP6'
        grd_file = 'D:\MODEL\DATA\NWP6\NWP6_grd.nc';
        scoord = [5 0.4 50 20];
    case 'west_sea'
        grd_file = '/Roms_tools/ecs/westsea/westsea_grd.nc';
        scoord = [5 0.4 50 30];
    case 'spirit'
        grd_file = 'D:\MODEL\DATA\spirit\spirit_grd_mask.nc';
        scoord = [5 0.4 5 20];
    case 'jeju'
        grd_file = 'D:\MODEL\DATA\jeju_kor\jeju_new_10lev\jeju90_grd.nc';
        scoord = [7 0.1 50 10];
    case 'ecs12'
        grd_file = '\Roms_tools\ecs\ecs12\ecs12_grd.nc';
        scoord = [5 0.4 50 30];
    case 'ecs12_mod'
        grd_file = '\Roms_tools\ecs\ecs12\ecs12_mod2_grd.nc';
        scoord = [5 0.4 50 30];
    case 'roms12'
        grd_file = '\Roms_tools\ecs\roms12\roms12_grd.nc';
        scoord = [5 0.4 50 30];
    case 'navy'
        grd_file = 'D:\MODEL\DATA\kangjung\navy4\Navy_mod4_grd.nc';
        scoord = [7 0.1 10 10];

end

disp(' ')
disp([ 'Loading ROMS grd for application: ' location])
disp([ 'using grid file ' grd_file])
disp(' ')
grd = roms_get_grid(grd_file,scoord);

  
