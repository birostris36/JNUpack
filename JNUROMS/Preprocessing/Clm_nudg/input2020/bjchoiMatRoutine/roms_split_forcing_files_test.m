Iname='D:\MODEL\DATA\NWP4\NWP4_fr_Y2004_rv.nc';
% origin roms_split_forcing_files.m from bj-choi's tool
TIDES=0;
RIVERS=1;
WINDS=1;
HEAT=1;

Rname='NWP4_fr_2004_river.nc';
Wname='NWP4_fr_Y2004_wind.nc';
Hname='NWP4_coads_heat.nc';

%  Extract and write tides.

if (TIDES),
  f=nc_read(Iname,'tide_period'); s=nc_write(Tname,'tide_period',f);
  f=nc_read(Iname,'tide_Ephase'); s=nc_write(Tname,'tide_Ephase',f);
  f=nc_read(Iname,'tide_Eamp'  ); s=nc_write(Tname,'tide_Eamp'  ,f);
  f=nc_read(Iname,'tide_Cphase'); s=nc_write(Tname,'tide_Cphase',f);
  f=nc_read(Iname,'tide_Cangle'); s=nc_write(Tname,'tide_Cangle',f);
  f=nc_read(Iname,'tide_Cmin'  ); s=nc_write(Tname,'tide_Cmin'  ,f);
  f=nc_read(Iname,'tide_Cmax'  ); s=nc_write(Tname,'tide_Cmax'  ,f);
end,

if (RIVERS),
  time=nc_read(Iname,'river_time');
  Ntime=length(time);
  f=nc_read(Iname,'river'          ); s=nc_write(Rname,'river',f);
  f=nc_read(Iname,'river_Xposition'); s=nc_write(Rname,'river_Xposition',f);
  f=nc_read(Iname,'river_Eposition'); s=nc_write(Rname,'river_Eposition',f);
  f=nc_read(Iname,'river_direction'); s=nc_write(Rname,'river_direction',f);
  f=nc_read(Iname,'river_flag'     ); s=nc_write(Rname,'river_flag'     ,f);
  f=nc_read(Iname,'river_Vshape'   ); s=nc_write(Rname,'river_Vshape'   ,f);
  for n=1:Ntime,
    f=nc_read(Iname,'river_time',     n); s=nc_write(Rname,'river_time'     ,f,n);
    f=nc_read(Iname,'river_transport',n); s=nc_write(Rname,'river_transport',f,n);
    f=nc_read(Iname,'river_temp',     n); s=nc_write(Rname,'river_temp'     ,f,n);
    f=nc_read(Iname,'river_salt',     n); s=nc_write(Rname,'river_salt'     ,f,n);
  end,  
end,

if (WINDS),
  time=nc_read(Iname,'sms_time');
  Ntime=length(time);
  for n=1:Ntime,  
    f=nc_read(Iname,'sms_time',n); s=nc_write(Wname,'sms_time',f,n);
    f=nc_read(Iname,'sustr'   ,n); s=nc_write(Wname,'sustr'   ,f,n);
    f=nc_read(Iname,'svstr'   ,n); s=nc_write(Wname,'svstr'   ,f,n);
  end,
end,


if (HEAT),
  f=nc_read(Iname,'srf_time'); s=nc_write(Hname,'srf_time',f);
  f=nc_read(Iname,'sst_time'); s=nc_write(Hname,'sst_time',f);
  f=nc_read(Iname,'shf_time'); s=nc_write(Hname,'shf_time',f);
  f=nc_read(Iname,'swf_time'); s=nc_write(Hname,'swf_time',f);
  f=nc_read(Iname,'sss_time'); s=nc_write(Hname,'sss_time',f);
  f=nc_read(Iname,'swrad'   ); s=nc_write(Hname,'swrad'   ,f);
  f=nc_read(Iname,'SST'     ); s=nc_write(Hname,'SST'     ,f);
  f=nc_read(Iname,'dQdSST'  ); s=nc_write(Hname,'dQdSST'  ,f);
  f=nc_read(Iname,'shflux'  ); s=nc_write(Hname,'shflux'  ,f);
  f=nc_read(Iname,'swflux'  ); s=nc_write(Hname,'swflux'  ,f);
  f=nc_read(Iname,'SSS'     ); s=nc_write(Hname,'SSS'     ,f);
end,

