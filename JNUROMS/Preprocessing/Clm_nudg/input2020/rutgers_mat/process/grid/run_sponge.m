clear

%  sponge (Gname, factor, Nfilter, Lplot, Lwrite)
Gname='/home/jhlee/model/ROMS_latest/project/NWP4/input/NWP4_grd_3_test.nc';
factor=5;
Nfilter=10;
Lplot=true;
Lwrite=true;
S = sponge (Gname, factor, Nfilter, Lplot, Lwrite);