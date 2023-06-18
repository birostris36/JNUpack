function  create_forcing_bulk(frcname,grdname,title,coadst,...
                         coadsc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf forcing file
%       frcname: name of the forcing file
%       grdname: name of the grid file
%       title: title in the netcdf file  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc=netcdf(grdname);
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
result=close(nc);
Lp=L+1;
Mp=M+1;

nw = netcdf(frcname, 'clobber');
result = redef(nw);

%
%  Create dimensions
%

nw('xi_u') = L;
nw('eta_u') = Mp;
nw('xi_v') = Lp;
nw('eta_v') = M;
nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('xi_psi') = L;
nw('eta_psi') = M;
nw('time') = length(coadst);
%
%  Create variables and attributes
%
nw{'coads_time'} = ncdouble('time');
nw{'coads_time'}.long_name = ncchar('forcing observations time');
nw{'coads_time'}.long_name = 'forcing observations time';
nw{'coads_time'}.units = ncchar('days');
nw{'coads_time'}.units = 'days';
nw{'coads_time'}.cycle_length = coadsc;

nw{'Uwind'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Uwind'}.long_name = ncchar('COADS (10m) u winds [m/s] ');
nw{'Uwind'}.long_name = 'COADS (10m) u winds [m/s] ';
nw{'Uwind'}.units = ncchar('meter second-1');
nw{'Uwind'}.units = 'meter second-1';
nw{'Uwind'}.time  = ncchar('coads_time');

nw{'Vwind'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Vwind'}.long_name = ncchar('COADS (10m) v winds [m/s] ');
nw{'Vwind'}.long_name = 'COADS (10m) v winds [m/s]';
nw{'Vwind'}.units = ncchar('meter second-1');
nw{'Vwind'}.units = 'meter second-1';
nw{'Vwind'}.time  = ncchar('coads_time');

nw{'Tair'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Tair'}.long_name = ncchar('sea level air temperature');
nw{'Tair'}.long_name = 'sea level air temperature';
nw{'Tair'}.units = ncchar('Celsius');
nw{'Tair'}.units = 'Celsius';
nw{'Tair'}.time  = ncchar('coads_time');

nw{'Qair'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Qair'}.long_name = ncchar('relative humidity');
nw{'Qair'}.long_name = 'relative humidity';
nw{'Qair'}.units = ncchar('percentage');
nw{'Qair'}.units = 'percentage';
nw{'Qair'}.time  = ncchar('coads_time');

nw{'Pair'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'Pair'}.long_name = ncchar('sea level air presure');
nw{'Pair'}.long_name = 'sea level air presure';
nw{'Pair'}.units = ncchar('mbar');
nw{'Pair'}.units = 'mbar';
nw{'Pair'}.time  = ncchar('coads_time');

nw{'swrad'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'swrad'}.long_name = ncchar('solar shortwave radiation');
nw{'swrad'}.long_name = 'solar shortwave radiation';
nw{'swrad'}.units = ncchar('Watts meter-2');
nw{'swrad'}.units = 'Watts meter-2';
nw{'swrad'}.positive = ncchar('downward flux, heating');
nw{'swrad'}.positive = 'downward flux, heating';
nw{'swrad'}.negative = ncchar('upward flux, cooling');
nw{'swrad'}.negative = 'upward flux, cooling';
nw{'swrad'}.time  = ncchar('coads_time');

nw{'swflux'} = ncdouble('time', 'eta_rho', 'xi_rho');
nw{'swflux'}.long_name = ncchar('surface freshwater flux (E-P)');
nw{'swflux'}.long_name = 'surface freshwater flux (E-P)';
nw{'swflux'}.units = ncchar('centimeter day-1');
nw{'swflux'}.units = 'centimeter day-1';
nw{'swflux'}.positive = ncchar('net evaporation');
nw{'swflux'}.positive = 'net evaporation';
nw{'swflux'}.negative = ncchar('net precipitation');
nw{'swflux'}.negative = 'net precipitation';
nw{'swflux'}.time  = ncchar('coads_time');

result = endef(nw);

%
% Create global attributes
%

nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.grd_file = ncchar(grdname);
nw.grd_file = grdname;
nw.type = ncchar('ROMS forcing file');
nw.type = 'ROMS forcing file';

%
% Write time variables
%

nw{'coads_time'}(:) = coadst;

close(nw);
