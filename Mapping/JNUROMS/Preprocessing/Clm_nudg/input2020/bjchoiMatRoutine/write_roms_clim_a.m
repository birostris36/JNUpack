function write_roms_clim(init_file,roms)
%
% append data to an existing zeros clim file
%
%

nc = netcdf(init_file,'write');
if isempty(nc)
  error([ 'Problem opening ' init_file ]);
end

% append to an existing empty initial file
time_variable_name = nc{'temp'}.time(:);
nc_tindex = length(nc{time_variable_name}(:))+1;

nc{'ocean_time'}(nc_tindex) = roms.time;

for varlist = { 'temp','salt','u','v'}
  varname = char(varlist);
  data = getfield(roms,varname);
  nc{varname}(nc_tindex,:,:,:) = data;
end

for varlist = { 'zeta','ubar','vbar'}
  varname = char(varlist);
  data = getfield(roms,varname);
  nc{varname}(nc_tindex,:,:) = data;
end

close(nc)

