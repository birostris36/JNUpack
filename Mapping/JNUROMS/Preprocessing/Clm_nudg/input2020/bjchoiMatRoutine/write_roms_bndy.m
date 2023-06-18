function write_roms_bndy(bndy_file,roms)
  
nc = netcdf(bndy_file,'write');
if isempty(nc)
  error([ 'Problem opening ' bndy_file ]);
end

% append to an existing boundary file
time_variable_name = nc{'temp_west'}.time(:);
nc_tindex = length(nc{time_variable_name}(:))+1;

nc{'ocean_time'}(nc_tindex) = roms.time;

for varlist = { 'temp','salt','u','v'}
  varname = char(varlist);
  data = getfield(roms,varname);
  nc{[varname '_west']}(nc_tindex,:,:) = data(:,:,1);
  nc{[varname '_east']}(nc_tindex,:,:) = data(:,:,end);
  nc{[varname '_north']}(nc_tindex,:,:) = data(:,end,:);
  nc{[varname '_south']}(nc_tindex,:,:) = data(:,1,:);
end

for varlist = { 'zeta','ubar','vbar'}
  varname = char(varlist);
  data = getfield(roms,varname);
  nc{[varname '_west']}(nc_tindex,:,:) = data(:,1);
  nc{[varname '_east']}(nc_tindex,:,:) = data(:,end);
  nc{[varname '_north']}(nc_tindex,:,:) = data(end,:);
  nc{[varname '_south']}(nc_tindex,:,:) = data(1,:);
end

close(nc)

