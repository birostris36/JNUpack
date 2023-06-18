function write_roms_frc(frc_file,roms)
  
nc = netcdf(frc_file,'write');
if isempty(nc)
  error([ 'Problem opening ' frc_file ]);
end

% the variable name must be passed with the roms data structure
theVarname = roms.frcvarname;

% append to an existing file
time_variable = nc{theVarname}.time(:);
nc_tindex = length(nc{time_variable}(:))+1;

% enter the values
nc{time_variable}(nc_tindex) = roms.time;
data = getfield(roms,theVarname);
nc{theVarname}(nc_tindex,:,:) = data;

close(nc)

if any(isnan(data))
  disp(' ')
  warning([ 'There were NaNs in the data passed to write_roms_frc'])
end	

    
    
