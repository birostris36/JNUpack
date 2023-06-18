function [Z]=roms_z(cdf)
%ROMS_Z  Returns 3D z(i,j,k) from ROMS netcdf file.


ncmex('setopts',0);
ncid=ncmex('open',cdf,'nowrite');
if(ncid==-1),
  disp(['file ' cdf ' not found'])
  return
end

% Get all the variables needed to compute z.
% The equation is z = zeta * (1 + s) + hc*s + (h - hc)*C(s)

[s_rho_dimid, status] = ncmex ( 'dimid', ncid, 's_rho' );
if ( status == -1 )
	fprintf ( 2, 'Could not get s_rho dimid from %s.\n', cdf );
	ncmex ( 'close', ncid );
	return;
end
[dimname, s_rho_length, status] = ncmex ( 'diminq', ncid, s_rho_dimid );
if ( status == -1 )
	fprintf ( 2, 'Could not get s_rho length from %s.\n', cdf );
	ncmex ( 'close', ncid );
	return;
end

[sc, status] = ncmex ( 'varget', ncid, 'sc_r', [0], [-1] );

[hc_varid, status] = ncmex ( 'varid', ncid, 'hc' );
if ( status == -1 )
	fprintf ( 2, 'Could not get hc varid from %s.\n', cdf );
	ncmex ( 'close', ncid );
	return;
end
[hc, status] = ncmex ( 'varget1', ncid, hc_varid, [0] );
if ( status == -1 )
	fprintf ( 2, 'Could not get hc from %s.\n', cdf );
	ncmex ( 'close', ncid );
	return;
end
[h, status] = ncmex ( 'varget', ncid, 'h', [0 0], [-1 -1] );
if ( status == -1 )
	fprintf ( 'scrum_zslice:  could not get ''h'' in %s.', cdf );
	return;
end
h = h';
[Cs_r, status] = ncmex ( 'varget', ncid, 'Cs_r', [0], [-1] );
if ( status == -1 )
	fprintf ( 'scrum_zslice:  could not get ''Cs_r'' in %s.', cdf );
	return;
end


% get the free surface if you want to be really accurate about it
% [zeta, status] = ncmex ( 'varget', ncid, 'zeta', [timestep 0 0], [1 -1 -1] );
% zeta = zeta';
zeta = zeros(size(h));

% Construct the depth.
n = length(sc);
i = 1;
zi = zeta * (1+sc(i)) + hc*sc(i) + (h - hc)*Cs_r(i);
z = zi;

% 
% if it looks upside-down, replace with following line 
%for i = n:-1:1
for i = 2:n
	zi = zeta * (1+sc(i)) + hc*sc(i) + (h - hc)*Cs_r(i);
	z = cat ( 3, z, zi );
end

z = permute ( z, [2 1 3] );
Z = z;
