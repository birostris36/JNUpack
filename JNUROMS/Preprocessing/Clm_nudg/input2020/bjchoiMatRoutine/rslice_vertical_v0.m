function [x,z,z_yx_data,dist] = rslice_xslice (pt1, pt2, grid, arg )
% RSLICE_XSLICE:  arbitrary vertical slicing for ROMS output
%
% PARAMETERS:
% Input:
%    pt1, pt2:
%        Structures with two fields, x and y.  This could either
%        be in geographic coordinates or a projected coordinate
%        system.
%     grid:
%        Structure with four fields, x, y, z, and mask.
% Output:

x = [];
y = [];
u = [];
dist = [];

if nargin == 4
	gdata = arg;
end

if gdata.verbose > 0
	rslice_message ( sprintf ( '%s:  pt1 = [%f %f]', mfilename, pt1.x, pt1.y ) );
	rslice_message ( sprintf ( '%s:  pt2 = [%f %f]', mfilename, pt2.x, pt2.y ) );
end

if isempty ( gdata.z_type )
	msg = sprintf ( '%s:  variable %s does not have a depth dimension\n', mfilename, gdata.varname );
	rslice_warning ( msg );
	return
end

horizontal_resolution = 100;

xinc = (pt2.x - pt1.x) / horizontal_resolution;
x_transect = [pt1.x:xinc:pt2.x]';
transect_length = length(x_transect);

yinc = (pt2.y - pt1.y) / horizontal_resolution;
y_transect = [pt1.y:yinc:pt2.y]';



for j = 1:length(x_transect)


	%
	% Assume the projection is ok to do this.
	d = sqrt ( (grid.x - x_transect(j)).^2 + (grid.y - y_transect(j)).^2 );
	ind = find(d == min(d(:)) );
	closest_point_index(j,1) = ind(1);

end

[r,c] = size ( grid.x );
closest_row_index = mod ( closest_point_index - 1, r ) + 1;
closest_col_index = floor ( (closest_point_index - 1) / r ) + 1;

%
% Where is the mask equal to zero?  Use this information to NaN out
% land spots down the water column.
land_inds = find(grid.mask == 0);
grid.mask(land_inds) = NaN;


%
% retrieve the data at the current timestep
data = rslice_retrieve_data ( gdata );

%
% Nan out land indices
nz = size(data,1);
mask = repmat ( grid.mask, [1 1 nz] );
mask = permute ( mask, [3 1 2] );
data = data .* mask;

%
% Now interpolate the data.
switch ( gdata.interpolation_method )
case 'nearest_neighbor'
	%
	% Just turn the closest point 2D indices into 3D indices.
	% This extracts out the entire water column at a given point,
	% and does it for each point.
	nz = size(data,1);
	xy_inds = repmat ( closest_point_index' - 1, nz, 1 ) * nz + 1;

	%
	% xy_inds gives a 2D array where each column is the same value
	% and points to the same location in the 1st z hyperslab.  By
	% adding [0:nz-1] to that column, we then get the entire water column.
	%

	offset_inds = repmat ( [0:nz-1]', 1, transect_length );

	inds = xy_inds + offset_inds;

	z_yx_data = data(inds);




	%
	% Now reshape it into a 2D array.
	z_yx_data = reshape ( z_yx_data, nz, transect_length );
	z_yx_data = flipud ( z_yx_data );


	%
	% Now create the X array.
	if isempty ( gdata.proj ) 

		%
		% Geographic coordinates.
		% find the cumulative distance along the transect in
		% geographic coordinates.
		r = m_lldist ( x_transect, y_transect );
		x = [0; cumsum(r)];

	else
		
		%
		% Assume the projection is ok to do this.
		%transect_resolution = sqrt ( xinc^2 + yinc^2 );
		%x = transect_resolution * [0:transect_length-1];
		r = m_xydist ( x_transect, y_transect );
		x = [0; cumsum(r)];

	end

	x = repmat ( x', nz, 1 );

	dist = x(1,:);

	%
	% Now create the y array (for display purposes).  In real physical
	% terms, it is actually z.
	[z,sc,Cs]=scoord_xy( gdata.bathymetry, gdata.theta_s, gdata.theta_b, ...
	                     gdata.Tcline, length(gdata.sc_r), gdata.z_type, ...
			     closest_point_index );


    varname = gdata.varname;
    save slice_data.mat x z z_yx_data varname;

otherwise
	msg = sprintf ( '%s:  unhandled interpolation method.\n' );
	rslice_error ( msg );
end




