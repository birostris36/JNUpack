function coast_points_index = roms_find_coastal_points(grd)
% Find roms grid points that are adjacent to masked land so 
% they can be made uniform depth if desired
% 
% Usage:
%   coast_points_index = roms_find_coastal_points(grd)
% 
% where 
%   grd is a grd structure (from roms_get_grid)
%
% Use ind2sub to convert these to i,j indices if required.
%
% John Wilkin

% land/sea mask
m = grd.mask_rho_nan;

% convolve
k = ones(3,3);
tmp = conv2(m,k,'same');

% tmp is now NaN at any point that had land for a neighbour
% The zero padding of the convolution causes the correct behaviour 
% at the boundary 

% now change to zero the values that were land to start 
land = find(isnan(m)==1);
tmp(land) = 0;

% this leaves only the grid points that are *adjacent* to masked 
% land, but are water
coast_points_index = find(isnan(tmp)==1);

