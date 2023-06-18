function [x,iz,idata,izeta,dist] = rslice_vertical_interp0z(pt1,pt2,grid,gdata,zeta);
% interp1: "bilinear_interpolation''
% it was modified from John Evans's rslice_xslice.m
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
%
% John Evans, BJ Choi

x = [];
y = [];
u = [];
dist = [];

% define sampling points along a horizontal transect 

horizontal_resolution = 100;

xinc = (pt2.x - pt1.x) / horizontal_resolution;
if( abs(xinc) > 1.e-10 )
 x_transect = [pt1.x:xinc:pt2.x]';
else
 x_transect = ones(horizontal_resolution+1,1).*pt1.x;
end
transect_length = length(x_transect);

yinc = (pt2.y - pt1.y) / horizontal_resolution;
if( abs(yinc) > 1.e-10 ) 
 y_transect = [pt1.y:yinc:pt2.y]';
else
 y_transect = ones(horizontal_resolution+1,1).*pt1.y;
end


% find the closest points from the sampling points

for j = 1:length(x_transect)
	% Assume the projection is ok to do this.
	d = sqrt ( (grid.x - x_transect(j)).^2 + (grid.y - y_transect(j)).^2 );
	ind = find(d == min(d(:)) );
	closest_point_index(j,1) = ind(1);
end

[r,c] = size ( grid.x );
closest_row_index = mod ( closest_point_index - 1, r ) + 1;
closest_col_index = floor ( (closest_point_index - 1) / r ) + 1;

% Where is the mask equal to zero?  Use this information to NaN out
% land spots down the water column.
land_inds = find(grid.mask == 0);
grid.mask(land_inds) = NaN;

% Nan out land indices
zeta = zeta .* grid.mask;

nz = size(gdata,1);
mask = repmat ( grid.mask, [1 1 nz] ); %REPMAT Replicate and tile an array.
mask = permute ( mask, [3 1 2] );
data = gdata .* mask;


% Now interpolate the data.
% Now create the y array (for display purposes).  In real physical
% terms, it is actually z.
       
% switch ( gdata.interpolation_method )
%case 'nearest_neighbor'

    idata=zeros(nz,length(x_transect));
    izeta=zeros(length(x_transect),1);
    iz=zeros(nz,length(x_transect));
    zzr=zeros(nz);

    % transformation from s-coordinate to z-coordinate
    %z0r = (grid.sc_r-grid.Cs_r).*grid.hc + grid.Cs_r.*grid.h;
    %zzr = z0r + squeeze(zeta).*(1.0 + z0r./grid.h);

    for j = 1:length(x_transect)  
          jj=closest_row_index(j);
          ii=closest_col_index(j);
          izeta(j)=zeta(jj,ii);
          z0r=(grid.sc_r-grid.Cs_r).*grid.hc + grid.Cs_r.*grid.h(jj,ii);
          zzr=z0r+zeta(jj,ii).*(1.0 + z0r./grid.h(jj,ii));
          z_top=max( zzr );
          z_bottom=min( zzr );
          idz=(z_top-z_bottom)/(nz-1);
          iz(:,j)=[z_bottom:idz:z_top]';
          idata(:,j)=interp1(zzr(:),data(:,jj,ii),iz(:,j),'linear');     
    end

% Geographic coordinates.
% find the cumulative distance along the transect in
% geographic coordinates.
  r = m_lldist ( x_transect, y_transect );
  x = [0; cumsum(r)];
  x = repmat ( x', nz, 1 );
  dist = x(1,:);
