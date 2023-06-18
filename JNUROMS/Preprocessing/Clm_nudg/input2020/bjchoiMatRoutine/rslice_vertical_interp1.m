function [x,iz,idata,dist] = rslice_vertical_interp1(pt1,pt2,grid,gdata);
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


% find 4 closest points from the sampling points

for j = 1:length(x_transect)
	% Assume the projection is ok to do this.
	d = sqrt ( (grid.x - x_transect(j)).^2 + (grid.y - y_transect(j)).^2 );
        d_temp = d;
        ind=[];

  
 	ind_temp = find( d == min(d_temp(:)) );

        if( grid.mask(ind_temp(1) ) < 1 )

          ind = [ind ind_temp(1)];
         
        else

          % find the first closest point 
          d_temp( ind_temp ) = 1.e20;
          ind = [ind ind_temp];

          % find the second, third and fourth closest points
          for kk=2:4
 	   ind_temp = find( d == min(d_temp(:)) );
           d_temp( ind_temp ) = 1.e20;
           if( grid.mask(ind_temp) > 0 )
            ind = [ind ind_temp];
           end
          end

        end  

        Npts=length(ind);
        if( Npts < 1 )
 	 ind_temp = find( d == min(d(:)) );
         ind(1:4)=ind_temp;
        elseif( Npts < 4 )
         ind_temp = ind(1);
         ind(Npts+1:4)=ind_temp;
        end
        closest_point_index(j,1:4) = ind(1:4);
end

[r,c] = size ( grid.x );
closest_row_index = mod ( closest_point_index - 1, r ) + 1;
closest_col_index = floor ( (closest_point_index - 1) / r ) + 1;

% calculate linear weights based on distance between points

for j = 1:length(x_transect)
     xx0=x_transect(j);
     yy0=y_transect(j); 
     for kk=1:4
         jj=closest_row_index(j,kk);
         ii=closest_col_index(j,kk);
         xx=grid.x(jj,ii);
         yy=grid.y(jj,ii);
         dis(kk)=m_lldist([xx0  xx],[yy0 yy]);
     end
     sum_dis=sum( dis(1:4) );
     closest_point_weight(j,1:4) = dis(1:4)./sum_dis;
end

% Where is the mask equal to zero?  Use this information to NaN out
% land spots down the water column.
land_inds = find(grid.mask == 0);
grid.mask(land_inds) = NaN;

% Nan out land indices
nz = size(gdata,1);
mask = repmat ( grid.mask, [1 1 nz] ); %REPMAT Replicate and tile an array.
mask = permute ( mask, [3 1 2] );
data = gdata .* mask;


% Now interpolate the data.
% Now create the y array (for display purposes).  In real physical
% terms, it is actually z.
       
% switch ( gdata.interpolation_method )
% case 'bilinear_interpolation'
	
    idata=zeros(nz,length(x_transect));
    iz=zeros(nz,length(x_transect));
    z=zeros(nz,1);
   %zeta=0;
    for j = 1:length(x_transect)  
      for kk=1:4  
          jj=closest_row_index(j,kk);
          ii=closest_col_index(j,kk);
          idata(:,j)=idata(:,j) + data(:,jj,ii).*closest_point_weight(j,kk);
          z(:)=grid.z_r(:,jj,ii)*closest_point_weight(j,kk);
          iz(:,j)   =iz(:,j)    + z(:).*closest_point_weight(j,kk);
      end
    end

% Geographic coordinates.
% find the cumulative distance along the transect in
% geographic coordinates.
  r = m_lldist ( x_transect, y_transect );
  x = [0; cumsum(r)];
  x = repmat ( x', nz, 1 );
  dist = x(1,:);
