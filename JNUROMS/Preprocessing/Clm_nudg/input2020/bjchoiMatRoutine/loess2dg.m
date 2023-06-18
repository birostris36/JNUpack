function zi = loess2dg(xd,yd,dd,xi,yi,rmax,method)
% generic 2-dimensional loess filter estimate at xi,yi of data [xd,yd,dd] 
%
% zi = loess2d(xd,yd,dd,xi,yi,rmax,method)
%
% Inputs:
%   xd,yd are the coordinates
%   dd is the data
%   xi,yi are the target coordinates
%   Depending on 'method', rmax is either
%      loess filter radius in same units as xi,yi 
%   If method is:
%      'number' then rmax is the number of points to use (adaptive width)
%      'distance' then rmax is the half width of the filter (can be a vector 
%      if different scales are to be applied to x,y directions)
%
% John Wilkin

if nargin < 6
  method = 'distance';
end

switch method
  case 'number'    
  case 'distance'
    if length(rmax)==1
      % use the same scale in both dimensions
      rmax = [rmax rmax];
    end    
  otherwise
    error('invalid method option')
end

% preallocate output
zi = NaN*ones(size(xi));

% hold the indicies of valid data points - this will be clipped so that we
% only determine a smoothed value for points that were valid unsmoothed

datalocations = 1:length(dd(:));

% check for NaNs in the data
nodata = find(isnan(dd));
if ~isempty(nodata)
  dd(nodata) = [];
  xd(nodata) = [];
  yd(nodata) = [];
  datalocations(nodata) = [];
end
if isempty(dd)
  bell;disp('There are no data');return
end

% loop over index of input coordinates
for j=1:length(xi(:))
  
  if rem(j,100)==0
    disp([ 'Doing ' int2str(j) ' of ' int2str(length(xi(:)))])
  end
  
  xj = xi(j);
  yj = yi(j);

  % centre coordinates
  x = xd - xj;
  y = yd - yj;

  switch method
    
    case 'number'
      % find the nearest rmax points to the grid location
      r = abs(x+sqrt(-1)*y); 
      [tmp,Q] = sort(r);
      Q = Q(1:rmax);
      
    case 'distance'      
      % normalized distance metric
      r = abs(x/rmax(1)+sqrt(-1)*y/rmax(2)); 
      Q = find(r<=1);
  end

  r = r(Q);
  x = x(Q);
  y = y(Q);
  d = dd(Q);

  % convert to vectors
  d = d(:);
  x = x(:); 
  y = y(:); 
  r = r(:); 

  if 1
    %  if (min([length(find(x>0)) length(find(x<0))]) > 3 & ...
    %      min([length(find(y>0)) length(find(y<0))]) > 3 )
    % this checks that there are enough data to do the quadratic fit
    
    % compute local weighted least squares quadratic estimate (loess)
    
    % form matrix A of data coordinates
    A = [ones(size(x)) x y x.*y x.^2 y.^2];
    
    % calculate weights w
    w = (1-r.^3).^3;
    
    % apply weights
    A = w(:,ones([1 6])).*A;
    rhs = w.*d;
    
    % c is the vector of coefficients of the quadratic fit:
    % fit = c(1) + c(2)*x + c(3)*x^2 + c(4)*x*y +c(5)*y +c(6)*y^2
    
    % solve least-squares fit of A*c = d
    c = A\rhs;
    
    % evaluate fitted value:
    % formally, zj = [1 xj yj xj*yj xj^2 yj^2]*c, but since we moved
    % origin to xj the answer is just c(1)
    zj = c(1);
    
    zi(j) = zj;

  end
  
end

