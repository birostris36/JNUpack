function r=roms_edge_ramp(d,m)
% r = roms_edge_ramp(d,m)
%
% Linearly ramp the edges of matrix d from interior value to zero at edges
% The edge zone is m points wide (Default m = 5)

if nargin < 2
  m = 5;
end
n = floor((m+1)/2);
m = 2*n-1;

% start with ones
b = ones(size(d));

% make edges of mask b zero
% this border needs to be wide enough to remain zero after averaging
b([1:n size(b,1)+[(-n+1):0]],:) = 0;
b(:,[1:n size(b,2)+[(-n+1):0]]) = 0;

% average points to complete the mask
h = fspecial('average',2*n-1);
b = filter2(h,b);

% apply the ramp mask to the data
r = d.*b;
