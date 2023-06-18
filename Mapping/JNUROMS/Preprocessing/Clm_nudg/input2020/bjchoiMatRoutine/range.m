function tmp = range(vect)
%
% RANGE   RANGE(X) returns the minimum and maximum of X ignoring NaNs
%
tmp = [min(vect(:)) ;max(vect(:))];
