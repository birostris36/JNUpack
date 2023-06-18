function [sum_noNAN] = sumnan(data1);
%
%  function [sum_noNAN] = sumnan(data1);
%
%  This function will calculate the sum of a series, excluding
%  the NaN's.
%  
%  i.e.  data1 = [1.0 2.1 NaN 0.6 1.1 NaN 0.5 0.9 1.4 2.0]; 
%
%         [x] = sumnan(data1)
%
%        will calculate the sum of ...
%                [1.0 2.1 ... 0.6 1.1 ... 0.5 0.9 1.4 2.0]
%
%           x = 9.6000
%

g1 = find(isfinite(data1));

sum_noNAN = sum(data1(g1));

