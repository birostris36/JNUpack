function [mean_noNAN] = meannan(data1); 
% 
%  function [mean_noNAN] = meannan(data1); 
% 
%  This function will calculate the mean of a series, excluding 
%  the NaN's. 
%   
%  i.e.  data1 = [1.0 2.1 NaN 0.6 1.1 NaN 0.5 0.9 1.4 2.0];  
% 
%         [x] = meannan(data1) 
% 
%        will calculate the mean of ... 
%                [1.0 2.1 ... 0.6 1.1 ... 0.5 0.9 1.4 2.0] 
% 
%           x = 1.2000 
% 
 
g1 = find(isfinite(data1)); 
 
if ~isempty(g1); 
mean_noNAN = mean(data1(g1)); 
else; 
mean_noNAN = NaN; 
end; 
 
