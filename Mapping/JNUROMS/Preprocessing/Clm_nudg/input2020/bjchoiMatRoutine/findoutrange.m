function index = findoutrange(x,lim,flag)
% index = findoutrange(x,lim,flag)
%
% The reverse of findinrange:
%
% Returns the indicies of elements of matrix X that fall within the  
% ranges in vector LIM, i.e.   lim(1) <= x <= lim(2)
%
% X is assumed to be M by N where the row dimension M is the number of
% data points, and the column dimension N is the number of coordinates, e.g.
% N=2 if X has columns for latitude and longitude
%
% LIM is a vector with 2*N elements being pairs of [min max] for each
% coordinate
%
% If input FLAG is present the test is lim(1) <= x < lim(2)
%                                                 ***
%
% Examples: 
%
% To find the indicies of LONGITUDE that fall in the range Xmin to Xmax:
% index = findinrange(LONGITUDE,[Xmin Xmax]);
%
% To find the indicies of LATITUDE and LONGITUDE that fall in the range 
% Xmin to Xmax, Ymin to Ymax:
% index = findinrange([LATITUDE LONGITUDE],[Ymin Ymax Xmin Xmax]);
%
% John Wilkin 17 Jul 97

% In case x was given as a row vector, convert it to a column vector so that
% the rest of the code will work 
if min(size(x))==1
  x = x(:);
end

if nargin < 3
  switch size(x,2);
    case 1
      index = find(x<=lim(1)|x>=lim(2));
    case 2
      index = find(x(:,1)<=lim(1) | x(:,1)>=lim(2) | ...
	  x(:,2)<=lim(3) | x(:,2)>=lim(4));
    case 3
      index = find(x(:,1)<=lim(1) | x(:,1)>=lim(2) | ...
	  x(:,2)<=lim(3) | x(:,2)>=lim(4) | ...
	  x(:,3)<=lim(5) | x(:,3)>=lim(6) );
  end
else
  switch size(x,2);
    case 1
      index = find(x<=lim(1)|x>lim(2));
    case 2
      index = find(x(:,1)<=lim(1) | x(:,1)>lim(2) | ...
	  x(:,2)<=lim(3) | x(:,2)>lim(4));
    case 3
      index = find(x(:,1)<=lim(1) | x(:,1)>lim(2) | ...
	  x(:,2)<=lim(3) | x(:,2)>lim(4) | ...
	  x(:,3)<=lim(5) | x(:,3)>lim(6) );
  end
end
