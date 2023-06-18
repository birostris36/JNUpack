function data=ext_data_OGCM_yecs(nc,X,Y,vname,tndx,lon,lat,k,Roa,interp_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extrapole one horizontal ECCO (or Data) slice on a ROMS grid
%
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Contributions of P. Marchesiello (IRD) and J. Lefevre (IRD)
%
%  Updated    6-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% zeta=ext_data_OGCM_yecs(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method);
% dnum = datenum([2017,04,01])-datenum([2017,01,01]) + 1;
% nc=netcdf([OGCM_dir,OGCM_prefix,num2str(dnum,'%4.4i'),'.nc']);
% X = lonT;
% Y = latT;
% vname = 'zeta';
% tndx = 1; %time
% lon = lon;
% lat = lat;
% k = 16; % s - level?
% Roa = 0;
% interp_method =  'cubic'; 


%  lonU
%  latU
%  vname = 'u_eastward';
%  tndx = 1;
%  lon
%  lat
%  k

% ubar1d=ext_data_OGCM_yecs(
% nc
% X = lonU;
% Y = latU;
% vname = 'ubar';
% tin = 1;
% lon = squeeze(lon(1:2,:));
% lat = squeeze(lat(1:2,:));
% k = 1;
% Roa = 0;
% interp_method = 'cubic';

%
% extrapolation parameters
%
default=0;
if strcmp(vname,'SAVE') | strcmp(vname,'salt')
  default=34.6;
end
%
% Get the ROMS grid extension + a little margin (~ 2 data grid points)
%
dx=max(abs(gradient(X)));
dy=max(abs(gradient(Y)));
dl=2*max([dx dy]);
%
lonmin=min(min(lon))-dl;
lonmax=max(max(lon))+dl;
latmin=min(min(lat))-dl;
latmax=max(max(lat))+dl;
%
% Extract a data subgrid
%
% j=find(Y>=latmin & Y<=latmax);
% i1=find(X-360>=lonmin & X-360<=lonmax);
% i2=find(X>=lonmin & X<=lonmax);
% i3=find(X+360>=lonmin & X+360<=lonmax);

aa1 = (X - lonmin).^2 + (Y - latmax).^2;
% ia1 = find(aa1>0);
[j1,i1] = find(aa1 == min(min(aa1)));
aa2 = (X - lonmax).^2 + (Y - latmin).^2;
% ia2 = find(aa2>0);
[j2,i2] = find(aa2 == min(min(aa2)));


% if ~isempty(i2)
%   x=X(i2);
% else
%   x=[];
% end
% if ~isempty(i1)
%   x=cat(2,X(i1)-360,x);
% end
% if ~isempty(i3)
%   x=cat(2,x,X(i3)+360);
% end

% y=Y(j);

%
%  Get dimensions
%
% ndims=length(dim(nc{vname}));
%
% Get data (Horizontal 2D matrix)
%
% iid = find(X>=lonmin & X<=lonmax & Y>=latmin & Y<=latmax);

% a1 = find(j == i2);

% if ~isempty(i2)
%   if ndims==2
%     data=squeeze(nc{vname}(j,i2));
%   elseif ndims==3
%     data=squeeze(nc{vname}(tndx,j,i2));
%     
%   elseif ndims==4
%     data=squeeze(nc{vname}(tndx,k,j,i2));
%   else
%     error(['Bad dimension number ',num2str(ndims)])
%   end
% else
%   data=[];
% end
% if ~isempty(i1)
%   if ndims==2
%     data=cat(2,squeeze(nc{vname}(j,i1)),data);
%   elseif ndims==3
%     data=cat(2,squeeze(nc{vname}(tndx,j,i1)),data);
%   elseif ndims==4
%     data=cat(2,squeeze(nc{vname}(tndx,k,j,i1)),data);
%   else
%     error(['Bad dimension number ',num2str(ndims)])
%   end
% end
% if ~isempty(i3)
%   if ndims==2
%     data=cat(2,data,squeeze(nc{vname}(j,i3)));
%   elseif ndims==3
%     data=cat(2,data,squeeze(nc{vname}(tndx,j,i3)));
%   elseif ndims==4
%     data=cat(2,data,squeeze(nc{vname}(tndx,k,j,i3)));
%   else
%     error(['Bad dimension number ',num2str(ndims)])
%   end
% end

I1 = i1-20;
I2 = i2+20;
J2 = j2-20;
J1 = j1+20;

%
%  Get dimensions
%

ndims=length(dim(nc{vname}));
%
% Get data (Horizontal 2D matrix)
%

if ndims==2
    data=squeeze(nc{vname}(J2:J1,I1:I2));
elseif ndims==3
    data=squeeze(nc{vname}(tndx,J2:J1,I1:I2));
elseif ndims==4
    data=squeeze(nc{vname}(tndx,k,J2:J1,I1:I2));
else
    error(['Bad dimension number ',num2str(ndims)])
end
% if isequal(vname ,'zeta'|'ubar'|'vbar')
% %     disp('yes')
% data=squeeze(nc{vname}(tndx,J2:J1,I1:I2)); % scha
% else
% %     disp('yes1')
% data=squeeze(nc{vname}(tndx,k,J2:J1,I1:I2)); % scha
% end
x = X(J2:J1,I1:I2);
y = Y(J2:J1,I1:I2);

% min(min(lonT(J2:J1,I1:I2)))
% max(max(lonT(J2:J1,I1:I2)))
% min(min(latT(J2:J1,I1:I2)))
% max(max(latT(J2:J1,I1:I2)))
%
% Perform the extrapolation
%
% [data,interp_flag]=get_missing_val(x,y,data,NaN,Roa,default);
data(data>10000) = nan;
% [data,interp_flag]=get_missing_val(x,y,data,NaN,Roa,default);

idn = find(~isnan(data));
%
% Interpolation on the ROMS grid
%
data=griddata(x(idn),y(idn),data(idn),lon,lat,interp_method);
% if interp_flag==0
%   data1=interp2(x,y,data,lon,lat,'nearest');
% else
%   data1=interp2(x,y,data,lon,lat,interp_method);
% end
%
return
