function grd = roms_zcalc(scoord,h);
% grd = roms_zcalc(scoord,h);
%
% Computes ROMS z coordinates from the model depth and s-coordinate 
% parameters
%
% The results are returned in a structure
%
% If the first input SCOORD is a history file name, the function attempts to
% compute Z from the parameters stored in the history file
%
% If the first input is not a file name it must a 4-element vector of
% s-coordinate parameters [theta_s theta_b Tcline N] and there must be 
% second argument specifying the model depths
%
% John Wilkin

if isstr(scoord)
  file = scoord;
  nc = netcdf(file);

  grd.h = nc{'h'}(:);
  grd.theta_s = nc{'theta_s'}(:);
  grd.theta_b = nc{'theta_b'}(:);
  grd.Tcline = nc{'Tcline'}(:);
  
  % this is ugly but I can't find a better way to get around the case of
  % having not written all the depths to the history file (like in mcc)
  thedims = dim(nc);
  grd.N = thedims{9}(:);
  
else
  grd.theta_s = scoord(1);
  grd.theta_b = scoord(2);
  grd.Tcline = scoord(3);
  grd.N = scoord(4);
  grd.h = h;
end

% calculate the z depths of the s-coordinate points
% rho-points
h = grd.h;
[z,grd.sc_r,grd.Cs_r,grd.hc] = scoord(h(:),grd.theta_s,grd.theta_b,...
    grd.Tcline,grd.N,0,1,1);
z_r = reshape(z,[size(h) grd.N]);
grd.z_r = permute(z_r,[3 1 2]); 

% w-points
[z,grd.sc_w,grd.Cs_w,grd.hc] = scoord(h(:),grd.theta_s,grd.theta_b,...
    grd.Tcline,grd.N,1,1,1);
z_w = reshape(z,[size(h) grd.N+1]); 
grd.z_w = permute(z_w,[3 1 2]); 

