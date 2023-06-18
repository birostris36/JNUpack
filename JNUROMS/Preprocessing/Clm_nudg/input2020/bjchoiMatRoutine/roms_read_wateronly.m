function f = roms_read_wateronly(fname,vname,grd,tindex,kindex);
% really ugly hack of Hernan's nc_read
% John Wilkin

%----------------------------------------------------------------------------
% Inquire information from NetCDF file.
%----------------------------------------------------------------------------

% Inquire about file dimensions.

[dnames,dsizes]=nc_dim(fname);

for n=1:length(dsizes),
  name=deblank(dnames(n,:));
  switch name
    case 'xi_rho',
      Lr=dsizes(n);
    case 'xi_u',
      Lu=dsizes(n);
    case 'xi_v',
      Lv=dsizes(n);
    case 'eta_rho',
      Mr=dsizes(n);
    case 'eta_u',
      Mu=dsizes(n);
    case 'eta_v',
      Mv=dsizes(n);
    case 's_rho',
      Nr=dsizes(n);
    case 's_w',
      Nw=dsizes(n);
  end,
end,  
 
% Inquire about requested variable.

[vdnames,vdsizes,igrid]=nc_vinfo(fname,vname);

% Check if data is only available at water points.

is2d=0;
is3d=0;
water=0;
if (~isempty(vdsizes)),
  for n=1:length(vdsizes),
    name=deblank(vdnames(n,:));
    switch name
      case 'xy_rho',
        msknam='_rho';
        is2d=1; Im=Lr; Jm=Mr;
      case 'xy_u',
        msknam='_u';
        is2d=1; Im=Lu; Jm=Mu;
      case 'xy_v',
        msknam='_v';
        is2d=1; Im=Lv; Jm=Mv;
      case 'xyz_rho',
        msknam='_rho';
        is3d=1; Im=Lr; Jm=Mr; Km=Nr;
      case 'xyz_u',
        msknam='_u';
        is3d=1; Im=Lu; Jm=Mu; Km=Nr;
      case 'xyz_v',    
        msknam='_v';
        is3d=1; Im=Lv; Jm=Mv; Km=Nr;
      case 'xyz_w',    
        msknam='_rho';
        is3d=1; Im=Lr; Jm=Mr; Km=Nw;
      otherwise
	msknam = [];
    end,    
    lonnam = ['lon' msknam];
    latnam = ['lat' msknam];
    msknam = ['mask' msknam];
  end,
end,
water=is2d | is3d;

% If water data only, read in Land/Sea mask.

if (water)
  mask = getfield(grd,msknam)';
  lon = getfield(grd,lonnam)';
  lat = getfield(grd,latnam)';
end

%----------------------------------------------------------------------------
% Read in requested variable.
%----------------------------------------------------------------------------

% keyboard

if (water),

  % possibly do something here to avoid reading the whole record, but I
  % don't think it can be done.
  %
  %  if nargin < 6
  % ind = find(mask==1);
  %  else
  %    ind = find(mask==1 & ...
  %	lon>=ax(1) & lon<=ax(2) & lat>=ax(3) & lat<=ax(4) );
  %  end
  
  ind = find(mask==1);
  data = ones(size(mask));
  if is2d
    kindex = 1;
    len = vdsizes(2);
  else
    len = vdsizes(2)/Nr;    
  end  
  nc = netcdf(fname);
  wet_data = nc{vname}(tindex,(kindex-1)*len+(1:len));  
  data(ind) = wet_data;
  f = data'; % switch back to same orientation as in grd
  close(nc)
  
else

  if (nargin < 3),
    f=ncread(fname,vname);
  else
    f=ncread(fname,vname,tindex);
  end,
  
end



