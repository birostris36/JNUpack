function make_inti_fromNCOMoutput
% ==========================================================================
% build INITIAL fields for a nested ROMS grid
% from data in a larger NCOM grid.
% Data file could be history or average file.
%
% Algorithm: 
% Get vertical distribution of a variable for a given grid cell,
% (j,i) nested grid cell.
% (1) find nearest four points in larger domain
% (2) calculate weight using distance
% (3) vertical interpolation at four points
% (4) horizontal interpolation base on wights in (2)
% BJ Choi, sept-15-2005
% BJ Choi, may-23-2005
%
% in_file      = larger domain roms data file (NCOM output)
% analysis_day = time of in_file (in days)
% out_file     = initial file for the nested domain (init.nc)
% base_date = [2000 1 1 0 0 0];
%
% USE: nc_create_roms_init.m
%      write_roms_init.m
%      grd.m
%      dailymean_ncom_grid.m
%      dailymean_ncom_output.m yearday.m
%
% ==========================================================================
% version 3 (Oct-07-2005)
% 
% (1) u and v are not true eastward and northward velocities. They are ROMS
% variables. Interpolate u and v on rho points and rotate velocity vectors
% to true eastward and northward.
%
% (2) volume conservation across open boundaries. 
% Because bottom depths are different in larger domain 
% and the nested domain, current speed across the boundary (ubar, vbar) should be 
% multiplied by depth ratio, i.e., U_nested = U_large * h_large / h_nested 
% Here, h_large / h_nested = 1 + (h_l - h_n) / h_n ~= 1 + alpha * (h_l - h_n) / h_n 
%       alpha = 0.30 or 1.0 
% Calculate amount of volume transport difference between large domain data and 
% nested domain estimated value, and add a constant speed to the nested domain 
% estimated speed for the compensation of total volume transport difference
%
% ==========================================================================
clear 
close all

disp(' read grid information of the nested smaller domain')

if exist('gn')~=1
    gn = grd('ctz3km');     % (2) Get Nested grid information
end

% (4) set time: time from the base_date in days.
% May  31, 2000 = 152
% June 30, 2000 = 182
% July 31, 2000 = 213
% Aug  31, 2000 = 244
% Sep   1, 2000 = 245
% analysis_day=152; % June  1, 2005 0 0 0 = 152

 analysis_year=2002;% year
 analysis_day=1;    % January 1, 2002
 analysis_day=90;   % April 1, 2002

 % analysis_year and analysis_day 
 % will be used in dailymean_ncom_output.m

% (5) output file, initial file for the nested grid

 out_file = 'ctz3km_init2002Jan01.nc';

% ================= End of Input Parameters =================================

% Information for the initial conditions file
% and create an empty initial file

 disp([' create an empty intial file = ' out_file])
 noclobber = 0;
 base_date = [2000 1 1 0 0 0];
 time_variable = 'ocean_time';
 time_variable_units = basedate2str(base_date);
 roms.grd = gn;
 grd_file = gn.grd_file;
 sourcestr = [ 'From NCOM-CCS: 9km resolution '];
 details = [ 'Quick initial file from NCOM-CCS ' ...
             'by script ' which(mfilename) ];
 donuts = 0;
 nc_create_roms_init  % create an nc file, write variables and close it.
 disp([' creating an empty initial file done ..... '])
 disp('  ')


% read grid information from NCOM output

  disp([' read grid information of the larger domain '])

  xmin=min(min(gn.lon_rho));
  xmax=max(max(gn.lon_rho));
  ymin=min(min(gn.lat_rho));
  ymax=max(max(gn.lat_rho));

  ncom.bad_val = -1.e33;

  dailymean_ncom_grid        % make ncom grid information
                             % ncom.grd
                             % ncom.grd.z (40x84x47)
                             % ncom.grd.bottom
                             % ncom.grd.longitude
                             % ncom.grd.latitude
                             % ncom.grd.lon_rho ncom.grd.lat_rho (84x47)

% read data from larger domain



  disp([' read data from the larger domain NCOM file '])
  % analysis_year + analysis_day ==>
  dailymean_ncom_output      % calculate daily mean ncom output 
                             % on analysis_year analysis_day
                             % ncom.u  ncom.ubar
                             % ncom.v  ncom.vbar
                             % ncom.zeta
                             % ncom.temp 
                             % ncom.salt


% extraplolate data (variables) horizonally
% into the land. 
% three land grid points next to the coast will have data.
% we assume that land boundary is located only in the east.

disp([' horizontal extrapolation of original data '])

 ncom_mask=ncom.grd.mask;
 ncom_z=ncom.grd.z;

 ncom_temp = ncom.temp;
 ncom_salt = ncom.salt;
 ncom_zeta = zeros( size(ncom.zeta) );
 ncom_ubar=ncom.ubar;
 ncom_vbar=ncom.vbar;
 ncom_u=ncom.u;
 ncom_v=ncom.v;

 mask_temp=ncom_mask;
 [r,c] = size ( ncom_mask );

  for numofextrapol=1:3

   Iland=find(  mask_temp == 0 );
   num_land_grid=length(Iland);
   for i=1:num_land_grid
   
     ind = Iland(i);
     row_index = mod ( ind - 1, r ) + 1;
     col_index = floor( (ind - 1) / r ) + 1;
     extflag=0;
  
     if(     (col_index > 2) && (ncom_mask(row_index,col_index-1) == 1) )
       oj=row_index; oi=col_index-1; extflag=1;
     elseif( (col_index < c) && (ncom_mask(row_index,col_index+1) == 1) )
       oj=row_index; oi=col_index+1; extflag=1;
     elseif( (row_index > 2) && (ncom_mask(row_index-1,col_index) == 1) )
       oj=row_index-1; oi=col_index; extflag=1;
     elseif( (row_index < r) && (ncom_mask(row_index+1,col_index) == 1) )
       oj=row_index+1; oi=col_index; extflag=1;
     end

     if( extflag )
       % 2D variables
       ncom_zeta(row_index,col_index)  = ncom_zeta(oj,oi);
       ncom_ubar(row_index,col_index)  = ncom_ubar(oj,oi);
       ncom_vbar(row_index,col_index)  = ncom_vbar(oj,oi);
       % 3D variables
       ncom_temp(:,row_index,col_index)   = ncom_temp(:,oj,oi);
       ncom_salt(:,row_index,col_index)   = ncom_salt(:,oj,oi);
       ncom_u(:,row_index,col_index) = ncom_u(:,oj,oi);
       ncom_v(:,row_index,col_index) = ncom_v(:,oj,oi);
       ncom_z(:,row_index,col_index) =ncom_z(:,oj,oi);
       % reset mask value
       mask_temp(row_index,col_index)=1; 
     end

   end % of for i=1:num_land_grid
   ncom_mask=mask_temp; 
  end  % of for numofextrapol=1:3



% two-dimensional variables

  zeta      = ncom_zeta;
  true_ubar = ncom_ubar;
  true_vbar = ncom_vbar;



% ======================================================
% objective analysis in x,y,z space --> make 3-dimensional structure
% ======================================================

  % ncom_z              std_depth
  % ncom_temp           temp
  % ncom_salt      ---> salt
  % ncom_u              oa_u
  % ncom_v              oa_v

  % objective analysis vertical depth
  
  vertical_depth = [0.5:1:11 12:2:20 25:5:50 60:10:100 125:25:250 300:50:500 ...
                    600:100:1000 1200:200:2000 2500:500:5000];
  vertical_depth = -vertical_depth;

  temp = zeros( length(vertical_depth), size(ncom_mask,1), size(ncom_mask,2) )*NaN; 
  salt = temp;
  oa_u = temp;
  oa_v = temp;

  for j=1:r
  for i=1:c

      if( ncom_mask(j,i) > 0 )
      temp(:,j,i)  = interp1( ncom_z(:,j,i), ncom_temp(:,j,i), vertical_depth );
      salt(:,j,i)  = interp1( ncom_z(:,j,i), ncom_salt(:,j,i), vertical_depth );
      oa_u(:,j,i)  = interp1( ncom_z(:,j,i), ncom_u(:,j,i), vertical_depth );
      oa_v(:,j,i)  = interp1( ncom_z(:,j,i), ncom_v(:,j,i), vertical_depth );
      end

      if ( ~finite(temp(1,j,i)) ) % top level at -0.5 meter.
       temp(:,j,i)  = ncom_temp(1,j,i);
       salt(:,j,i)  = ncom_salt(1,j,i);
       oa_u(:,j,i)  = ncom_u(1,j,i);
       oa_v(:,j,i)  = ncom_v(1,j,i);
      end
   
  end
  end

  figure; pcolorjw( squeeze( temp(1,:,:) ) ); colorbar
  title([' OA result, Temperature at ' num2str( vertical_depth(1) ) ' m'])
  figure; pcolorjw( squeeze( temp(2,:,:) ) ); colorbar
  title([' OA result, Temperature at ' num2str( vertical_depth(2) ) ' m'])


  % important procedure
  % extrapolate temp, salt, u and v in horizontally direction along z-coordinate
  % over the land mass

  for k=1:length(vertical_depth)

    for j=1:r
   
     Ip = find(  finite( temp(k,j,:) ) );
     Iq = find( ~finite( temp(k,j,:) ) );
     if( length(Iq) > 1 & Iq(1)-1 > 1 )
      temp(k,j,Iq) = temp(k,j,Iq(1)-1);
     end

     Ip = find(  finite( salt(k,j,:) ) );
     Iq = find( ~finite( salt(k,j,:) ) );
     if( length(Iq) > 1 & Iq(1)-1 > 1 )
      salt(k,j,Iq) = salt(k,j,Iq(1)-1);
     end

     Ip = find(  finite( oa_u(k,j,:) ) );
     Iq = find( ~finite( oa_u(k,j,:) ) );
     if( length(Iq) > 1 & Iq(1)-1 > 1 )
      oa_u(k,j,Iq) = oa_u(k,j,Iq(1)-1);
      oa_v(k,j,Iq) = oa_v(k,j,Iq(1)-1);
     end

    end

  end

  figure; pcolorjw( squeeze( temp(46,:,:) ) ); colorbar
  title([' OA result, Temperature at ' num2str( vertical_depth(46) ) ' m'])


% ======================================================
% interpolateion from x,y,z space to ROMS grid.
% ======================================================



% size of grids
[r,c] = size ( ncom.grd.mask );
   nl = length(vertical_depth);  
[M N] = size ( gn.lon_rho );
   nn = gn.N;

% find deepest depth
maxdepth=max([max(max(ncom.grd.bottom)) max(max(gn.h)) max(-vertical_depth) ])+500;

% fill in NaN values in vertically direction
% nl, r, c
%
 
  for j=1:r
  for i=1:c

     Ip = find(  finite( temp(:,j,i) ) );
     Iq = find( ~finite( temp(:,j,i) ) );
     if ( length(Iq) > 1 & ( Iq(1)-1) > 0 )
      temp(Iq,j,i)=temp( Iq(1)-1,j,i );
     elseif(  length(Ip) < 1 & ncom.grd.mask(j,i) > 0 )
      error([' temp: error at ', num2str(j), '  ' ,num2str(j) ])
     end 

     Ip = find(  finite( salt(:,j,i) ) );
     Iq = find( ~finite( salt(:,j,i) ) );
     if ( length(Iq) > 1 & ( Iq(1)-1) > 0 )
      salt(Iq,j,i)=salt( Iq(1)-1,j,i );
     elseif(  length(Ip) < 1 & ncom.grd.mask(j,i) > 0 )
      error([' salt: error at ', num2str(j), '  ' ,num2str(j) ])
     end 

     Ip = find(  finite( oa_u(:,j,i) ) );
     Iq = find( ~finite( oa_u(:,j,i) ) );
     if ( length(Iq) > 1 & ( Iq(1)-1) > 0 )
      oa_u(Iq,j,i)=oa_u( Iq(1)-1,j,i );
     elseif(  length(Ip) < 1 & ncom.grd.mask(j,i) > 0 )
      error([' u: error at ', num2str(j), '  ' ,num2str(j) ])
     end 

     Ip = find(  finite( oa_v(:,j,i) ) );
     Iq = find( ~finite( oa_v(:,j,i) ) );
     if ( length(Iq) > 1 &  ( Iq(1)-1) > 0 )
      oa_v(Iq,j,i)=oa_v( Iq(1)-1,j,i );
     elseif(  length(Ip) < 1 & ncom.grd.mask(j,i) > 0 )
      error([' v: error at ', num2str(j), '  ' ,num2str(j) ])
     end 

  end
  end
 
% find land mask from gl grid
ocean=ones(r,c);
land =ones(r,c)*1.e20;
Isea=find( ncom_mask > 0);
land(Isea)=ocean(Isea);    % 1 for ocean and 1.e20 for land
clear ocean Isea

% extrapolate data in vertical direction
disp([' vertical extrapolation of original data '])

extsalt(1,:,:)=salt(1,:,:);         % add bottom ( at -maxdepth)
extsalt(2:nl+1,:,:)=salt(1:nl,:,:); % data
extsalt(nl+2,:,:)=salt(nl,:,:);     % add top    ( 20 m above sea level)

exttemp(1,:,:)=temp(1,:,:);         % add bottom ( at -maxdepth)
exttemp(2:nl+1,:,:)=temp(1:nl,:,:); % data
exttemp(nl+2,:,:)=temp(nl,:,:);     % add top    ( 20 m above sea level)

extu(1,:,:)=oa_u(1,:,:);         % add bottom ( at -maxdepth)
extu(2:nl+1,:,:)=oa_u(1:nl,:,:); % data
extu(nl+2,:,:)=oa_u(nl,:,:);     % add top    ( 20 m above sea level)

extv(1,:,:)=oa_v(1,:,:);         % add bottom ( at -maxdepth)
extv(2:nl+1,:,:)=oa_v(1:nl,:,:); % data
extv(nl+2,:,:)=oa_v(nl,:,:);     % add top    ( 20 m above sea level)


% initailize inital data file
 disp([' initailize roms structure (inital data file) '])

 roms.time = yearday(analysis_year,analysis_day)-datenum(base_date);
 roms.temp = zeros([gn.N size(gn.lon_rho)]);
 roms.salt = zeros([gn.N size(gn.lon_rho)]);
 roms.zeta = zeros(size(gn.h));
 roms.u = zeros([gn.N size(gn.lon_u)]);
 roms.v = zeros([gn.N size(gn.lon_v)]);
 roms.vbar = zeros(size(gn.lon_v));
 roms.ubar = zeros(size(gn.lon_u));

 ntrue_ubar = zeros(size(gn.h));
 ntrue_vbar = zeros(size(gn.h));
 ntrue_u = zeros([gn.N size(gn.lon_rho)]);
 ntrue_v = zeros([gn.N size(gn.lon_rho)]);

 write_roms_init(out_file,roms) % append zeros to an existing empty initial file
                               % and close it.
 disp([' appending default values to initial file done ..... '])
 disp('  ')

% =======================================================================
% vertical and horizonatal interpolation of variables 
% such as temp, salt, zeta, true_u, true_v, true_ubar and true_vbar
% from larger domain to the nested domain.
% =======================================================================

disp([' ================================================ '])
disp([' interpolating temp and salt on the nested domain '])
disp([' wait .....                              '])
disp(['                                         '])

for i=1:N        % grid point in the nested (smaller) domain
for j=1:M

 if ( gn.mask_rho(j,i) > 0 ) % sea; we works on a cell under water

   % find 4 nearest points
   % ind --> index of the 4 points
   % Assume the projection is ok to do this.
   d = sqrt ( (ncom.grd.lon_rho - gn.lon_rho(j,i)).^2 + (ncom.grd.lat_rho - gn.lat_rho(j,i)).^2 ); 
   d = d.*land;
   d_temp = d;
   ind=[];
   while length( ind ) < 4
      ind_temp = find( d == min(d_temp(:)) );
      ind = [ind ind_temp(1)];   
      d_temp( ind_temp(1) ) = 1.e20; 
   end

   % closest points row and column indice
   row_index = mod ( ind - 1, r ) + 1;
   col_index = floor( (ind - 1) / r ) + 1;

   % calculate linear weights based on distance between points
   xx0=gn.lon_rho(j,i); 
   yy0=gn.lat_rho(j,i); 
   for kk=1:4
    jj=row_index(kk);
    ii=col_index(kk);
    xx=ncom.grd.lon_rho(jj,ii);
    yy=ncom.grd.lat_rho(jj,ii);
    dis(kk)=m_lldist([xx0  xx],[yy0 yy]);
   end
   sum_dis=sum( dis(1:4) );
   weight(1:4) = dis(1:4)./sum_dis;

   % transformation from s-coordinate to z-coordinate
   %z0r = (grid.sc_r-grid.Cs_r).*grid.hc + grid.Cs_r.*grid.h;
   %zzr = z0r + squeeze(zeta).*(1.0 + z0r./grid.h);

   % vertical interpolation and extrapolation 
   % interpolate zeta horizontally
   izeta=0;              % interpolated sea surface, zeta
   ihl=0;                % interpolated depth of water, h
   for kk=1:4  
    jj=row_index(kk);
    ii=col_index(kk);
    zzr=vertical_depth;
    izeta=izeta+zeta(jj,ii).*weight(kk);
    ihl  =ihl  +ncom.grd.bottom(jj,ii).*weight(kk);
    ntrue_ubar(j,i)=ntrue_ubar(j,i)+true_ubar(jj,ii).*weight(kk);
    ntrue_vbar(j,i)=ntrue_vbar(j,i)+true_vbar(jj,ii).*weight(kk);
   end

   % apply volume flux conservation across open boundary (vertical factor)
   % vfactor = 1 if ihl = gn.h(j,i)
   vfactor=(ihl+izeta)/(gn.h(j,i)+izeta);    
   ntrue_ubar(j,i)=ntrue_ubar(j,i) * vfactor ;
   ntrue_vbar(j,i)=ntrue_vbar(j,i) * vfactor ;

   iz0r=(gn.sc_r-gn.Cs_r).*gn.hc + gn.Cs_r.*gn.h(j,i);
   izzr(1:nn)=iz0r+izeta.*(1.0 + iz0r./gn.h(j,i));

   % add extra top level at 20 m above sea level and bottom level at maxdepth.
   extzzr=ones(nl+2)*NaN;
   extzzr(1)=20;                          % add top
   extzzr(2:nl+1)=zzr(1:nl);
   extzzr(nl+2)=-maxdepth;                % add bottom

   % vertical interpolation of variable at larger domain grid and
   % horizontal interpolation to the nested grid
   for kk=1:4
    jj=row_index(kk);
    ii=col_index(kk);
    itempdata=interp1(extzzr(1:nl+2),exttemp(1:nl+2,jj,ii),izzr(1:nn),'linear');
    isaltdata=interp1(extzzr(1:nl+2),extsalt(1:nl+2,jj,ii),izzr(1:nn),'linear');
    iudata=interp1(extzzr(1:nl+2),extu(1:nl+2,jj,ii),izzr(1:nn),'linear');
    ivdata=interp1(extzzr(1:nl+2),extv(1:nl+2,jj,ii),izzr(1:nn),'linear');
    roms.temp(1:nn,j,i)=roms.temp(1:nn,j,i)+(itempdata.*weight(kk))';
    roms.salt(1:nn,j,i)=roms.salt(1:nn,j,i)+(isaltdata.*weight(kk))';
    ntrue_u(1:nn,j,i)=ntrue_u(1:nn,j,i)+(iudata.*weight(kk))';
    ntrue_v(1:nn,j,i)=ntrue_v(1:nn,j,i)+(ivdata.*weight(kk))';
   end

   % post-processing for NaN values
   Ip=find(  finite(squeeze(roms.salt(:,j,i))) );
   Iq=find( ~finite(squeeze(roms.salt(:,j,i))) ); 
   
   if ( length(Ip) < 1)
    error(['error at j=' num2str(j) ' ,  i = ' num2str(i) ' no data!'])
   elseif( length(Iq) > 1)
    error(['fix NaN value at j=' num2str(j) ' ,  i = ' num2str(i)])
   end

   roms.zeta(j,i)=izeta; % updata zeta values

 end  % of if ( gn.mask_rho(j,i) > 0 ) ; sea - we works on a cell under water

end % of j 
            disp(['done at  i = ' num2str(i) ',    j = ' num2str(j)])
end % of i  grid point in the nested (smaller) domain


%(3)roate true eastward and northward velocities to roms u and v
%
%   unit of g.angle is radian ( 0.74 = 42 degree )
%
%  | true_u | = | cos(angle) -sin(angle) | | u |
%  | true_v |   | sin(angle)  cos(angle) | | v |
%  where (u,v)' is roms space vector
%
%  | u | = |  cos(angle) +sin(angle) | | true_u |
%  | v |   | -sin(angle)  cos(angle) | | true_v |


nubar= cos(gn.angle).*ntrue_ubar + sin(gn.angle).*ntrue_vbar;
nvbar=-sin(gn.angle).*ntrue_ubar + cos(gn.angle).*ntrue_vbar;

clear cos_angle3D sin_angle3D

cos_angle3D=repmat( cos(gn.angle) ,[1 1 nn ]);
cos_angle3D=permute(cos_angle3D, [3 1 2]);
sin_angle3D=repmat( sin(gn.angle) ,[1 1 nn ]);
sin_angle3D=permute(sin_angle3D, [3 1 2]);

nu= cos_angle3D.*ntrue_u + sin_angle3D.*ntrue_v;
nv=-sin_angle3D.*ntrue_u + cos_angle3D.*ntrue_v;

%(4)interpolate nubar and nvbar on rho-points to ubar and vbar on velocity points
% rho    points 82x42 
% ubar   points 82x41
% vbar   points 81x42

clear ubar vbar u v

[M N]=size(nubar);
roms.ubar = ( nubar(:,1:N-1)+nubar(:,2:N) ) * 0.5;
roms.vbar = ( nvbar(1:M-1,:)+nvbar(2:M,:) ) * 0.5;
roms.u = ( nu(:,:,1:N-1)+nu(:,:,2:N) ) * 0.5;
roms.v = ( nv(:,1:M-1,:)+nv(:,2:M,:) ) * 0.5;

maxtemp = max(max(max(roms.temp)));
mintemp = min(min(min(roms.temp(roms.temp>0))));
maxsalt = max(max(max(roms.salt)));
minsalt = min(min(min(roms.salt(roms.salt>0))));
disp(['================================================================'])
disp([' max. temp = ' num2str(maxtemp) ' min. temp = ' num2str(mintemp)])
disp([' max. salt = ' num2str(maxsalt) ' min. salt = ' num2str(minsalt)])

maxubar = max(max(roms.ubar));
minubar = min(min(roms.ubar));
maxvbar = max(max(roms.vbar));
minvbar = min(min(roms.vbar));
disp(['================================================================'])
disp([' max. ubar = ' num2str(maxubar) ' min. ubar = ' num2str(minubar)])
disp([' max. vbar = ' num2str(maxvbar) ' min. vbar = ' num2str(minvbar)])

maxu = max(max(max(roms.u)));
minu = min(min(min(roms.u)));
maxv = max(max(max(roms.v)));
minv = min(min(min(roms.v)));
disp(['================================================================'])
disp([' max. u = ' num2str(maxu) ' min. u = ' num2str(minu)])
disp([' max. v = ' num2str(maxv) ' min. v = ' num2str(minv)])

% check if any NaN salt value in data

     %repmat Replicate and tile an array
     %B = repmat(A,M,N) creates a large matrix B consisting 
     %of an M-by-N tiling of copies of A.

  indexNaN1 = find( ~finite(roms.salt) );
  rho_mask3D= repmat(gn.mask_rho,[1,1,nn]);
  rho_mask3D= permute(rho_mask3D,[3,1,2]);
  indexNaN2 = find( rho_mask3D(indexNaN1) == 1);
  if ( length(indexNaN2) > 1 )
     beep;
     disp(' there are NaN salt values in data ')
     beep;
     error(' error! ')  
  end

plotsalt=1; % Do you want to plot temp/salt at top and bottom level? yes=1 no=0
if (plotsalt)
  figure
  subplot(2,2,1)
  pcolorjw( gn.lon_rho, gn.lat_rho, squeeze(roms.temp(nn,:,:)) )
  title(' temp at top level ')
  caxis([mintemp maxtemp])
  colorbar
  subplot(2,2,2)
  pcolorjw( gn.lon_rho, gn.lat_rho, squeeze(roms.temp(1,:,:)) )
  title(' temp at bottom level ')
  caxis([mintemp maxtemp])
  colorbar
  subplot(2,2,3)
  pcolorjw( gn.lon_rho, gn.lat_rho, squeeze(roms.salt(nn,:,:)) )
  title(' salt at top level ')
  caxis([minsalt maxsalt])
  colorbar
  subplot(2,2,4)
  pcolorjw( gn.lon_rho, gn.lat_rho, squeeze(roms.salt(1,:,:)) )
  title(' salt at bottom level ')
  caxis([minsalt maxsalt])
  colorbar
end


% check if any NaN u value in data

     %repmat Replicate and tile an array
     %B = repmat(A,M,N) creates a large matrix B consisting 
     %of an M-by-N tiling of copies of A.

  indexNaN1 = find( ~finite(roms.u) );
  u_mask3D= repmat(gn.mask_u,[1,1,nn]);
  u_mask3D= permute(u_mask3D,[3,1,2]);
  indexNaN2 = find( u_mask3D(indexNaN1) == 1);
  if ( length(indexNaN2) > 1 )
     beep;
     disp(' there are NaN u values in data ')
     beep;
     error(' error! ')  
  end


% check if any NaN v value in data

     %repmat Replicate and tile an array
     %B = repmat(A,M,N) creates a large matrix B consisting 
     %of an M-by-N tiling of copies of A.

  indexNaN1 = find( ~finite(roms.v) );
  v_mask3D= repmat(gn.mask_v,[1,1,nn]);
  v_mask3D= permute(v_mask3D,[3,1,2]);
  indexNaN2 = find( v_mask3D(indexNaN1) == 1);
  if ( length(indexNaN2) > 1 )
     beep;
     disp(' there are NaN v values in data ')
     beep;
     error(' error! ')  
  end

  % check volume conservation across open boundaries
  %
  %  volume transport = sum ( ubar*h ) + sum( vbar*h)

  vt_north = roms.vbar(end,:).*(gn.h(end,:)+gn.h(end-1,:))*0.25.*(dx(end,:)+dx(end-1,:)); % + outflow
  vt_south = roms.vbar(1,:)  .*(gn.h(1,:)+gn.h(2,:)      )*0.25.*(dx(1,:)+dx(2,:));       % + inflow 
  vt_west  = roms.ubar(:,1)  .*(gn.h(:,1)+gn.h(:,2)      )*0.25.*(dy(:,1)+ dy(:,2));      % + inflow
  vt_east  = roms.ubar(:,end).*(gn.h(:,end)+gn.h(:,end-1))*0.25.*(dy(:,end)+dy(:,end-1)); % + outflow

  vt_total = -sum(vt_north) + sum(vt_south) ...
           +sum(vt_west)  - sum(vt_east)     ; % + volume increase in the domain

  disp([' transport from the nested domain (interpolated data) '])
  disp([' north = ' num2str(-sum(vt_north)) ' south = ' num2str(sum(vt_south)) ...
        ' west = ' num2str(sum(vt_west)) ' east = ' num2str(-sum(vt_east)) ])
  disp([' total volume increase through open boundaries = ' num2str(vt_total) ' m^3/s '])
  disp(['   '])


write_roms_init(out_file,roms) % append data to an existing zeros initial file
                               % and close it.
disp([' appending data to initial file done ..... '])
disp([' Finished ' which(mfilename) ])
disp('  ')


plotcompare=1;     % compare large and nested domains data
if (plotcompare)
  figure
  roms_zview(out_file,'temp',1,-10,gn,3,1,'k')
  ax=axis;
  colorbar
  cax=caxis;
  title('temp: nested domain interpolated data at 10 m')

  figure
  roms_zview(out_file,'salt',1,-50,gn,3,1.0,'k')
  ax=axis;
  colorbar
  cax=caxis;
  title('salt: nested domain interpolated data at 50 m')

  figure
  roms_zview(out_file,'temp',1,-250,gn,3,1,'k')
  ax=axis;
  colorbar
  cax=caxis;
  title('temp: nested domain interpolated data at 250 m')

  figure
  roms_zview(out_file,'temp',1,-400,gn,3,1,'k')
  ax=axis;
  colorbar
  cax=caxis;
  title('temp: nested domain interpolated data at 400 m')

  figure
  roms_zview(out_file,'salt',1,-400,gn,3,1,'k')
  ax=axis;
  colorbar
  cax=caxis;
  title('salt: nested domain interpolated data at 400 m')

end

return
