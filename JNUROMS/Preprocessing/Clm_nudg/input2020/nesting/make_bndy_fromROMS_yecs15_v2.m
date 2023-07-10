function make_bndy_fromROMS_yecs15_v2
% ==========================================================================
% build Boundary fields for a nested roms grid
% from data in a larger roms grid.
% Data file could be history or average file.
%
% Algorithm: 
% Get vertical distribution of a variable for a given grid cell,
% (j,i) nested grid cell.
% (1) find nearest four points in larger domain
% (2) calculate weight using distances
% (3) vertical interpolation at four points
% (4) horizontal interpolation base on wights in (2)
% BJ Choi, sept-15-2005
%
% in_file      = larger domain roms data file (his.nc or avg.nc)
% analysis_day = time of in_file (in days)
% out_file     = boundary file for the nested domain (init.nc)
%
% USE: nc_create_roms_bndy.m  ( create empty file )
%      write_roms_bndy.m      ( append data )
%      grd.m
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
% ==========================================================================
% clear all
clear
close all

disp(' read grid information of larger domain and the nested smaller domain')

if exist('gl')~=1
   gl = grd('yecs15');         % (1) Get Larger grid information
end

if exist('gn')~=1
    gn = grd('jeju2'); % (2) Get Nested grid information
end

% (3) roms data file from a larger domain
% please enter one of data file with full path of directory
%in_file='/data4/roms/ecs10/ecs10_avg-2006-1-n_6400.nc';
% in_file='G:\case\NWP4\year_10\avg\avg_NWP4_10year_0001.nc';
in_file='/data/JHMoon/model/ROMS839/project/yecs15_auto/output/yecs15_avg_0091.nc';
mm=length(in_file);
%    files: (year 2000)
%      ocean_avg_0002.nc  ...  ocean_avg_0367.nc
%      (day 1)                 (day 366)


% (4) output file, boundary file for the nested grid
%out_file = 'ctz_lg_bndy2000June.nc';
%  out_file = 'D:\MODEL\DATA\ecsy12\ecsy12_nest_bndy_365_v2.nc';
 out_file = './JEJU2_bry_by_yecs15_M04.nc';

% (5) start date and end date
% May  31, 2000 = 152
% June 30, 2000 = 182
% July 31, 2000 = 213
% Aug  31, 2000 = 244
% Sep   1, 2000 = 245
% startday=152;  % May 31, 2000 
% endday  =183;  % July 1, 2000 
%  startday=182; %195;  % July  15, 2006
%  endday  =238;  % Aug   27, 2006
startday=90   ;  % Mar  31 2017
endday  =120 ;  % Aug  1 2007
% Information for the boundary data file
% and create an empty boundary file
disp([' create an empty boundary file = ' out_file])
noclobber = 0;

base_date = [1993 1 1 0 0 0];
time_variable = 'ocean_time';
time_variable_units = basedate2str(base_date);
roms.grd = gn;
grd_file = gn.grd_file;
sourcestr = [ 'From NWP4 : 1/12 esolution, daily '];
details =   [ 'Quick boundary file from NWP4 daily mean ' ...
              'by script ' which(mfilename) ];
donuts = 0;

nc_create_roms_bndy_by_yecs  % create an nc file, write variables and close it.
disp([' creating an empty bndy file done ..... '])
disp('  ')

% ================= End of Input Parameters ==========================


% size of grids
[r,c] = size ( gl.lon_rho );
   nl = gl.N;  
[M N] = size ( gn.lon_rho );
   nn = gn.N;

% find deepest depth
maxdepth=max([max(max(gl.h)) max(max(gn.h))])+500;

% define grid spacing dx and dy
dx = 1./gn.pm;
dy = 1./gn.pn;

for indexday=startday:endday

 day=num2str(indexday);
 lday=length(day);
% in_file(mm-3-lday:mm-3)=day+5;
 in_file(mm-3-lday+1:mm-3)=day;
 s=dir(in_file);  % s.bytes ==> file size

% if (s.bytes > 10000000 )   
%   analysis_day=indexday-0.5; % middle of day
ocn_day = datenum([2017,01,01]) - datenum(base_date);
  analysis_day=ocn_day+indexday-0.5; % middle of day

  % read data from larger domain
  disp([' read data from the larger domain roms file = ' in_file])
  disp([' wait ..... '])
  nc=netcdf(in_file,'read');
  temp=nc{'temp'}(:);
  salt=nc{'salt'}(:);
  zeta=nc{'zeta'}(:);
  ubar_dump=nc{'ubar'}(:);
  vbar_dump=nc{'vbar'}(:);
  u_dump=nc{'u'}(:);
  v_dump=nc{'v'}(:);
  close(nc)
  disp([' reading done ..... '])
  disp('  ')
%--------------------------------------------------------------------------
% for above ROMS-3.2, you have to control the FillValue !!!!
% Added By Peter CNU, 17th of Feb. 2009
t_index=find(temp>500); temp(t_index)=0;
t_index=find(salt>500); salt(t_index)=0;
t_index=find(zeta>500); zeta(t_index)=0;
t_index=find(ubar_dump>500); ubar_dump(t_index)=0;
t_index=find(vbar_dump>500); vbar_dump(t_index)=0;
t_index=find(u_dump>500); u_dump(t_index)=0;
t_index=find(v_dump>500); v_dump(t_index)=0;
clear t_index 
%--------------------------------------------------------------------------
  %(1)interpolate ubar_dump and vbar_dump on rho-points
  % rho    points 514x510 
  % ubar   points 514x509
  % vbar   points 513x510

  [M N]=size(ubar_dump);
  ubar =ones(M,N+1)*NaN;
  ubar(:,2:N) = ( ubar_dump(:,1:N-1)+ubar_dump(:,2:N) ) * 0.5;
  ubar(:,1)   = ubar_dump(:,1);
  ubar(:,N+1) = ubar_dump(:,N);

  [M N]=size(vbar_dump);
  vbar =ones(M+1,N)*NaN;
  vbar(2:M,:) = ( vbar_dump(1:M-1,:)+vbar_dump(2:M,:) ) * 0.5;
  vbar(1,:)   = vbar_dump(1,:);
  vbar(M+1,:) = vbar_dump(M,:);
  

  [L M N]=size(u_dump);
  u =ones(L,M,N+1)*NaN;
  u(:,:,2:N) = ( u_dump(:,:,1:N-1)+u_dump(:,:,2:N) ) * 0.5;
  u(:,:,1)   = u_dump(:,:,1);
  u(:,:,N+1) = u_dump(:,:,N);

  [L M N]=size(v_dump);
  v =ones(L,M+1,N)*NaN;
  v(:,2:M,:) = ( v_dump(:,1:M-1,:)+v_dump(:,2:M,:) ) * 0.5;
  v(:,1,:)   = v_dump(:,1,:);
  v(:,M+1,:) = v_dump(:,M,:);

  %(2)roate them to get true eastward and northward velocities.
  %   unit of g.angle is radian ( 0.74 = 42 degree )
  %
  %  | true_u | = | cos(angle) -sin(angle) | | u |
  %  | true_v |   | sin(angle)  cos(angle) | | v |
  %  where (u,v)' is roms space vector

  true_ubar=cos(gl.angle).*ubar - sin(gl.angle).*vbar;
  true_vbar=sin(gl.angle).*ubar + cos(gl.angle).*vbar;

  cos_angle3D=repmat( cos(gl.angle) ,[1 1 L ]);
  cos_angle3D=permute(cos_angle3D, [3 1 2]);
  sin_angle3D=repmat( sin(gl.angle) ,[1 1 L ]);
  sin_angle3D=permute(sin_angle3D, [3 1 2]);

  true_u=cos_angle3D.*u - sin_angle3D.*v;
  true_v=sin_angle3D.*u + cos_angle3D.*v;


  % size of grids
  [r,c] = size ( gl.lon_rho );
     nl = gl.N;  
  [M N] = size ( gn.lon_rho );
     nn = gn.N;


  % extraplolate data (variables) horizonally
  % into the land. 
  % three land grid points next to the coast will have data.
  disp([' horizontal extrapolation of original data '])

  mask=gl.mask_rho;
  mask_temp=mask;
  for numofextrapol=1:3
   Iland=find(  mask_temp == 0 );
   num_land_grid=length(Iland);
   for i=1:num_land_grid
   
     ind = Iland(i);
     row_index = mod ( ind - 1, r ) + 1;
     col_index = floor( (ind - 1) / r ) + 1;
     extflag=0;

     if(     (col_index > 2) && (mask(row_index,col_index-1) == 1) )
       oj=row_index; oi=col_index-1; extflag=1;
     elseif( (col_index < c) && (mask(row_index,col_index+1) == 1) )
       oj=row_index; oi=col_index+1; extflag=1;
     elseif( (row_index > 2) && (mask(row_index-1,col_index) == 1) )
       oj=row_index-1; oi=col_index; extflag=1;
     elseif( (row_index < r) && (mask(row_index+1,col_index) == 1) )
       oj=row_index+1; oi=col_index; extflag=1;
     end

     if( extflag )
       % 2D variables
       zeta(row_index,col_index)       = zeta(oj,oi);
       true_ubar(row_index,col_index)  = true_ubar(oj,oi);
       true_vbar(row_index,col_index)  = true_vbar(oj,oi);
       % 3D variables
       temp(:,row_index,col_index)   = temp(:,oj,oi);
       salt(:,row_index,col_index)   = salt(:,oj,oi);
       true_u(:,row_index,col_index) = true_u(:,oj,oi);
       true_v(:,row_index,col_index) = true_v(:,oj,oi);
       % reset mask value
       mask_temp(row_index,col_index)=1; 
     end

   end % of for i=1:num_land_grid
   mask=mask_temp; 
  end  % of for numofextrapol=1:2


  % find land mask from gl grid
  ocean=ones(r,c);
  land =ones(r,c)*1.e20;
  Isea=find( mask > 0);
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

  extu(1,:,:)=true_u(1,:,:);         % add bottom ( at -maxdepth)
  extu(2:nl+1,:,:)=true_u(1:nl,:,:); % data
  extu(nl+2,:,:)=true_u(nl,:,:);     % add top    ( 20 m above sea level)

  extv(1,:,:)=true_v(1,:,:);         % add bottom ( at -maxdepth)
  extv(2:nl+1,:,:)=true_v(1:nl,:,:); % data
  extv(nl+2,:,:)=true_v(nl,:,:);     % add top    ( 20 m above sea level)

  % initailize roms structure (data file)
  disp([' initailize roms structure  '])

  roms.time = analysis_day;
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

  % ====================================================
  % vertical and horizonatal interpolation of variables 
  % such as temp, salt, zeta, u, v, ubar and vbar
  % from larger domain to the nested domain.
  % ====================================================

  disp([' ================================================ '])
  disp([' interpolating zeta, temp, salt, u, v, ubar and vbar on the nested domain '])
  disp([' wait .....                              '])
  disp(['                                         '])

  % west  and east
  % north and south
  for i=1:N 

    if( i==1 || i==2 || i==(N-1) || i==N )
     jindex=[1:1:M];                         % for west  and east
    else 
     jindex=[1 2 (M-1) M];                   % for north and south
    end 
  
    for j=jindex

    % disp([' index (j,i) ' num2str(j) '  ' num2str(i) ])

    % =================================================== %

     if ( gn.mask_rho(j,i) > 0 ) % sea; we works on a cell under water

       % find 4 nearest points
       % ind --> index of the 4 points
       % Assume the projection is ok to do this.
       d = sqrt ( (gl.lon_rho - gn.lon_rho(j,i)).^2 + ...
                  (gl.lat_rho - gn.lat_rho(j,i)).^2 ); 
       d = d.*land; % set d = 1.e20 for land points
       d_temp = d;
       ind=[];
       while ( length( ind ) < 4  )
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
        xx=gl.lon_rho(jj,ii);
        yy=gl.lat_rho(jj,ii);
        dis(kk)=m_lldist([xx0  xx],[yy0 yy]);
       end
       sum_dis=sum( dis(1:4) );
       weight(1:4) = dis(1:4)./sum_dis;

       % transformation from s-coordinate to z-coordinate
       %z0r = (grid.sc_r-grid.Cs_r).*grid.hc + grid.Cs_r.*grid.h;
       %zzr = z0r + squeeze(zeta).*(1.0 + z0r./grid.h);

       % vertical interpolation and extrapolation 
       % interpolate zeta, ntrue_ubar and ntrue_vbar horizontally
       izeta=0;              % interpolated sea surface, zeta
       ihl=0;                % interpolated depth of water, h
       for kk=1:4  
        jj=row_index(kk);
        ii=col_index(kk);
%%    for classical Song and Haidvogel version 
%%    This is for parent model
        z0r=(gl.sc_r-gl.Cs_r).*gl.hc + gl.Cs_r.*gl.h(jj,ii);
        zzr(1:nl,kk)=z0r+zeta(jj,ii).*(1.0 + z0r./gl.h(jj,ii));
%%    for the new vertical transformation from Shechepetkin 2005        
%%    This is for parent model
%        z0r=(gl.hc.*gl.sc_r+gl.Cs_r.*gl.h(j,i))./(gl.h(j,i) + gl.hc);
%        zzr(1:nl,kk)=zeta+(zeta+gl.h(j,i)).* z0r;
       
        izeta=izeta+zeta(jj,ii).*weight(kk);
        ihl  =ihl  +gl.h(jj,ii).*weight(kk);
        ntrue_ubar(j,i)=ntrue_ubar(j,i)+true_ubar(jj,ii).*weight(kk);
        ntrue_vbar(j,i)=ntrue_vbar(j,i)+true_vbar(jj,ii).*weight(kk);
       end

       % apply volume flux conservation across open boundary
       % vfactor = 1 if ihl = gn.h(j,i)
       vfactor=(ihl+izeta)/(gn.h(j,i)+izeta);    
       ntrue_ubar(j,i)=ntrue_ubar(j,i) * vfactor ;
       ntrue_vbar(j,i)=ntrue_vbar(j,i) * vfactor ;
%%    for classical Song and Haidvogel version
       % for using the same level as parents
%        z0r=(gn.sc_r-gn.Cs_r).*gn.hc + gn.Cs_r.*gn.h(j,i);
%        izzr(1:nn)=z0r+izeta.*(1.0 + z0r./gn.h(j,i));  
       % for different level
%        iz0r=(gn.sc_r-gn.Cs_r).*gn.hc + gn.Cs_r.*gn.h(j,i);
%        izzr(1:nn)=iz0r+izeta.*(1.0 + iz0r./gn.h(j,i)); 

%%    for the new vertical transformation from Shechepetkin 2005
       % for using the same level as parents
%        z0r=(gn.hc.*gn.sc_r+gn.Cs_r.*gn.h(j,i))./(gn.h(j,i) + gn.hc);
%        izzr(1:nn)=izeta+(izeta+gn.h(j,i)).* z0r;
       % for different level as parent
       iz0r=(gn.hc.*gn.sc_r+gn.Cs_r.*gn.h(j,i))./(gn.h(j,i) + gn.hc);
       izzr(1:nn)=izeta+(izeta+gn.h(j,i)).* iz0r;
       
       % add extra top level at 20 m above sea level and bottom level at maxdepth.
       extzzr=ones(nl+2,4)*NaN;
       extzzr(1,1:4)=-maxdepth;                          % add bottom
       extzzr(2:nl+1,1:4)=zzr(1:nl,1:4);
       extzzr(nl+2,1:4)=20;                              % add top

       % horizontal interpolation
       for kk=1:4
        jj=row_index(kk);
        ii=col_index(kk);
        itempdata=interp1(extzzr(1:nl+2,kk),exttemp(1:nl+2,jj,ii),izzr(1:nn),'linear');
        isaltdata=interp1(extzzr(1:nl+2,kk),extsalt(1:nl+2,jj,ii),izzr(1:nn),'linear');
        iudata=interp1(extzzr(1:nl+2,kk),extu(1:nl+2,jj,ii),izzr(1:nn),'linear');
        ivdata=interp1(extzzr(1:nl+2,kk),extv(1:nl+2,jj,ii),izzr(1:nn),'linear');
        roms.temp(1:nn,j,i)=roms.temp(1:nn,j,i)+(itempdata.*weight(kk))';
        roms.salt(1:nn,j,i)=roms.salt(1:nn,j,i)+(isaltdata.*weight(kk))';
        ntrue_u(1:nn,j,i)=ntrue_u(1:nn,j,i)+(iudata.*weight(kk))';
        ntrue_v(1:nn,j,i)=ntrue_v(1:nn,j,i)+(ivdata.*weight(kk))';
       end

       % post-processing for NaN values
       Ip=find(  isfinite(squeeze(roms.salt(:,j,i))) );
       Iq=find( ~isfinite(squeeze(roms.salt(:,j,i))) ); 
   
       if ( length(Ip) < 1)
        error(['error at j=' num2str(j) ' ,  i = ' num2str(i) ' no data!'])
       elseif( length(Iq) > 1)
        error(['fix NaN value at j=' num2str(j) ' ,  i = ' num2str(i)])
       end

       roms.zeta(j,i)=izeta; % updata zeta values

     end % of if-statement ( gn.mask_rho(j,i) > 0 )

     % =================================================== %

    end % of for j=jindex

  end % of i=1:N 


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

  % apply mask

  roms.ubar = roms.ubar .* gn.mask_u;
  roms.vbar = roms.vbar .* gn.mask_v;


  write_roms_bndy(out_file,roms) % append boundary values to an existing initial file
                                 % and close it.
  disp([' appending boundary values to bndy file done ..... '])
  disp('  ')


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

  plotvt=0; % do you want to plot?  yes=1 no=0
  if (plotvt)

   figure(21)                               % plot ubar 
   subplot(2,1,1)
   plot( gn.lat_u(:,1), roms.ubar(:,1) )
   title(' Western boundary ubar ')
   subplot(2,1,2)
   plot( gn.lat_u(:,end), roms.ubar(:,end) )
   title(' Eastern boundary ubar ')

   figure(22)                               % plot vbar 
   subplot(2,1,1)
   plot( gn.lon_v(end,:), roms.vbar(end,:) )
   title(' Northern boundary vbar ')
   subplot(2,1,2)
   plot( gn.lon_v(1,:), roms.vbar(1,:) )
   title(' Southern boundary vbar ')

   figure(23)
   subplot(2,1,1)
   pcolorjw( gn.lon_v(end,:), [-gn.N:1:-1], squeeze( roms.v(:,end,:) ) )
   title(' Northern boundary v ')
   colorbar
   subplot(2,1,2)
   pcolorjw( gn.lon_v(1,:), [-gn.N:1:-1], squeeze( roms.v(:,1,:) ) )
   title(' Southern boundary v ')
   colorbar

   figure(31)       % plot vt along western and eastern boundary?
   subplot(2,1,1)
   plot( gn.lat_u(:,1), vt_west )
   title(' Western boundary volume transport ')
   subplot(2,1,2)
   plot( gn.lat_u(:,end), vt_east )
   title(' Eastern boundary  volume transport ')

   figure(32)       % plot vt along northern and southern boundary?
   subplot(2,1,1)
   plot( gn.lon_v(end,:), vt_north )
   title(' Northern boundary  volume transport ')
   subplot(2,1,2)
   plot( gn.lon_v(1,:), vt_south )
   title(' Southern boundary  volume transport  ')
  end

  plotcompare=0;  % compare volume transports between original data and interpolated data
  if (plotcompare)
  
    % northern boundary
    j=size(gn.lat_v,1); i=1;
    ddd = sqrt ( (gl.lon_rho - gn.lon_rho(j,i)).^2 + (gl.lat_rho - gn.lat_rho(j,i)).^2 ); 
    ind = find( ddd == min(min(ddd)) );
    row_index1 = mod ( ind - 1, r ) + 1;
    col_index1 = floor( (ind - 1) / r ) + 1;

    j=size(gn.lat_v,1); i=size(gn.lat_v,2);
    ddd = sqrt ( (gl.lon_rho - gn.lon_v(j,i)).^2 + (gl.lat_rho - gn.lat_v(j,i)).^2 ); 
    ind = find( ddd == min(min(ddd)) );
    row_index2 = mod ( ind - 1, r ) + 1;
    col_index2 = floor( (ind - 1) / r ) + 1;

    row=row_index2;
    co1=col_index1;
    co2=col_index2;
    vt_north_lg = vbar_dump(row,co1:co2) .*( gl.h(row,co1:co2)+gl.h(row-1,co1:co2) ) ...
                         .* ( 1./gl.pm(row,co1:co2)+ 1./gl.pm(row-1,co1:co2) )*0.25;  % + outflow

    % southern boundary
    j=1; i=1;
    ddd = sqrt ( (gl.lon_rho - gn.lon_v(j,i)).^2 + (gl.lat_rho - gn.lat_v(j,i)).^2 ); 
    ind = find( ddd == min(min(ddd)) );
    row_index3 = mod ( ind - 1, r ) + 1;
    col_index3 = floor( (ind - 1) / r ) + 1;

    j=1; i=size(gn.lat_v,2);
    ddd = sqrt ( (gl.lon_rho - gn.lon_v(j,i)).^2 + (gl.lat_rho - gn.lat_v(j,i)).^2 ); 
    ind = find( ddd == min(min(ddd)) );
    row_index4 = mod ( ind - 1, r ) + 1;
    col_index4 = floor( (ind - 1) / r ) + 1;

    row=row_index3;
    co1=col_index3;
    co2=col_index4;
    vt_south_lg = vbar_dump(row,co1:co2) .*( gl.h(row,co1:co2)+gl.h(row-1,co1:co2) ) ...
                    .* ( 1./gl.pm(row,co1:co2)+ 1./gl.pm(row-1,co1:co2) )*0.25;  % + inflow

    % eastern boundary
    ro1=row_index3;
    ro2=row_index1;
    col=col_index1;
    vt_west_lg = ubar_dump(ro1:ro2,col) .*( gl.h(ro1:ro2,col)+gl.h(ro1:ro2,col+1) ) ...
                    .* ( 1./gl.pn(ro1:ro2,col)+ 1./gl.pn(ro1:ro2,col+1) )*0.25;  % + inflow

    figure(33)
    subplot(2,1,1)
    plot( gl.lon_rho(row_index2,col_index1:col_index2), vt_north_lg )
    title(' Northern boundary volume tranport from large model ')
    subplot(2,1,2)
    plot( gl.lon_rho(row_index3,col_index3:col_index4), vt_south_lg )
    title(' Southern boundary volume tranport from large model ')

    disp(' transport from larger domain (original data) ')
    disp([' north = ' num2str(-sum(vt_north_lg)) ' south = ' num2str(sum(vt_south_lg)) ])
    disp([' west  = ' num2str(sum(vt_west_lg)) ])

  end % of if-statement (plotcompare)


% end % of if-statement (s.bytes > 10000000 ) 

end % of for-statement; indexday

disp([' Finished ' which(mfilename) ])
disp('  ')

return
