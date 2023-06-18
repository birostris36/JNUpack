function [volume,heat,salinity] = ij_transport_onefile(filename,g,selectdepth);
% ====================================================================
% [volume,heat,salinity] = ij_transport_onefile(filename,g,selectdepth);
% calculate transports of  volume(in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
%
% across a line segment (slice along a constant I or J)
%
% 'filename'       = history or average file name
% 'grid_nick_name' = grid information such as 'eas' 'hudson' 'latte'
% 'selectdepth' = surface to which depth (m) such as 100, 500, 50000
%
% keep in mind that
% vertical coordinate changes in time = h + zeta(t) in ROMS
% ====================================================================
% version 1 for single   input file  and single segment.
% version 2 for multiple input files and single segment.
% onefile   for single   input file  and single segment.

  if( nargin < 1)
    %filename='his_eas_005_0014.nc';
    filename='jes_clm_his_0149.nc';
    g = grd('jes');
    selectdepth=5000;
    %g = grd('eas');
    %g = grd('hudson');
    selectdepth  = input('surface to which depth (m) = ');
  end

% *************************************************************
%
%   END OF USER DEFINED VARIABLES
%
% *************************************************************

% size of grids
[r,c] = size ( g.lon_rho );

% transport from surface to which depth (m)

    if ( selectdepth > 0 ) 
     selectdepth = selectdepth*-1;
    end

%  graphical input or text input

qq = input('Do you want to use the keyboard or the mouse? ','s');

if isempty(qq)
  qq = 'm';
end

if strcmpi(qq(1),'m');
  disp(' ');
  f = input('Which Figure? (Return to create new figure) ');
  
  if isempty(f)
     f = figure;
     figure(f)
     set(f,'Position',[300 300 900 600]);
     pcolorjw(g.lon_rho,g.lat_rho,g.mask_rho_nan);
%    pnc;
  else;
     figure(f)
  end;
  
  for n=1:2;
    figure(f)
    [endpt_lon(n),endpt_lat(n)]=ginput(1);
  end;

  for cpts=1:2 % corner points
    dist = sqrt(  ( g.lon_rho - endpt_lon(cpts) ).*( g.lon_rho - endpt_lon(cpts) ) + ...
                  ( g.lat_rho - endpt_lat(cpts) ).*( g.lat_rho - endpt_lat(cpts) ) );
    ind=find( min( dist(:) ) == dist );
    % closest points row and column indice
    row_index = mod ( ind - 1, r ) + 1;
    col_index = floor( (ind - 1) / r ) + 1;
    endpt_lon(cpts)=col_index;
    endpt_lat(cpts)=row_index;
  end

else;

    %choose two grid point (i,j)
    disp(' ============================ ')
    disp(' Enter two grid points (i,j): i = longitude index, j = latitude index')
    disp(' along fixed i or j           ') 
    disp(['   1 < i < ',num2str(size(g.h,2))])
    disp(['   1 < j < ',num2str(size(g.h,1))])

    endpt_lon(1) = input('Starting i = ');
    endpt_lat(1) = input('Starting j = ');
    endpt_lon(2) = input('  Ending i = ');
    endpt_lat(2) = input('  Ending j = ');

end;

    if ( ( endpt_lon(1) ~= endpt_lon(2) ) &   ...
         ( endpt_lat(1) ~= endpt_lat(2) )    )
       disp(' either i or j should be same ')
       diflon = abs( endpt_lon(1) - endpt_lon(2) );
       diflat = abs( endpt_lat(1) - endpt_lat(2) );
       if( diflon > diflat )
         endpt_lat(1) = round( ( endpt_lat(1) + endpt_lat(2) )*0.5 ) ;
         endpt_lat(2) = endpt_lat(1);
         disp(' latitude (j) is adjusted ')
       else
         endpt_lon(1) = round( ( endpt_lon(1) + endpt_lon(2) )*0.5 ) ;
         endpt_lon(2) = endpt_lon(1);
         disp(' longitude (i) is adjusted ')
       end
    end
    if(endpt_lon(1) >  endpt_lon(2))
     ttt=endpt_lon(1);
     endpt_lon(1)=endpt_lon(2);
     endpt_lon(2)=ttt;
    end
    if(endpt_lat(1) >  endpt_lat(2))
     ttt=endpt_lat(1);
     endpt_lat(1)=endpt_lat(2);
     endpt_lat(2)=ttt;
    end

    xx=g.lon_rho(endpt_lat,endpt_lon);
    yy=g.lat_rho(endpt_lat,endpt_lon);


%   grid information

    figure(f)
    %pcolorjw(g.lon_rho,g.lat_rho,g.h.*g.mask_rho_nan) 
    hold on
    plot(xx,yy,'x-k','LineWidth',2)
    colorbar
    title('bottom topography (m)')

   %N  is the number of vertical levels
   %hz is thickness  of each level
    N = g.N;  
    [M L]=size(g.h);
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));

%   load model output file
 
    nc=netcdf(filename,'read');
    disp([' opening your data file: ', filename])
    zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
    u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
    v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
    temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
    salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
    close(nc)

%   vertical coordinate changes in time 
%   because sea surface height changes in time.
%   thickness of each layer changes propotional to total water thicknes.

    h_total = g.h + zeta;       %total water thickness
    for level=1:N               %thickness of each layer
     Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
    end

% *******************************************************************
% FOR FIXED I (or longitude)
% *******************************************************************

  if  ( endpt_lon(1) == endpt_lon(2) )  

    % along a longitude
    %    _ __ 
    %     |
    %     | u  Hz(k+1,j,ii) 
    %    _|__
    %     |
    %     | u  Hz(k  ,j,ii)
    %    _|__

    ii=endpt_lon(1);    
    j1=endpt_lat(1);
    j2=endpt_lat(2);

    % average Hz to  Arakawa-C u points 

    Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
     z_u(N,:,:)=-Hz_u(N,:,:);             % z @ bottom of each layer
    for k=N-1:-1:1
     z_u(k,:,:)=z_u(k+1,:,:)-Hz_u(k,:,:);
    end

    % actual calculation of transport

    sum_segment=0;
    for j=j1:j2
      Ip=find( z_u(:,j,ii) > selectdepth );
      lastlevel=min([N; min(Ip)]);
      sum=0;
      for k=N:-1:lastlevel         % whole cell transport
          sum=sum+u(k,j,ii)*Hz_u(k,j,ii);
      end
      if (lastlevel == N)          % partial cell transport 
           partialH=zeta(j,ii)-selectdepth;
           sum=sum+u(N,j,ii)*partialH;
      elseif (lastlevel ~= 1)     
          partialH=z_u(k,j,ii)-selectdepth;
          sum=sum+u(k-1,j,ii)*partialH;
      end
      sum_segment=sum_segment+sum*dy_u(j,ii);
    end

% *******************************************************************
% FOR FIXED J (or latitude)
% *******************************************************************

  elseif ( endpt_lat(1) == endpt_lat(2) )   

    % along a latitude
    %    _ __ 
    %     |
    %     | v  Hz(k+1,jj,i) 
    %    _|__
    %     |
    %     | v  Hz(k  ,jj,i)
    %    _|__
   
    jj=endpt_lat(1);
    i1=endpt_lon(1); 
    i2=endpt_lon(2); 

    % average Hz to  Arakawa-C v points 

    Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
     z_v(N,:,:)=-Hz_v(N,:,:);             % z @ bottom of each layer
    for k=N-1:-1:1
     z_v(k,:,:)=z_v(k+1,:,:)-Hz_v(k,:,:);
    end

    % actual calculation of transport

    sum_segment=0;
    for i=i1:i2
      Ip=find( z_v(:,jj,i) > selectdepth );
      lastlevel=min([N; min(Ip)]);
      sum=0;
      for k=N:-1:lastlevel         % whole cell transport
          sum=sum+v(k,jj,i)*Hz_v(k,jj,i);
      end
      if (lastlevel == N)          % partial cell transport 
          partialH=zeta(jj,i)-selectdepth;
          sum=sum+v(N,jj,i)*partialH;
      elseif (lastlevel ~= 1)         
          partialH=z_v(k,jj,i)-selectdepth;
          sum=sum+v(k-1,jj,i)*partialH;
      end
      sum_segment=sum_segment+sum*dx_v(jj,i);
    end

  else
  
      error(' either i or j must be same ! ')
   
  end

  disp([' volume transport = ', num2str(sum_segment) ,' m^3/s '])
  disp(['                  = ', num2str(sum_segment/1e+6) ,' x 10^6m^3/s '])
