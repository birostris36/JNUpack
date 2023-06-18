function [volume,heat,salinity] = ij_transport(filename,grid_nick_name);
% ====================================================================
% calculate transports of  volume(in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
%
% across a line segment (slice along a constant I or J)
%
% 'filename'       = history or average file name
% 'grid_nick_name' = grid information such as 'eas' 'hudson' 'latte'
%
% keep in mind that
% vertical coordinate changes in time = h + zeta(t) in ROMS
% ====================================================================


%   choose two grid point (i,j)

    disp(' ============================ ')
    disp(' Enter two grid points (i,j): ')
    disp(' along fixed i or j           ') 

    endpt_lon(1) = input('Starting i = ');
    endpt_lat(1) = input('Starting j = ');
    endpt_lon(2) = input('  Ending i = ');
    endpt_lat(2) = input('  Ending j = ');

    selectdepth  = input('surface to which depth (m) = ');
    if ( selectdepth > 0 ) 
     selectdepth = selectdepth*-1;
    end

%   grid information

    g = grd(grid_nick_name);
   %g = grd('hudson'); g = grd('latte'); g = grd('eas');

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
%     temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
%     salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
    close(nc)

%   vertical coordinate changes in time 
%   because sea surface height changes in time.
%   thickness of each layer changes propotional to total water thicknes.

    h_total = g.h +squeeze(zeta(1,:,:));       %total water thickness
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
      lastlevel=min(Ip);
      sum=0;
      for k=N:-1:lastlevel         % whole cell transport
          sum=sum+u(k,j,ii)*Hz_u(k,j,ii);
      end
      partialH=z_u(k,j,ii)-selectdepth;
      sum=sum+u(k+1,j,ii)*partialH;% partial cell transport 
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
      lastlevel=min(Ip);
      sum=0;
      for k=N:-1:lastlevel         % whole cell transport
          sum=sum+v(k,jj,i)*Hz_v(k,jj,i);
      end
      partialH=z_v(k,jj,i)-selectdepth;
      sum=sum+u(k+1,jj,i)*partialH;% partial cell transport 
      sum_segment=sum_segment+sum*dx_v(jj,i);
    end

  else
  
      error(' either i or j must be same ! ')
   
  end

  disp([' volume transport = ', num2str(sum_segment) ,' m^3/s '])
