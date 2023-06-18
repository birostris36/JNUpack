% function [total_volume,total_salt,total_fresh] = xy_transport_onefile(filename,g,selectdepth);
% ====================================================================
% [total_volume,total_salt,total_fresh]= xy_transport_onefile(filename,g,selectdepth);
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
%
% USE: xy_transport_function.m
%
% ====================================================================
% onefile   for single   input file  and single segment.
% BJ Choi, Marhch07, 2006.

clear; close all

D_total_vol=[];
D_total_salt=[];
D_total_fresh=[];
g = grd('ecsy12');
selectdepth=70;
%     selectdepth  = input('surface to which depth (m) = ');

%     s_max = 32.76; % maximum salinity for the calculation of
%                    % freshwater flux - Hudson River, NY Bight.

% s_max = 34.00; % maximum salinity for the calculation of
s_max = 34.50; % maximum salinity for the calculation of
% freshwater flux - Hudson River, NY Bight.

% *************************************************************
%
%   END OF USER DEFINED VARIABLES
%
% *************************************************************

% size of grids
[r,c] = size ( g.lon_rho );
mask3d_rho=repmat(g.mask_rho,[1 1 g.N]);
mask3d_rho=permute(mask3d_rho,[3 1 2]);

% transport from surface to which depth (m)

if ( selectdepth > 0 )
    selectdepth = selectdepth*-1;
end

%  graphical input or text input
%
% qq = input('Do you want to use the keyboard or the mouse? ','s');
%
% if isempty(qq)
%     qq = 'm';
% end
%
% if strcmpi(qq(1),'m');
%     disp(' ');
%     f = input('Which Figure? (Return to create new figure) ');
%
%     if isempty(f)
%         f = figure;
%         figure(f)
%         set(f,'Position',[300 300 900 600]);
%         pcolorjw(g.lon_rho,g.lat_rho,g.h.*g.mask_rho_nan);
%     else;
%         figure(f)
%     end;
%
%     for n=1:2;
%         figure(f)
%         [endpt_lon(n),endpt_lat(n)]=ginput(1);
%     end;
%
% else;
%
%     disp(' ');
%     f = input('Which Figure? (Return to create new figure) ');
%
%     if isempty(f)
%         f = figure;
%         figure(f)
%         set(f,'Position',[300 300 900 600]);
%         pcolorjw(g.lon_rho,g.lat_rho,g.h.*g.mask_rho_nan);
%     else;
%         figure(f)
%     end;
%
%     %choose two grid points (lon, lat)
%     disp(' ============================ ')
%     disp(' Enter two grid points (lon, lat) ')
%     disp(' Example: -74.20  40.10  -73.50   40.10   for New Jesery line')
%     disp('   or     -73.50  40.70  -73.50   40.10   for Long Island line')
%     disp('   or      -73.6953  40.6177 -74.0103  40.3981 for Raritan Bay line ')
%     disp('   or      -73.8809  40.9249  -73.8611  40.9086 for Hudson source point ')
%     % disp('   or     -73.6643  40.6091 -73.9906  40.3817 for Raritan Bay line ')
%
%     endpt_lon(1) = input('Starting lon = ');
%     endpt_lat(1) = input('Starting lat = ');
%     endpt_lon(2) = input('  Ending lon = ');
%     endpt_lat(2) = input('  Ending lat = ');
%
% end;

% endpt_lon(1) = 122.2;
% endpt_lat(1) = 32;
% endpt_lon(2) = 126;
% endpt_lat(2) = 34.3;

% endpt_lon(1) = 123;
% endpt_lat(1) = 35;
% endpt_lon(2) = 125;
% endpt_lat(2) = 32;

endpt_lon(1) = 123;
endpt_lat(1) = 35;
endpt_lon(2) = 125;
endpt_lat(2) = 33.5;

disp(' Enter two grid points (lon, lat) ')
disp([ num2str(endpt_lon(1)) ' to ' num2str(endpt_lon(2))])
disp([ num2str(endpt_lat(1)) ' to ' num2str(endpt_lat(2))])

for cpts=1:2 % corner points
    dist = sqrt(  ( g.lon_rho - endpt_lon(cpts) ).*( g.lon_rho - endpt_lon(cpts) ) + ...
        ( g.lat_rho - endpt_lat(cpts) ).*( g.lat_rho - endpt_lat(cpts) ) );
    ind=find( min( dist(:) ) == dist );
    % closest points row and column indice
    row_index = mod ( ind - 1, r ) + 1;
    col_index = floor( (ind - 1) / r ) + 1;
    corner_endpt_col(cpts)=col_index;
    corner_endpt_row(cpts)=row_index;
end


% my xy_transport_onefile works only if corner_endpt_row(2) < corner_endpt_row(1).
if( corner_endpt_row(2) > corner_endpt_row(1)  )
    tmp_col=corner_endpt_col(2);
    tmp_row=corner_endpt_row(2);
    corner_endpt_col(2)=corner_endpt_col(1);
    corner_endpt_row(2)=corner_endpt_row(1);
    corner_endpt_col(1)=tmp_col;
    corner_endpt_row(1)=tmp_row;
    beep
    disp(' === switching two end points === ')
end

% longitude and latitude coordinate.

for i=1:length(endpt_lat)
    xx(i)=g.lon_rho(corner_endpt_row(i),corner_endpt_col(i));
    yy(i)=g.lat_rho(corner_endpt_row(i),corner_endpt_col(i));
end
distance_r = m_lldist ( xx, yy );

% figure(f)
% hold on
% plot(xx,yy,'x-k','LineWidth',2)
% colorbar
% title('bottom topography (m)')

%  transect information

% delj = j increasment
if( corner_endpt_col(2) >= corner_endpt_col(1) )
    delj=1;
    west2east_transect=1; % previously zonaltransect
else
    delj=-1;
    west2east_transect=0; % previously meridionaltransect
end

% deli = i increasment
if( corner_endpt_row(2) > corner_endpt_row(1) )
    deli=1;
else
    deli=-1;
end


%               i j
%             row col
% g.lon_rho: [142x254 double] for latte grid
% g.lon_u:   [142x253 double]
% g.lon_v:   [141x254 double]

xzero=g.lon_rho( corner_endpt_row(1), corner_endpt_col(1) );
yzero=g.lat_rho( corner_endpt_row(1), corner_endpt_col(1) );
xone=g.lon_rho( corner_endpt_row(2), corner_endpt_col(2) );
yone=g.lat_rho( corner_endpt_row(2), corner_endpt_col(2) );
slope=( yone-yzero) / (xone - xzero);
% A x + B y + C = 0;
A=slope;
B=-1;
C=-slope*xzero+yzero;
D=sqrt( A*A + B*B );
% distance = abs( A x + B y + C ) / D

%   grid information

%N  is the number of vertical levels
%hz is thickness  of each level
N = g.N;
[M L]=size(g.h);
hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
dx = 1./g.pm;
dy = 1./g.pn;
dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);




f_path=['/disk3/ecsy12_2009_2012_wt1/2012/'];

for day=1:365
    filename=[f_path 'avg_ecsy12_td_2012_kkl_v2_hycom_a4_wt1_' num2str(day,'%4.4i') '.nc'];

    %   load model output file

    nc=netcdf(filename,'read');
    disp([' opening your data file: ', filename])
    zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
    u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
    v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
    temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
    salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
    close(nc)

    salt = salt.*mask3d_rho; % zero for land, ** very important **

    szero = ones(size(salt)).*s_max;
    fresh = ( szero - salt ) ./ szero; % freshwater fraction
    fresh = fresh.*mask3d_rho; % zero for land, ** very important **

    %   vertical coordinate changes in time
    %   because sea surface height changes in time.
    %   thickness of each layer changes propotional to total water thicknes.

    h_total = g.h + zeta;       %total water thickness
    for level=1:N               %thickness of each layer
        Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
    end

    % average Hz to  Arakawa-C u points

    Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
    z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
    for k=2:+1:N
        z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
    end

    temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
    salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
    fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point

    % average Hz to  Arakawa-C v points

    Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
    z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
    for k=2:+1:N
        z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
    end

    temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
    salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
    fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point


    %   ====================================================================================
    %   find path from corner_endpt(1) to corner_endpt(2)

    icount=1;
    col_index(icount)=corner_endpt_col(1);
    row_index(icount)=corner_endpt_row(1);
    on_vpoint(icount)=0;
    vpoint=0;
    xpoint=g.lon_u( row_index(icount), col_index(icount) );
    ypoint=g.lat_u( row_index(icount), col_index(icount) );

    signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
    dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
    tmp_dist=dist(icount);
    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
    flag_approach = 1;

    %  icount
    %  row_index(icount)
    %  col_index(icount)
    %  vpoint
    %  signline(icount)
    %  dist2endpoint(icount)

    % plot(  xpoint,  ypoint , 'ro')


    if (west2east_transect)

        while ( dist2endpoint(icount) > avgdxdy  &&  flag_approach )

            icount=icount+1;

            if ( vpoint == 1 )

                col_index(icount)=col_index(icount-1)+delj;
                if ( on_vpoint(icount-1) == 1)
                    row_index(icount)=row_index(icount-1);
                else
                    row_index(icount)=row_index(icount-1)+deli;
                end
                xpoint=g.lon_v( row_index(icount), col_index(icount) );
                ypoint=g.lat_v( row_index(icount), col_index(icount) );
                signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                if ( signline(icount)*signline(icount-1) < 0  ...
                        ||   dist(icount) <= dist(icount-1)       ...
                        ||   dist(icount) <= tmp_dist                  )
                    tmp_dist=0;
                    on_vpoint(icount)=1;
                    %                 plot(  xpoint,  ypoint , 'ro')
                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
                else
                    tmp_dist=dist(icount);
                    vpoint=0;
                    icount=icount-1;
                    %                 plot(  xpoint,  ypoint , 'gx')
                end

            else % on upoint

                col_index(icount)=col_index(icount-1);
                if ( on_vpoint(icount-1) == 0)
                    row_index(icount)=row_index(icount-1)+deli;
                else
                    row_index(icount)=row_index(icount-1);
                end
                xpoint=g.lon_u( row_index(icount), col_index(icount) );
                ypoint=g.lat_u( row_index(icount), col_index(icount) );
                signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                if (      signline(icount)*signline(icount-1) < 0 ...
                        ||   dist(icount) <= dist(icount-1)   ...
                        ||   dist(icount) <= tmp_dist                       )
                    tmp_dist=0;
                    on_vpoint(icount)=0;
                    %                 plot(  xpoint,  ypoint , 'ro')
                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
                else
                    tmp_dist=dist(icount);
                    vpoint=1;
                    icount=icount-1;
                    %                 plot(  xpoint,  ypoint , 'gx')
                end

            end % if ( on_vpoint == 1 )

            if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                flag_approach = 0;
            end

            % icount
            % row_index(icount)
            % col_index(icount)
            % vpoint
            % signline(icount)
            % dist2endpoint(icount)

        end % while

    else % if (west2east_transect)

        while ( dist2endpoint(icount) > avgdxdy  &&  flag_approach )

            icount=icount+1;

            if ( vpoint == 1 )

                if ( on_vpoint(icount-1) == 1)
                    col_index(icount)=col_index(icount-1)+delj;
                    row_index(icount)=row_index(icount-1);
                else
                    col_index(icount)=col_index(icount-1);
                    row_index(icount)=row_index(icount-1)+deli;
                end

                xpoint=g.lon_v( row_index(icount), col_index(icount) );
                ypoint=g.lat_v( row_index(icount), col_index(icount) );
                signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                if ( signline(icount)*signline(icount-1) < 0  ...
                        ||   dist(icount) <= dist(icount-1)       ...
                        ||   dist(icount) <= tmp_dist                  )
                    tmp_dist=0;
                    on_vpoint(icount)=1;
                    %                 plot(  xpoint,  ypoint , 'ro')
                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
                else
                    tmp_dist=dist(icount);
                    vpoint=0;
                    icount=icount-1;
                    %                 plot(  xpoint,  ypoint , 'gx')
                end

            else % on upoint

                if ( on_vpoint(icount-1) == 0)
                    col_index(icount)=col_index(icount-1);
                    row_index(icount)=row_index(icount-1)+deli;
                else
                    col_index(icount)=col_index(icount-1)+delj;
                    row_index(icount)=row_index(icount-1);
                end


                xpoint=g.lon_u( row_index(icount), col_index(icount) );
                ypoint=g.lat_u( row_index(icount), col_index(icount) );
                signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                if (      signline(icount)*signline(icount-1) < 0 ...
                        ||   dist(icount) <= dist(icount-1)   ...
                        ||   dist(icount) <= tmp_dist                       )
                    tmp_dist=0;
                    on_vpoint(icount)=0;
                    %                 plot(  xpoint,  ypoint , 'ro')
                    dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone]);
                else
                    tmp_dist=dist(icount);
                    vpoint=1;
                    icount=icount-1;
                    %                 plot(  xpoint,  ypoint , 'gx')
                end

            end % if ( on_vpoint == 1 )

            if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                flag_approach = 0;
            end

            %icount
            % row_index(icount)
            % col_index(icount)
            % vpoint
            % signline(icount)
            % dist2endpoint(icount)

        end % while

    end % if (west2east_transect)


    total_volume=0;
    total_salt=0;
    total_fresh=0;

    for index=1:icount

        vpoint=on_vpoint(index);
        %% Subroutine ! ! !
        xy_transport_function

        if( west2east_transect == 0 && vpoint == 1 )
            total_volume = total_volume - sum_segment;
            total_salt   = total_salt   - sum_segment_salt;
            total_fresh  = total_fresh  - sum_segment_fresh;
            freshtransect(index)=-sum_segment_fresh;
            voltransect(index)=-sum_segment;
            % Peter added start
            salttransect(index)=-sum_segment_salt;
            % Peter added end
        else
            total_volume = total_volume + sum_segment;
            total_salt   = total_salt   + sum_segment_salt;
            total_fresh  = total_fresh  + sum_segment_fresh;
            freshtransect(index)=sum_segment_fresh;
            voltransect(index)=sum_segment;
            % Peter added start
            salttransect(index)=sum_segment_salt;
            % Peter added end
        end


    end

    if( west2east_transect )
        disp(' from West to East transect ')
    else
        disp(' from East to West (reverse) transect ')
    end
    disp(' ')

    disp([' volume transport = ', num2str(total_volume) ,' m^3/s '])
    disp(['                  = ', num2str(total_volume/1e+6) ,' x 10^6m^3/s '])

    disp([' salt volume transport = ', num2str(total_salt) ,' g/s '])
    disp(['                  = ', num2str(total_salt/1000) ,' x Kg/s '])

    disp([' freshwater volume transport = ', num2str(total_fresh) ,' m^3/s '])
    disp(['                  = ', num2str(total_fresh/1e+6) ,' x 10^6m^3/s '])

    disp([' transect distance = ', num2str( distance_r/1000 ),' km '])
    disp([' first  point, lon, lat = ',  num2str( xx(1) ),'  ',num2str( yy(1) )])
    disp([' second point, lon, lat = ',  num2str( xx(2) ),'  ',num2str( yy(2) )])

    D_total_vol(day)=total_volume;
    D_total_salt(day)=total_salt;
    D_total_fresh(day)=total_fresh;
end

% save salt_transport_2012_re.mat D_total_fresh D_total_salt D_total_vol
save salt_transport_2012_test.mat D_total_fresh D_total_salt D_total_vol
iplot=0;
if (iplot)
    figure
    subplot(4,1,1)
    plot([1:icount],voltransect/100,'g-x')
    xlabel('index')
    ylabel('volume trans')
    grid
    subplot(4,1,2)
    plot([1:icount],freshtransect/100,'g-x')
    xlabel('index')
    ylabel('freshwater trans')
    grid
    % Peter added start
    subplot(4,1,3)
    plot([1:icount],salttransect/100,'g-x')
    xlabel('index')
    ylabel('salt trans')
    grid
    % Peter added end
    subplot(4,1,4)
    plot([1:icount],on_vpoint,'ro')
    xlabel('index')
    legend('1 for vpoint, 0 for upoint')
    grid
end


clf
subplot(311)
plot(1:365,D_total_vol,'r')
hold on
plot([0 365],[0 0],'k')
axis([1 365 -3e6 3e6])

subplot(312)
plot(1:365,D_total_salt,'r')
hold on
plot([0 365],[0 0],'k')
axis([1 365 -1e8 1e8])

subplot(313)
plot(1:365,D_total_fresh,'r')
hold on
plot([0 365],[0 0],'k')
axis([1 365 -25e4 25e4])

vol=mean(D_total_vol);
salt=mean(D_total_salt);
fresh=mean(D_total_fresh);

disp([' volume transport = ', num2str(vol) ,' m^3/s '])
disp(['                  = ', num2str(vol/1e+6) ,' x 10^6m^3/s '])

disp([' salt volume transport = ', num2str(salt) ,' g/s '])
disp(['                  = ', num2str(salt/1000) ,' x Kg/s '])

disp([' freshwater volume transport = ', num2str(fresh) ,' m^3/s '])
disp(['                  = ', num2str(fresh/1e+6) ,' x 10^6m^3/s '])


break

clf
subplot(311)
plot(152:243,D_total_vol(152:243),'r')
hold on
plot([152 243],[0 0],'k')
axis([152 243 -3e6 3e6])

subplot(312)
plot(152:243,D_total_salt(152:243),'r')
hold on
plot([152 243],[0 0],'k')
axis([152 243 -1e8 1e8])

subplot(313)
plot(152:243,D_total_fresh(152:243),'r')
hold on
plot([152 243],[0 0],'k')
axis([152 243 -25e4 25e4])

break
save salt_transport_2010.mat D_total_fresh D_total_salt D_total_vol