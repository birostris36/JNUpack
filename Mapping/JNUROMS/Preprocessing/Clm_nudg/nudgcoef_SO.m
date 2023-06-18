%
%  D_NUDGCOEF:  Driver script to create a ROMS nudging coefficients file.
%
%  This a user modifiable script that can be used to prepare ROMS
%  nudging inverse time scales NetCDF file.  It sets-up  all the
%  necessary parameters and variables. USERS can use this as a
%  prototype for their application.
%
%  Nudging to climatology can be used in ROMS for various purposes:
%
%  (1) Improve the behavior of open boundary conditions.
%  (2) Used in conjunction with sponges.
%  (3) Minimize numerical diapycnal mixing of tracers over steep
%      bathymetry (improve T-S properties in deep water masses).  For
%      example, we can nudge to T-S climatology is areas depeer than
%      1500 m.
%
%  The inverse nudging coefficients have units of 1/time.  The default
%  input units in ROMS is 1/day but 1/second is also possible. The
%  routine 'get_nudgcoef.F' will check the 'units' attribute to compute
%  the conversion factor for 1/second. Users need to be sure that the
%  'units' variable attribute is consistent with the written data.
%
%  The variable names for the nudging coefficients is as follows:
%
%     M2_NudgeCoef       for 2D momentum
%     M3_NudgeCoef       for 3D momentum
%     temp_NudgeCoef     for potential temperature
%     salt_NudgeCoef     for salinity
%     ...
%     NO3_NudgeCoef      for nitrate
%     ...
%     tracer_NudgeCoef   for any generic tracer
%
%  They are all defined at RHO-points. If the nudging coefficients for
%  a specific tracer are available in the NetCDF, ROMS will read that
%  NetCDF variable. If NOT and the generic coefficients 'tracer_NudgeCoef'
%  are available, ROMS will process those values instead.
%
%  Notice that the input swicth 'LnudgeTCLM(itrc,ng)' in ROMS input
%  script 'ocean.in' will control which tracer to nudge in the desired
%  grid.
%
%  Currently, the nudging coefficients are time invariant in ROMS.  The
%  same scales are used for the entire simulation.
%

% svn $Id: d_nudgcoef.m 832 2017-01-24 22:07:36Z arango $
%=========================================================================%
%  Copyright (c) 2002-2017 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
% Set input/output NetCDF files.
% 
% restoredefaultpath
% startup

 my_root = '/data4/base158/Warehouse01/';

 GRDname = fullfile(['/data4/base158/Warehouse01/', 'Grd_SO_05d_sponge.nc']);
%  INIname = fullfile(my_root, 'JEJU2_ini_HYCOM_Y2017M01D01.nc');
 INIname = fullfile(['/data4/base158/Warehouse01/soda/', 'Ini_soda_05d_jhlee_198002.nc']);
 
  NUDname = [my_root 'NUDG_05d_SO.nc'];
 
% Get grid structure.

G = get_roms_grid(GRDname, INIname);

spherical = G.spherical;

[Lr,Mr] = size(G.h);
% Nr = length(G.s_rho);
Nr=50;
% Nr=20;

% Set switches for state variables to nudge.

LnudgeM2CLM    = false;           % nudging 2D momentum
LnudgeM3CLM    = true;           % nudging 3D momentum
LnudgeTCLM     = true;           % nudging tracers (usually T-S)
LnudgeTgeneric = true;           % nudging generic tracers

% Set NetCDF variables to process.  Initialize inverse nudging
% coefficients with zeros.  Recall that division by zero is not
% defined and we will give "Inf".  Therefore, we just need to set
% only the values in the desired areas and the rest can be zero
% (for no nudging) because the nudging in ROMS is:
%
%      F(...,new) = F(...,new) +
%                   dt * F_nudgcoef * (Fclm - F(...,new))

VarNUD = [];
F = [];

if (spherical),
  VarNUD = {'lon_rho', 'lat_rho'};
else
  VarNUD = {'x_rho', 'y_rho'};
end

if (LnudgeM2CLM),
  VarNUD = [VarNUD, 'M2_NudgeCoef'];
  F.M2_NudgeCoef = zeros(Lr,Mr);                % RHO-points
end

if (LnudgeM3CLM),
  VarNUD = [VarNUD, 'M3_NudgeCoef'];
  F.M3_NudgeCoef = zeros(Lr,Mr,Nr);             % RHO-points
end

if (LnudgeTCLM),
  VarNUD = [VarNUD, 'temp_NudgeCoef', 'salt_NudgeCoef'];
  F.temp_NudgeCoef = zeros(Lr,Mr,Nr);
  F.salt_NudgeCoef = zeros(Lr,Mr,Nr);
end

if (LnudgeTgeneric),
  VarNUD = [VarNUD, 'tracer_NudgeCoef'];
  F.tracer_NudgeCoef = zeros(Lr,Mr,Nr);
end

%--------------------------------------------------------------------------
% Create Nudging coefficients NetCDF file: build creation parameters
% structure, S.
%--------------------------------------------------------------------------

S            = [];
Tindex       = [];
ReplaceValue = NaN;
PreserveType = true;
Unlimited    = false;                   % time dimension is umlimited
nctype       = 'nc_double';             % input data is in double precision

mode = netcdf.getConstant('CLOBBER');                    % overwrite!!!
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

% The strategy here is to build manually the NetCDF metadata structure to
% facilitate creating several NetCDF file in a generic and compact way.
% This structure is similar to that returned by "nc_inq" or native Matlab
% function "ncinfo".
%
% Notice that we call the "roms_metadata" function to create the fields
% in the structure.  Then, we call "check_metadata" for fill unassigned
% values and to check for consistency.

ncname = NUDname;

disp(blanks(1));
disp(['** Creating NetCDF file: ', ncname,' **']);
disp(blanks(1))

S.Filename = ncname;

S.Attributes(1).Name      = 'type';
S.Attributes(1).Value     = 'Nudging Coeffcients file';

S.Attributes(2).Name      = 'title';
S.Attributes(2).Value     = ['North Atlantic Damee #4, 0.75 Resolution'];

S.Attributes(3).Name      = 'grd_file';
S.Attributes(3).Value     = GRDname;

S.Attributes(4).Name      = 'ini_file';
S.Attributes(4).Value     = INIname;

S.Attributes(5).Name      = 'history';
S.Attributes(5).Value     = ['Nudging coefficient file created from ',  ...
                             'd_nudgcoef.m: ', date_stamp];

S.Dimensions(1).Name      = 'xi_rho';
S.Dimensions(1).Length    = Lr;
S.Dimensions(1).Unlimited = false;

S.Dimensions(2).Name      = 'eta_rho';
S.Dimensions(2).Length    = Mr;
S.Dimensions(2).Unlimited = false;

S.Dimensions(3).Name        = 's_rho';
S.Dimensions(3).Length      = Nr;
S.Dimensions(3).Unlimited   = false;

S.Variables(1) = roms_metadata('spherical');
if (spherical),
  S.Variables(2) = roms_metadata('lon_rho');
  S.Variables(3) = roms_metadata('lat_rho');
else
  S.Variables(2) = roms_metadata('x_rho');
  S.Variables(3) = roms_metadata('y_rho');
end

% Process inverse nudging coefficient variables.

for n = 3:length(VarNUD),

  Vname  = char(VarNUD{n});
  
  i = n + 1;
  S.Variables(i) = roms_metadata(Vname, spherical, nctype, Unlimited);

end
  
% Check ROMS metadata structure.  Fill unassigned fields.

S = check_metadata(S);
  
% Create forcing NetCDF files.  Write grid coordinates.

ncid = nc_create(S.Filename, mode, S);    % create a new NetCDF file

% status = nc_write(S.Filename, 'spherical', 'T');

status = nc_write(S.Filename, 'spherical', int32(spherical));
status = nc_write(S.Filename, 'lon_rho',   G.lon_rho);
status = nc_write(S.Filename, 'lat_rho',   G.lat_rho);

%--------------------------------------------------------------------------
% Set inverse time scales.
%--------------------------------------------------------------------------

IstrR = 0;
IendR = Lr-1;
JstrR = 0;
JendR = Mr-1;

% In the North Atlantic DAMEE_4 application the nudging is done in the
% southern and northern domain edges over a 8-point linearly tapered
% nudging scales of 5 to 60 days.

inner1 = 1/180;                        % 60 days at interior limit
outer1 = 1/2;                         % 1.5 days at boundary
width1 = 25;                           % 25 points
inner2 = 0;                        % 60 days at interior limit
outer2 = 0;                         % 1.5 days at boundary
width2 = 0;                           % 20 points


work  = zeros(Lr,Mr);

% for j=JstrR: width1,                   % Southern boundary
%   for i=IstrR+width2:IendR -width2,
%     work(i+1,j+1) = inner1 + (width1 - j) * (outer1 - inner1) / width1;
%   end
% end

for j=JendR-width1:JendR,             % Northern boundary
  for i=IstrR + width2:IendR - width2,
    work(i+1,j+1) = outer1 + (JendR - j) * (inner1 - outer1) / width1;
  end
end


% ii = 0;               
% for i=IendR-width1-width2:IendR-width2,        % Eastern boundary1
%     
%    for j=JstrR+width1-ii:JendR-width1+ii,
% %    for j=JstrR+ii:JendR-ii,
%     work(i+1,j+1) = outer1 + (IendR - i - width2) * (inner1 - outer1) / width1;
% 
%    end
%   ii = ii + 1;
% end
% 
% ii = 0;
% for i=IstrR+width2:width1+width2,                   % Western boundary1
%     
%   for j=JstrR+ii:JendR-ii,
% 
%     work(i+1,j+1) = inner1 + (width1 - i + width2) * (outer1 - inner1) / width1;
% 
%   end
%   ii = ii + 1;
% end


% ii = 0;               
% for i=IendR-width2:IendR,        % Eastern boundary2
%     
%    for j=JstrR:JendR,
% %    for j=JstrR+ii:JendR-ii,
%     work(i+1,j+1) = outer2 + (IendR - i) * (inner2 - outer2) / width2;
% 
%    end
%   ii = ii + 1;
% end
% % 
% ii = 0;
% for i=IstrR:width2,                   % Western boundary2
%     
%   for j=JstrR:JendR,
% 
%     work(i+1,j+1) = inner2 + (width2 - i) * (outer2 - inner2) / width2;
% 
%   end
%   ii = ii + 1;
% end

% 
% fac = (7 * outer - inner)/6;
% 
% for j=74:80,                         % Mediterranean outflow
%   for i=102:106,
%     cff = sqrt((i-109)^2 + (j-77)^2);
%     work(i+1,j+1) = max(0, (fac + cff * (inner - outer) / 6));
%   end
% end

% Load values into structure

if (LnudgeM2CLM),
  F.M2_NudgeCoef = work;
end

if (LnudgeM3CLM),
  F.M3_NudgeCoef = repmat(work, [1 1 Nr]);
end

if (LnudgeTCLM),
  F.temp_NudgeCoef = repmat(work, [1 1 Nr]);
  F.salt_NudgeCoef = repmat(work, [1 1 Nr]);
end

if (LnudgeTgeneric),
  F.tracer_NudgeCoef = repmat(work, [1 1 Nr]);
end

%--------------------------------------------------------------------------
% Write out inverse nudging coefficients.
%--------------------------------------------------------------------------

Nfields = length(VarNUD);

for i=3:Nfields,
  field = char(VarNUD(i));
  status = nc_write(NUDname, field, F.(field));
end
