function startup

% startup -- User script configuration for Matlab.  It can set default
%            paths, define Handle Graphics defaults, or predefine
%            variables in your workspace.

% svn $Id: startup.m 711 2014-01-23 20:36:13Z arango $
%===========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

global IPRINT
IPRINT=0;

format long g

% Change "my_root" to the appropriate path were these matlab scripts are
% installed in your computer.

% my_home = getenv('HOME');
% my_root = strcat(my_home, '/ocean/repository');

my_root = strcat('D:/JNUpack/JNUROMS/Preprocessing/Clm_nudg/input2020/rutgers_mat/');
path(path, fullfile(my_root, 'process', '4dvar', ''))
path(path, fullfile(my_root, 'process', 'bathymetry', ''))
path(path, fullfile(my_root, 'process', 'boundary', ''))
path(path, fullfile(my_root, 'process', 'coastlines', ''))
path(path, fullfile(my_root, 'process', 'forcing', ''))
path(path, fullfile(my_root, 'process', 'grid', ''))
path(path, fullfile(my_root, 'process', 'initial', ''))
path(path, fullfile(my_root, 'process', 'landmask', ''))
path(path, fullfile(my_root, 'process', 'mex', ''))
path(path, fullfile(my_root, 'process', 'netcdf', ''))
path(path, fullfile(my_root, 'process', 'seagrid', ''))
path(path, fullfile(my_root, 'process', 'seagrid', 'presto', ''))
path(path, fullfile(my_root, 'process', 'seawater', ''))
path(path, fullfile(my_root, 'process', 't_tide', ''))
path(path, fullfile(my_root, 'process', 'tidal_ellipse', ''))
path(path, fullfile(my_root, 'process', 'utility', ''))

% addpath '/home/jhlee/Matlab_tools/roms_tools/mexnc'
% addpath '/home/jhlee/Matlab_tools/roms_tools/mexnc/tests'

% Load NetCDF Toolbox for OpenDAP support for versions 2008b or higher. 
% However, this is not needed if version 2012a or higher since Matlab
% native NetCDF interface supports OpenDAP.  Users need to change the
% paths for SNCTOOLS and JAVA.

% v = version('-release');
% vyear = str2num(v(1:4));
% load_toolbox = vyear >= 2008;
% if ((vyear == 2008 && v(5:5) == 'a') || vyear >= 2012),
%   load_toolbox = false;
% end
% 
% if (load_toolbox),
%   addpath (strcat(my_home, '/ocean/matlab/snctools'), '-end');
%   javaaddpath (strcat(my_home, '/ocean/matlab/classes/toolsUI-4.1.jar'), '-end');
%   javaaddpath (strcat(my_home, '/ocean/matlab/classes/netcdfAll-4.2.jar'), '-end');
%   javaaddpath (strcat(my_home, '/ocean/matlab/snctools/classes'), '-end');
%   setpref('SNCTOOLS','USE_JAVA', true);
% end