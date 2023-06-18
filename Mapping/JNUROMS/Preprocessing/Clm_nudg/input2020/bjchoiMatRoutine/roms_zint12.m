function [vint,thickness,dz] = roms_zint12(data,grd,z1,z2)
% vertical integral of variable 'data' at rho points from depth z1 to z2
%
% If inputs z1,z2 are absent (nargin==2) it is assumed the integral is over
% the whole water column.
%
% If z1 == 0 or [] it is assumed the integral is from the surface to z2
% If z2 == NaN or [] it is assumed the integral is from z1 to -h
%
% If z1 or z2 are 2-d arrays is it assumed they described a matrix of
% depths (negative) defined on the rho points grid, such as the depth of
% some property (mixed layer, an isopycnal surface) determined elseewhere.
% z1 could be zeta and z2 could be -h, in which we just get the integral
% over the full water column.
%
% Optional second output is the thickness of the layer, so that a vertical
% average can be calculated if required.
%
% Method is to compute a matrix of layer thicknesses, then mutliply the
% data element-by-element, then sum over the vertical dimension
%
% To integrate over depths not the full water column, the layer thickness
% matrix is modified accordingly by effectively implementing
%
% int_z2^z1 dz = int_-h^0 dz - int_z1^0 dz - int_-h^z2 dz

% int_-h^0 dz first:
% layer thicknesses
dz = diff(grd.z_w,1,1); 

% interpret the vertical limits to see if we are asked to apply special
% limits z1 or z2 not equal to 0 and -h, respectively

if nargin > 2
    if ~isempty(z1) & z1~=0
        % then adjust dz by -int_z1^0
        dz = dz - diff(max(grd.z_w,z1),1,1);
    end
end
if nargin > 3
    if ~isempty(z2) & ~isnan(z2)
        % determine whether z2 is a constant or a matrix
        if prod(size(z2)==1
            % simple constant
            dz = dz - diff(min(grd.z_w,z2),1,1);
        else
            if ndims(z2)~=2
                error('z2 is not a 2D matrix')
            elseif ~all(size(z2)==size(grd.h))
                error('z2 is the wrong shape - should be same as h')
            else
                % z2 is a 2D matrix
                z2 = reshape(repmat(z2,[grd.N+1 1]),size(z_w));
                dz = dz - diff(min(grd.z_w,z2),1,1);
            end
        end
    end
end

% data times dz
vint = data.*dz;

% vertical integral
vint = squeeze(sum(vint,1));

% thickness is sum of dz
thickness = squeeze(sum(dz,1));