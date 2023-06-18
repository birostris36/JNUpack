function [tmpsum,tmpvar]=roms_loadvarsum(ctl,ind,field,field2,grd)
% (R)oms (N)umerical (T)oolbox
%
% [fieldsum, fieldvarsum]=roms_loadvarsum(ctl,ind,field,field2)
%
% Extract the sum of the variable FIELD for the indicies specified
% in the array IND from a composite netcdf file which is defined 
% in the time controll arrays CTL.
%
% SEE rnt_timectl.m to generate a the CTL struct. array 
%      (it is easy ..no worry ok!)
%
% INPUT: example
%     ctl = rnt_timectl(files,timevar);
%     ind = [ 1:6 ] (get time indiceis 1:6 from composite field
%     field = 'temp'
% 
% OUTPUT:
%   fieldsum(x,y,z) = sum ( theField(x,y,z,ind) ) over ind
%     if you want the mean just reassign:
%     fieldmean = fieldsum / lenght(ind) 
%
%   fieldvar is the sum of the variances
%
% RNT - E. Di Lorenzo (edl@ucsd.edu)
%
% This version modified by John Wilkin from EDL's rnt_loadvarsum to
% retain the natural C-language style index order returned by nc.
%
% If field = 'ut' or 'vt' then tracer fluxes are computed prior to
% averaging where the tracer variable is given by input FIELD2
% (which defaults to 'temp' if not given in the argument list)


% The computation of variances is still disabled.

tmpsum = 0;
% find info about the variabe and initialize arrays
if isempty(ind)
  disp(['roms_loadvarsum - no time index match for ',field]);
  nc = netcdf(ctl.file{1});
  s = ncsize(nc{field});
  close(nc)
  s(1) = 1;
  tmpsum = zeros(s);
  tmpvar = tmpsum;
  return
end

% compute mean
% the sum for each file in the sequence is calculated separately, then
% accumulated 
for istep=1:length(ctl.segm)-1
  in = find ( ind>ctl.segm(istep) & ind<=ctl.segm(istep+1));
  in_extr = ctl.ind(ind(in));

  if ~isempty(in_extr)  
    nc = netcdf(ctl.file{istep});    
    
    if ~isempty(nc{field})
      % requested variable is in the file
      tmp2 = nc{field}(in_extr,:);
    else
      % requested variable is not in the file
      % (function roms_cgridpos is available to possibly make this
      % section of code more general by detecting the shape of u and v
      % variables) 
      switch lower(field)
	case 'ut'
	  if nargin < 4
	    field2 = 'temp';
	  end	  
	  tmpuv = nc{'u'}(in_extr,:);
	  tmptr = nc{field2}(in_extr,:);
	  tmptr = 0.5*(tmptr(:,:,:,1:end-1)+tmptr(:,:,:,2:end));
	  tmp2 = tmpuv.*tmptr;
	case 'vt'
	  if nargin < 4
	    field2 = 'temp';
	  end	  
	  tmpuv = nc{'v'}(in_extr,:);
	  tmptr = nc{field2}(in_extr,:);
	  tmptr = 0.5*(tmptr(:,:,1:end-1,:)+tmptr(:,:,2:end,:));
	  tmp2 = tmpuv.*tmptr;
	case 'dzdyut'
	  % product of cell face area dz*dy and flux u*tracer
	  % work in progress
	  % need to embed a local loop over in_extr times because zeta
	  % needs to be extraced on each step
	  if nargin < 4
	    field2 = 'temp';
	  end
	  for i = 1:length(in_extr)
	    ix = in_extr(i);
	    % dz*dy at u-points
	    %
	    % make this into a function later --------
	    % z_w coordinates (depend on zeta)
	    zeta = nc{'zeta'}(ix,:);
	    scmCshc_w = (grd.sc_w-grd.Cs_w)*grd.hc;
	    z_w = repmat(scmCshc_w,[1 length(grd.h(:))])...
		+grd.Cs_w*grd.h(:)';
	    z_w = z_w + scmCshc_w*[zeta(:)./grd.h(:)]'...
		+ (1+grd.Cs_w)*zeta(:)';
	    grd.z_w = reshape(z_w,[grd.N+1 size(h)]);
	    % average to u-points
	    z_wu = 0.5*(grd.z_w(:,:,1:end-1)+grd.z_w(:,:,2:end));
	    n_u  = 0.5*(grd.pn(:,1:end-1)+grd.pn(:,2:end));
	    h_u  = 0.5*(grd.h(:,1:end-1)+grd.h(:,2:end));
	    dzon_u = diff(z_wu,1)./...
		permute(repmat(n_u,[1 1 grd.N]),[3 1 2]);
	    % ----------------------------------------
	    tmpuv = nc{'u'}(ix,:);
	    tmptr = nc{field2}(ix,:);
	    tmptr = 0.5*(tmptr(:,:,:,1:end-1)+tmptr(:,:,:,2:end));
	    tmp2(i,:) = dzon_u.*tmpuv.*tmptr;
	  end
	case 'dzdxvt'
	  
	  % work in progress
	  % not sure about need for more squeezes
	  % check case of only one time level in in-extr
		  
	  % product of cell face area dz*dx and flux v*tracer
	  %
	  % need to embed a local loop over in_extr times because zeta
	  % needs to be extraced on each step
	  if nargin < 4
	    field2 = 'temp';
	  end
	  for i = 1:length(in_extr)
	    ix = in_extr(i);
	    % dz*dx at v-points
	    %
	    % make this into a function later --------
	    % z_w coordinates (depend on zeta)
	    zeta = nc{'zeta'}(ix,:);
	    scmCshc_w = (grd.sc_w-grd.Cs_w)*grd.hc;
	    z_w = repmat(scmCshc_w,[1 length(grd.h(:))])...
		+grd.Cs_w*grd.h(:)';
	    z_w = z_w + scmCshc_w*[zeta(:)./grd.h(:)]'...
		+ (1+grd.Cs_w)*zeta(:)';
	    grd.z_w = reshape(z_w,[grd.N+1 size(h)]);
	    % average to v-points
	    z_wv = 0.5*(grd.z_w(:,1:end-1,:)+grd.z_w(:,2:end,:));
	    m_v  = 0.5*(grd.pm(1:end-1,:)+grd.pm(2:end,:));
	    h_v  = 0.5*(grd.h(1:end-1,:)+grd.h(2:end,:));
	    dzom_v = diff(z_wv,1)./...
		permute(repmat(m_v,[1 1 grd.N]),[3 1 2]);
	    % ----------------------------------------
	    tmpuv = nc{'v'}(in_extr,:);
	    tmptr = nc{field2}(in_extr,:);
	    tmptr = 0.5*(tmptr(:,:,1:end-1,:)+tmptr(:,:,2:end,:));
	    tmp2(i,:) = dzom_v.*tmpuv.*tmptr;
	  end	  
	otherwise
	  errorstr = [ 'Variable ' field ' is not present in ' ...
		ctl.file{istep} ', nor is there any rule for ' ...
		'deriving it in ' which(mfilename)];
	  error(errorstr)
      end      
    end
    if length(in_extr)>1
      tmp2 = squeeze(sum(tmp2,1));
    end
    close(nc);
    % accumulate the sum for this segment
    tmpsum = tmpsum + tmp2;
  end  
end

%compute variance
if nargout > 1
  tmpmean = tmpsum/length(ind);
  tmpvar = zeros(size(tmpmean));
  for istep=1:length(ctl.segm)-1
    in = find ( ind>ctl.segm(istep) & ind<=ctl.segm(istep+1));
    in_extr = ctl.ind(ind(in));  
    if ~isempty(in_extr)
      nc=netcdf(ctl.file{istep});
      ii=length(in_extr);

      % this is broken - need to determine whether 2D or 3D variable
      
      % if length(ncsize(nc{field}))==2
      %	tmp2 = repmat(tmpmean,[ii 1 1 ]) - nc{field}(in_extr,:) ;
      %     else
      tmp3 = nc{field}(in_extr,:);
      % tmp2 = tmp3 - reshape(repmat(tmpmean,[ii 1 1]),size(tmp3));
      tmp2 = tmp3;
      %    end    
      tmpvar = tmpvar + squeeze(sum(tmp2.^2,1));
      close(nc);
    end  
  end
end

