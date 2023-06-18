function romstofli(varname,file,vover,frame_size,fli_file,frame_prefix)
% roms2fli(varname,file,vover,frame_size,fli_file,frame_prefix)
%
% e.g. roms2fli('salt','eauc_his_3.nc',0,[400 400])
%
% read a ROMS history file and animate surface values
% 
% Inputs:
%       varname = variable (e.g. 'temp')
%       file = netcdf history file (default roms_his.nc)
%       vover = structure:
%         vover.scale > 0 causes vectors to be plotted over image
%                  and scaled by the factor vover.scale
%         vover.decimate subsamples the vectors (for clarity)
%         vover.angle rotates vectors in the case x-axis is not east
%                  (you get these values from the grd_file)
%       frame_size (in pixels, default = [200 200])
%       fli_file = name for fli animation (default = VARNAME.fli)
%         If empty the fli-creation is skipped
%       frame_prefix = prefix for ppm frames to be combined into fli
%
% John Wilkin

if nargin < 1
  varname = 'temp';
end
if nargin < 2
  % default history file name
  file = 'roms_his.nc';
end
if nargin < 3
  % option to plot vectors over
  vover.option = 0;
end
if nargin < 4
  frame_size = [200 200];
end
if isempty(frame_size)
  frame_size = [200 200];
end
if nargin < 5
  fli_file = [varname '.fli'];
end
if nargin < 6
  % prefix for ppm file names  
  frame_prefix = 'frame';  
end

% open the netcdf file (requires Chuck Denham's toolbox)
nc = netcdf(file);

% get the coordinates
days = nc{'ocean_time'}(:)/86400;
if lower(nc{'spherical'}(:)) == 't'
  x_rho = nc{'lon_rho'}(:);
  y_rho = nc{'lat_rho'}(:);
else
  x_rho = nc{'x_rho'}(:);
  y_rho = nc{'y_rho'}(:);
  x_rho = x_rho/1000; % convert to metres
  y_rho = y_rho/1000;
end

% find the index of the surface vertical level
tmp = nc('s_rho');
s_surf = tmp(:);

% PLOT EACH FRAME AND CONVERT TO PPM BITMAP

% set the frame details (filename prefix, size, position)

display_area = [int2str(frame_size(1)) 'x' int2str(frame_size(2))];
frame_pos = [1500 1100 0 0]+[-frame_size frame_size];
set(gcf,'position',frame_pos)

if varname == 'zeta'
  data = nc{varname}(1,:,:);
else  
  data = nc{varname}(1,s_surf,:,:);  
end
cax = range(data);

for irec=1:length(days)
  
  filestr = strrep(file,'_','\_');
  tmp = length(filestr);
  if tmp > 25
    tmp = [max(1,tmp-25):tmp];
    filestr = ['...' filestr(tmp)];
  end
  infostr = [varname ' at day ' num2str(days(irec))];
  titlestr{1} = infostr;
  titlestr{2} = filestr;
  disp(infostr)

  % get all the data
  if varname == 'zeta'
    data = nc{varname}(irec,:,:);
  else  
    data = nc{varname}(irec,s_surf,:,:);  
  end
  
  pcolor(x_rho,y_rho,squeeze(data));
  shading flat
  
  if vover.option > 0
    % plot surface velocity vectors
    vax = vover.option;
    vsub = vover.decimate;
    angle = vover.angle;
    grd = roms_get_grid(nc.grd_file(:));
    angle = grd.angle;
    udata = vax*nc{'u'}(irec,N,:,:);
    vdata = vax*nc{'v'}(irec,N,:,:);
    hold on
    roms_quiver(x_rho,y_rho,udata,vdata,angle,vsub,'k');
    hold off
  end

  if lower(nc{'spherical'}(:)) == 't'
    % set DataAspectRatio to approx Mercator proportions at figure centre
    ylim = get(gca,'ylim');
    set(gca,'DataAspectRatio',[1 cos(mean(ylim)*pi/180) 1]);
  else
    axis equal
    axis tight
  end
  
  caxis(cax);
  colorbar
  hant = title(titlestr);
  post = get(hant,'pos');
  pos_new = [min(get(gca,'xlim')) post([2 3])];
  set(hant,'HorizontalAlignment','left','Position',pos_new,'Fontsize',14)

  drawnow

  % create a ppm bitmap of the image
  frame_name = [frame_prefix int2str(irec) '.ppm'];
  ppmwrite(frame_name); 

end

% combine the frames into a fli animation

if ~isempty(fli_file)
  ppm2fli(frame_prefix,fli_file,display_area);
end





