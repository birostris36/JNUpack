function han = roms_addvect_scale(pos,uv,uscale,label,varargin)
%  han = roms_addvect_scale(pos,uv,uscale,label,varargin)
%
%  pos = [x y] position for scale vector (can use, e.g. ginput(1)
%  uv = [u v] scale vector
%  uscale = same scale as used in the roms_quiver command
%  label = string to label vector, default is 'm/s'
%  varargin = arguments to quiver, e.g. vector color
%
% John Wilkin
  
% add a scale vector to the plot

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

hq = quiver(pos(1),pos(2),uscale*uv(1),uscale*uv(2),0,varargin{:});

if nargin < 4
  label = 'm/s';
end
if isempty(label)
  label = 'm/s';
end

% interactive text placement
% ht = gtext([num2str(uscale) ' m/s']);

umag = abs(uv(1)+sqrt(-1)*uv(2));
%text_string=num2str(umag,'%-5.2e');
text_string=num2str(umag,'%-5.0e');
switch( str2num(text_string(end-2:end)) )
  case(-1)
  new_string=[ text_string(1:end-4) 'x10^{-1}' ];
  case(-2)
  new_string=[ text_string(1:end-4) 'x10^{-2}' ];
  case(-3)
  new_string=[ text_string(1:end-4) 'x10^{-3}' ];
  case(-4)
  new_string=[ text_string(1:end-4) 'x10^{-4}' ];
  case(-5)
  new_string=[ text_string(1:end-4) 'x10^{-5}' ];
  case(-6)
  new_string=[ text_string(1:end-4) 'x10^{-6}' ];
  case(-6)
  new_string=[ text_string(1:end-4) 'x10^{-7}' ];
  case(1)
  new_string=[ text_string(1:end-4) 'x10^{1}' ];
  case(2)
  new_string=[ text_string(1:end-4) 'x10^{2}' ];
  case(3)
  new_string=[ text_string(1:end-4) 'x10^{3}' ];
  case(4)
  new_string=[ text_string(1:end-4) 'x10^{4}' ];
  case(5)
  new_string=[ text_string(1:end-4) 'x10^{5}' ];
  case(6)
  new_string=[ text_string(1:end-4) 'x10^{6}' ];
  case(7)
  new_string=[ text_string(1:end-4) 'x10^{7}' ];
  otherwise
  new_string=text_string;
end
%ht = text(pos(1),pos(2),[text_string ' ' label ' ']);
ht = text(pos(1),pos(2),[new_string ' ' label ' ']);
set(ht,'horizontalalignment','right')

if nargout > 0
  han = [hq];
  if exist('ht')==1
    han = [han; ht];
  end
end

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);

