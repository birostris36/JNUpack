function han = plotlattecoast(option,varargin)
% han = plotnenacoast(option)
% plots Northeast North American coast
%
% abs(option) = 0,1 gives intermediate resolution            
%             = 2 gives highest resolution coastline
%             = 3 is same as 2 without the lakes
%
% option < 0 changes axis limits to encompass whole dataset 
%             (useful for getting quick map)

% get plot state
nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted
set(gca,'nextplot','add')

if nargin == 0
   % load('/home/wilkin/roms/nena/mat/grid/natl_coast_int')
     load('/home/caledonia/bjchoi/roms/hudson/grid/natl_coast_int')
   % load('/home/caledonia/bjchoi/roms/hudson/grid/natl_coast_low')

else
    switch abs(option)
      case 3
	     %load('useast_coast_high_3')
              load('/home/bchoi/roms/hudson/grid/useast_coast_high')
      case 2
	     % load('/home/wilkin/roms/nena/mat/grid/useast_coast_high')
	     % load('useast_coast_high')
               load('/home/bchoi/roms/hudson/grid/useast_coast_high')
      case {1,0}
	     % load('natl_coast_int')
	     % load('/home/wilkin/roms/nena/mat/grid/natl_coast_int')
             load('/home/bchoi/roms/hudson/grid/useast_coast_high')
    end
    if option<0
      axis([min(lon(:)) max(lon(:)) min(lat(:)) max(lat(:))])
  end
end
  
ax = axis;
lon(find(lon<ax(1)|lon>ax(2))) = NaN;
lat(find(lat<ax(3)|lat>ax(4))) = NaN;
h = plot(lon,lat,varargin{:},'k','LineWidth',0.1);

% restore nextplotstate to what it was
set(gca,'nextplot',nextplotstatewas);

if nargout > 0
   han = h;
end
