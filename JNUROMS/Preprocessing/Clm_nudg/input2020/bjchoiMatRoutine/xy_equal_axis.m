function xy_equal_axis(centerlat)
%
%   make equal axis
%

    %centerlat = 44;
    ax_ratio=cos( centerlat*pi/180 );
    set(gca,'DataAspectRatio',[1 ax_ratio 1])
    
return
