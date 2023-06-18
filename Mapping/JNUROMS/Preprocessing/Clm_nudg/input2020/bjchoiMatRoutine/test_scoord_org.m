%
% following transformation from s-coordinate to z-coordinate
% was adapted from set_depth.F
%

% if exist('g') ~= 1
   %g=grd('hudson');
   %g=grd('ctz3km');
   g=grd('ctz3km40');
% end

 sc=g.sc_r;
 Cs=g.Cs_r;
 hc=g.hc;  
 x=zeros(length(sc),1);

% =============================

 zeta=0.0;
 figure
 h=50;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x,z,'ro')

 hold on
 h=500;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x+1,z,'kx')

 hold on
 h=1000;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x+2,z,'bx')

 hold on
 h=3000;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x+3,z,'bx')

 axis([-1 4 -h 1])
 grid on
 title(' for different bottom depth: 50 500 1000 3000 meter ')


 axis([-1 4 -50 0])
%=========================

 h=15;

 figure
 zeta=0.0;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x,z,'ro')

 hold on
 zeta=+0.5;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x+1,z,'kx')

 hold on
 zeta=-0.5;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x-1,z,'bx')

 axis([-3 3 -h 1])
 grid on
 title(' for different zeta: -0.5 0.0 +0.5 meter')


% =============================

 zeta=-.50;
 figure
 h=15;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x,z,'ro')

 hold on
 h=10;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x+1,z,'kx')

 hold on
 h=40;
 z =    sc.* hc         + Cs.* (h-hc) ...
      + sc.*(hc*zeta/h) + Cs.*((h-hc)*zeta/h) ...
      + zeta;
 plot(x-1,z,'bx')

 axis([-3 3 -h 1])
 grid on
 title(' for different bottom depth: 20 15 10 meter; zeta=-0.5 meter ')

