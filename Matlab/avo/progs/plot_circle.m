function hcirc = plotcircle(p,c,n,linw,lstyle,col)
% function hcirc = plotcircle(p,c,n,linw,lstyle,col)
% hcirc: its handle
% plot circle around the origin in 3D
% hcirc: its handle
% p: a 3D row vector, the direction of which specifies the normal of the
% plane of circle, its norm defines the radius of the circle
% c a row vector specifying its center
% n: the number of points defining the circle
% set to arbitrary position by shift of of xdat ydat zdat of hcirc 

if nargin<3
    n = 50;linw=1;lstyle='-';col='b';
end


normp=norm(p);
rho=normp;
p = p/normp;
theta=acos(p(3))/pi;

phi=atan2(p(2),p(1))/pi;
gamma=0;

alpha = linspace(0,2*pi,n)';

DATA = rho*[cos(alpha) sin(alpha) zeros(n,1)];
DATA=rotash(DATA,[phi theta gamma], c);




hcirc = plot3(DATA(:,1),DATA(:,2),DATA(:,3),'b');
set(hcirc,'LineWidth',linw,'color',col,'linestyle',lstyle);
