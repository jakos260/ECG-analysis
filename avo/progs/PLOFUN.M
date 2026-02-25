echo off;
clear
load a; 
load b;
%hold on
[X,Y]=meshgrid(-4:.4:1., -8:.2:-1);
z=1.
Z=    (a(1,1)*X+a(1,2)*Y+a(1,3)*z-b(1)).^2;
Z= Z+ (a(2,1)*X+a(2,2)*Y+a(2,3)*z-b(2)).^2;
Z= Z+ (a(3,1)*X+a(3,2)*Y+a(3,3)*z-b(3)).^2;
%mesh(X,Y,Z);
%steps=10: 10: 80;
steps=linspace(10,70,21);
%contour(X,Y,Z,steps);
%contour(X,Y,Z,40);
cont=contour(X,Y,Z,steps);
clabel(cont,'manual')
colorbar
xlabel('x-axis')
ylabel('y-axis')
t=sprintf('z=%0.5g',z)
text(2.0,-0.70,t)
%print plotje.ps
orient tall;





