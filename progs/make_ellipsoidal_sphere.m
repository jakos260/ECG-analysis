% make_ellipsoidal_sphere.m
% function [VER, ITRI]=make_ellipsoide(nz,nphi,axis_ratio);
% generate triangular mesh for prolate ellipsoidal sphere using equidistant z and phi values
% 2012-06-31; A. van Oosterom
% rho^2+z^2=1
clear all
nz=10;
nphi=8;
axis_ratio=0.5; % estimated from Figure 1 ; Schilling et al Circulation 1998;98:887-898.

b=1.8;
a=b*axis_ratio;
z=-b:2*b/(nz-1):b;
rho=a*sqrt(1-z.^2/b^2);
phi=1:nphi;
phi=(phi-1)*(2*pi)/nphi;

% closing cap
VER=[0 0 z(1)];
ADDVERS=[rho(2)*[cos(phi)' sin(phi)'] z(2)*ones(nphi,1)];
VER=[VER;ADDVERS]; 
lista=1; listb=lista+(1:nphi);
ITRI=make_peel(VER,lista,listb,1,2);
lista=1; listb=lista+(1:nphi);

for i=2:8;
    ADDVERS=[rho(i+1)*[cos(phi)' sin(phi)'] z(i+1)*ones(nphi,1)];
    VER=[VER;ADDVERS];
    lista=listb; listb=lista(end)+(1:nphi);
    ADD_TRIS=make_peel(VER,lista,listb,1);
    ITRI=[ITRI;ADD_TRIS];
end

VER=[VER; [0 0 b]];
lista=listb; listb=lista(end)+1;
ADD_TRIS=make_peel(VER,lista,listb,1,2);
ITRI=[ITRI;ADD_TRIS];
ITRI=ITRI(:,[2 1 3]);
    
figure(1)
clf
VALS=VER;
zebra=-20;

%VER=rotash(VER,[0.5 45/180*pi -0.5],zeros(1,3));

triplot_contour
axis equal
hold on

pause

% relabel vertices
oldlabs=[66 58:-8:2 59:-8:3  60:-8:4  61:-8:5 62:-8:6 63:-8:7 64:-8:8 65:-8:9   1];

newlabs=1:66;

[VER,ITRI]=tri_relabel(VER,ITRI,oldlabs, newlabs);
    
figure(1)
clf
VALS=VER;
zebra=-0.4;

%VER=rotash(VER,[0.5 45/180*pi -0.5],zeros(1,3));

triplot_contour
axis equal
hold on


for i=1:65,
    text(VER(i+1,1),VER(i+1,2),VER(i+1,3),num2str(i))
end

'taking cm as a unit:'
'balloon volume [ml]'
trivolume(VER,ITRI) 

axis_ratio

'distance between poles  [cm]'
norm(VER(1,:)-VER(66,:))

    
