% make_cylinder.m
% function [VER, ITRI]=make_cylinder(nphi,nz,);
% generates triangular mesh for cylinder of unit radius r, and unit length
% disk oriented along z-axis
% using nphi equidistant phi values and and nz z values;
% simple closure at 'caps'
% 2006-12-09; A. van Oosterom

function [VER, ITRI]=make_cylinder(nphi,nz)


z=(0:nz-1)/(nz-1);
phi=(0:nphi-1)/(nphi)*2*pi;
cosphi=cos(phi);
sinphi=sin(phi);
nver=nz*nphi+2;
VER=zeros(nver,3);
VER(nver,3)=1;


jend=1;
for i=1:nz,
    jbeg=jend+1;
    jend=jend+nphi;
    VER(jbeg:jend,1:3)=[cosphi' sinphi' ones(nphi,1)*z(i)];
end
 %ITRI=zeros(2*nver-4,3);   
% create bottem cap
ITRI(1:nphi,1)=1;
ITRI(1:nphi,2)=[3:nphi+2]';
ITRI(nphi,2)=2;
ITRI(1:nphi,3)=[2:nphi+1]';
jb=2;

% create middle part
for i=1:nz-1,
   for j=jb:jb+nphi-1;
       ITRI=[ITRI;[j j+nphi+1 j+nphi]];
       ITRI=[ITRI;[j j+1 j+nphi+1 ]];
   end
   ntri=length(ITRI);
   ITRI(ntri-1,2)=jb+nphi;
   ITRI(ntri,2:3)=[jb jb+nphi];
   jb=jb+nphi;
end

% create top cap
   for j=jb:jb+nphi-1;
       ITRI=[ITRI;[j j+1 nver]];
   end
   
   ntri=length(ITRI);
   ITRI(ntri,2)=jb;
   ITRI(:,[2 3])=ITRI(:,[3 2]);




