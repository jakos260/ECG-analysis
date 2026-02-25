% make_arrow.m
% function [VER, ITRI]=make_arrow(frhz,frhr,lz,r,nz,nphi);
% generates triangular mesh for drawing an arrow in 3D;  of length lz
% pointing upward along z-axis
% frhz (<1) is the fraction of the total length occupied by the arrowhead
% frhr (>1) is the fraction of the shaft radius specifying the maximum width
% of the arrowhead
% uses nphi equidistant phi nalue and and nz z values;
% simple closure at 'bottom'
% defaults:[0.1 3 1 0.01 3 16 ]
% 2007-02-21; A. van Oosterom


function [VER, ITRI]=make_arrow(parms);

if nargin==0, 
    
   frhz=0.1;
   frhr=3;
   lz=1;
   r=0.01;
   nz=3;
   nphi=16;
   lz=1;
   r=0.01;
   nz=3;
   nphi=16;
else,
    
   frhz=parms(1);
   frhr=parms(2);
   lz=parms(3);
   r=parms(4);
   nz=parms(5);
   nphi=parms(6);  
end


z=lz*(1-frhz)*(0:nz-1)/(nz-2);
phi=(0:nphi-1)/(nphi-1)*2*pi;
cosphi=cos(phi);
sinphi=sin(phi);
nver=nz*nphi+2;
VER=zeros(nver,3);
VER(nver,3)=lz;
jend=1;
for i=1:nz-1,
    jbeg=jend+1;
    jend=jend+nphi;
    VER(jbeg:jend,1:3)=[r*cosphi' r*sinphi' lz*ones(nphi,1)*z(i)];
end
% branch out to arrow head
jbeg=jend+1;
jend=jend+nphi;
VER(jbeg:jend,1:3)=[r*frhr*cosphi' r*frhr*sinphi' lz*ones(nphi,1)*z(i)];

% create bottom cap
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




