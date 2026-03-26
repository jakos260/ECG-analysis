% makeglobe.m
% function [VER, ITRI]=makeglobe(ntheta,nphi);
% generates triangular mesh for sphere using equidistant theta phi values
% 2003-06-31; A. van Oosterom

function [VER, ITRI]=makeglobe(ntheta,nphi);
ndat=ntheta*nphi;
r=1;
theta=1:ntheta;
phi=1:nphi;
theta=(theta-1)*pi/(ntheta-1);
phi=(phi-1)*(2*pi)/nphi;
X=r*sin(theta)'*cos(phi);
Y=r*sin(theta)'*sin(phi);
Z=r*cos(theta)'*ones(1,nphi);
X=X';Y=Y';Z=Z';
x=X(:);
y=Y(:);
z=Z(:);
x([2:nphi ndat-nphi+2:ndat])=[];
y([2:nphi ndat-nphi+2:ndat])=[];
z([2:nphi ndat-nphi+2:ndat])=[];
VER=[x y z];
nver=length(VER);

% create arctic cap; north pole
ITRI(1:nphi,1)=1;
ITRI(1:nphi,2)=[3:nphi+2]';
ITRI(nphi,2)=2;
ITRI(1:nphi,3)=[2:nphi+1]';
jb=2;

% create middle part
for i=1:ntheta-3,
   for j=jb:jb+nphi-1;
       ITRI=[ITRI;[j j+nphi+1 j+nphi]];
       ITRI=[ITRI;[j j+1 j+nphi+1 ]];
   end
   ntri=length(ITRI);
   ITRI(ntri-1,2)=jb+nphi;
   ITRI(ntri,2:3)=[jb jb+nphi];
   jb=jb+nphi;
end

% create antarctic cap; south pole
   for j=jb:jb+nphi-1;
       ITRI=[ITRI;[j j+1 nver]];
   end
   
   ntri=length(ITRI);
   ITRI(ntri,2)=jb;





