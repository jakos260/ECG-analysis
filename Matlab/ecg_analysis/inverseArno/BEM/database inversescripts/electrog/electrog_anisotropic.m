% electrog_anisotropic.m
% A. van Oosterom; 20080725

% all trianglulated version of an annulus
% generated in the x;y plane
% inner radius a; outer radius b;
% assign anisotropic velocities: vr (radieel, transmural) and vl longitudinal
% along its edges,
% initiate activation at any node ( in particular either endo- or epicarial stimulus site)
% plot the resulting isochrones

clear all

cmap='tims.mcm';

a =0.6
b =1.
nr=20
delr=(b-a)/nr;
nl=round(2*pi*a/delr)

[VER,ITRI]=make_annulus(a,b,nr,nl,0);
VER=[VER(:,3) VER(:,2) VER(:,1)];

next=0;
if next==1,
   figure(1)
   clf
   VALS=VER;
   %zebra=10;
   triplot_contour
end

nver=size(VER,1)

vl=0.5; vr=0.5;

delr=delr/vr;

% setup adjacency matrix
A=zeros(nver,nver);

n=1;
% set distances/(travel times) along edges
for i=1:nr,
    nbeg=(i-1)*nl+1;
    dell=norm3d(VER(nbeg,:)-VER(nbeg+1,:))/vl;
    deloblique=sqrt(delr^2+dell^2);
    for j=1:nl,
         if j==nl, 
             k=n+1-nl;
         else,
             k=n+1;
         end
         A(n,k)=dell;
         if i < nr,
             A(n,n+nl)=delr;
             kp=n+nl+1; if j==nl, kp=kp-nl; end
             A(n,kp)=deloblique;
             km=n+nl-1; if j==1, km=km+nl; end
             A(n,km)=deloblique;             
         end
         n=n+1;
    end
end

A=A+A';

DIST=graphdist(A);

% interesting nodes:
%                    1                          (epi)   say, anterior wall 
%                    1+(nr-1)*nl                (endo)  say, anterior wall
%                    1+round(nl/2)              (epi)   say, posterior wall
%                    1+round(nl/2+(nr-1)*nl)    (endo)   say, posterior wall

innos=[1  1+(nr-1)*nl 1+round(nl/2) 1+round(nl/2+(nr-1)*nl)]  % 

next=1;
if next==1,
   figure(1)
   clf
   zebra=-0.2;
   VALS=DIST;
   triplot_contour
end   





