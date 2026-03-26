% make_annulus.m
% function [VER,ITRI]=make_annulus(a,b,nr,nphi,es);

% generate the trianglulated version of an annulus
% generated in the x;y plane
% a: inner radius b; outer radius;
% nr: total number of involved radii nphi: total number of individual phi values
% es=1 creates extra shift in phi values every other anulus
% see also: make_disk, or make_diabolo if aiming for (near)-equal edge lengths

function [VER,ITRI]=make_annulus(a,b,nr,nphi,es);
if nargin==4, es=0; end
delphi=2*pi/nphi;

% create VER data
nver=0;
for j=1:nr;
    r=b-(b-a)*(j-1)/(nr-1);
    eshift=es*rem(j-1,2)/2;
    for i=1:nphi,
        phi=(i-1+eshift)*delphi;
        nver=nver+1;
        VER(nver,1:3)=[r*sin(phi) r*cos(phi) 0];
    end
end

% create triangle specs
ntri=0;
for j=1:nr-1;
    for i=1:nphi-1,
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ (j-1)*nphi+i  (j-1)*nphi+i+1 j*nphi+i];
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ (j-1)*nphi+i+1 j*nphi+i+1 j*nphi+i];
    end
    ntri=ntri+1;
    ITRI(ntri,1:3)=[ j*nphi (j-1)*nphi+1  (j+1)*nphi];
    ntri=ntri+1;
    ITRI(ntri,1:3)=[ (j-1)*nphi+1 j*nphi+1 (j+1)*nphi];
end

    

