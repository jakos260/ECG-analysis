% make_disk
% program [VER,ITRI]=make_disk(r,nr,nphi,es);
% generates a trianglulated version of a disk
% generated in the x;y plane

% if r is scalar: a disk is made with radius r; with
% nr total number of
% (non-zero) radii 

% if r is a vector: elements of r are used as subsequent radii, nr is replaced by max(size(r)).
% If r(1) is non_zero: an annulus is produced

% total number of individual phi values involved: nphi

% es=1 creates extra shift in phi values every other anulus
% see also make_annulus
% A. van Oosterom; 20140216

function [VER,ITRI]=make_disk(r,nr,nphi,es);

dimr=max(size(r));

if nargin==3, es=0; end
if es~=0, es=1; end
delphi=2*pi/nphi;

% create VER data
nver=0;

if dimr==1,
    a=r(1);
    rj=a*(1-((1:nr+1)-1)/nr);
else,
    rj=r;
    rj=sort(rj,'descend');
end

%pause


for j=1:nr-1,
     eshift=es*rem(j+2,2)/2;
     for i=1:nphi,
            phi=(i-1+eshift)*delphi;
            nver=nver+1;
            VER(nver,1:3)=[rj(j)*sin(phi) rj(j)*cos(phi) 0];
      end 
end

if min(rj)==0,
   nver=nver+1;
   VER(nver,1:3)=[0 0 0];
end

ntri=0;
    
for j=1:nr-2;
    js=es*rem(j+1,2);
    jr=(j-1)*nphi;
    for i=1:nphi,
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ jr+icyc(i-js,nphi) jr+icyc(i+1-js,nphi) jr+nphi+icyc(i,nphi) ];
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ jr+icyc(i+1-js,nphi) jr+nphi+icyc(i+1,nphi) jr+nphi+icyc(i,nphi)];
    end
end

if min(rj)==0,
   for i=1:nphi-1,
       ntri=ntri+1;
       ITRI(ntri,1:3)=[ nver-nphi+i-1  nver-nphi+i nver];
   end
   ntri=ntri+1;
   ITRI(ntri,1:3)=[ nver-nphi+i nver-nphi nver];
   
end



