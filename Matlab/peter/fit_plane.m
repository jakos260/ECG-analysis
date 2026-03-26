% fit_sphere.m
% function [radius,orig]=fit_sphere(PNTS);
% least squares fit in 3D of a sphere, applied to points (PNTS; size(nn,3)
% for nn=2: vertices at the poles; for nn=3: great circle around the vertices
% uses the general expression for a sphere in 3D:
% x^2+y^2+z^2 + ax + by + cz + d =0;; linear in [a b c d]
% origin: -[a b c]/2;   radius: sqrt((a^2+b^2+c^2)/4-d)
% A. van Oosterom
% date: 0090318


function ABCD=fit_plane(VER)

PNTS=[];
used=zeros(size(VER,1),1);
P=ones(size(VER,1),1);
next=1;
delta =min(diff(range(VER)))/1.4;

while sum(used) < length(used)
    dist=norm3d(VER-P*VER(next,:));
    dist(used==1)=10000;
    a=find(dist< delta );
    used(a)=1;
    PNTS=[PNTS;mean(VER(a,:),1)];
    
    dist(used==1)=10000;
    next = find(dist==min(dist));
   
    
end


PNTS;

Y= [PNTS ones(size(PNTS,1),1)];
[U,S,V]=svd(Y);

ABCD=V(:,4);
ABCD = ABCD ./ norm3d(ABCD(1:3));
