% fit_sphere.m
% function [radius,orig]=fit_sphere(PNTS);
% least squares fit in 3D of a sphere, applied to points (PNTS; size(nn,3)
% for nn=2: vertices at the poles; for nn=3: great circle around the vertices
% uses the general expression for a sphere in 3D:
% x^2+y^2+z^2 + ax + by + cz + d =0; linear in [a b c d]
% origin: -[a b c]/2;   radius: sqrt((a^2+b^2+c^2)/4-d)
% A. van Oosterom
% date: 20140929

function [radius,orig]=fit_sphere(PNTS);
nn=size(PNTS,1);
meanpnts=mean(PNTS); 
% scale after shift to zero mean to maximze numerical stability
R=PNTS-ones(nn,1)*meanpnts;
scal=norm(R(:));
R=R/scal;
A=[R ones(nn,1)];
b=-sum(R.^2,2);
x=pinv(A)*b;
orig=-scal*x(1:3)'/2 + meanpnts; %  rescaled; origin reshifted  

radius=scal*sqrt((x(1)^2+x(2)^2+x(3))/3-x(4));

% res=norm3d(PNTS-ones(nver,1)*orig);
% radius=rms(res);
% least squares estimate of radius of the sphere
% produces identical result




