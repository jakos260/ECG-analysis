% fit_circle.m
% function [radius,orig]=fit_circle(PNTS);
% least squares fit in 2D of a circle, applied to points (PNTS; size(nn,2)
% for nn=2: vertices at the poles; for nn=3: great circle around the vertices
% uses the general expression for a circle in 2D:
% x^2+y^2 + ax + by  + c =0; linear in [a b c ]
% origin: -[a b ]/2;   radius: sqrt((a^2+b^2)/4-c)
% use: fit_sphere for fitting a circe in 3D
% A. van Oosterom
% date: 2013

function [radius,orig]=fit_circle(PNTS)
nn=size(PNTS,1);
meanpnts=mean(PNTS); 
% scale after shift to zero mean to maximze numerical stability
R=PNTS-ones(nn,1)*meanpnts;
scal=norm(R(:));
R=R/scal;
A=[R ones(nn,1)];
b=-sum(R.^2,2);
if rank(A)>2,
 x=pinv(A)*b;
orig=-scal*x(1:2)'/2 + meanpnts; %  rescaled; origin reshifted  
radius=scal*sqrt((x(1)^2+x(2)^2)/4-x(3));
else
    radius=inf; orig=[inf inf];
end



