% curvature.m
% curv=curvature(y,ds,x)
% for nargin==2, curv=y''/(1+(y')^2)^(3/2); y=y(s)
% for nargin==3: curv=-(x''*y'-x'*y'')/((x')^2+(y')^2)^(3/2); y=y(s), x=x(s);

% assumes equal sampling intervals ,ds, along the curves!!
% 20130908; A. van Oosterom

% for non_equidistant samples along a curve in 3D use  
% fit_sphere(POS), applied to subsequent M>=4 of observations
% POS(i-1:i+1),:); POS=[x y zeros(nsamp,1)]
% for frequent use of 2D curves, the function fit_sphere(POS)  may
% be replaced by fit_circle(POS)

function curv=curvature(y,ds,x)
if size(y,1)>1, y=y'; end

yp=diffrows(y)/ds;
ypp=diffrows(yp)/ds;

if nargin<3,
    curv=ypp./(1+yp.^2).^(3/2);
else,
    if size(x,1)>1, x=x'; end
    xp=diffrows(x)/ds;
    xpp=diffrows(xp)/ds;
    curv=-(xpp.*yp-ypp.*xp)./(xp.^2+yp.^2).^(3/2);
end


