% curvature.m
% curv=curvature(y,x,ds)
% for nargin==2: curv=-(x''*y'-x'*y'')/((x')^2+(y')^2)^(3/2); y=y(s), x=x(s); 
% for nargin==1, curv=y''/(1+(y')^2)^(3/2); this assumes equal sampling intervals: ds
% 20050802; A. van Oosterom

function curv=curvature(y,ds,x)
yp=diffrows(y)/ds;
ypp=diffrows(yp)/ds;

if nargin<3,
    
curv=ypp./(1+yp.^2).^(3/2);
else,
    xp=diffrows(x)/ds;
    xpp=diffrows(xp)/ds;
    curv=-(xpp.*yp-ypp.*xp)./(xp.^2+yp.^2).^(3/2);
end

    
