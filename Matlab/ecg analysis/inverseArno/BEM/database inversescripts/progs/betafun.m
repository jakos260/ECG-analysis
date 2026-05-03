% betafun.m
% function [y,ypeak,modex,meanx,sdx]=betafun(x,alpha,beta);
% beta func:y=x^alpha*(1-x)^beta; treated as a probibility density
% modex  (location peak) at x=alpha/(alpha+beta)*xmax
% meanx=(alpha+1)/(alpha+beta+2)*xmax
% sdx=sqrt((alpha+1)*(beta+1)/(alpha+beta+2)^2*(alpha+beta+3))*xmax
% Abram page 930

function [y,ypeak,modex,meanx,sdx]=betafun(x,alpha,beta);
xmax=max(x);
x=x/xmax;
y=(x.^alpha).*((1-x).^beta);
modex=alpha/(alpha+beta);
ypeak=(modex^alpha)*((1-modex)^beta);
ypeak=ypeak/sum(y);
y=y/sum(y);
modex=modex*xmax;
meanx=(alpha+1)/(alpha+beta+2)*xmax;
sdx=sqrt((alpha+1)*(beta+1)*(alpha+beta+3))/(alpha+beta+2)*xmax;

