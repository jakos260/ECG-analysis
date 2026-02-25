% rflogistprod_new.m
% function [y,G]=rflogistprod(x,p,mode);
% models negative derivative of the TMP as:
% function: constant*(constant + logist_up)*logist_down:
% y=(p(1)*(p(2)+1/(1+exp(p(3)*(x-p(5)))*1/(1+exp(p(4)*(x-p(5)))
% p(1)  an overall scaling factor of Tdom
% p(2)  the initial value of the Tdom like function; = approx: y(0)
% -p(3) determines the (positive) slope leading up to the apex
% -p(4) determines the (negative) slope following the apex
% p(5)  the timing of the apex of Tdom)
% TMP, with normalized upstroke, may be computed as 1-cumsum(y)/sum(y)
% if mode==0, just the function is computed, else
% also the derivatives of the function with respect of its parameters 

% funtype=6;

% A. van Oosterom; 20050606

function [y,G]=rflogistprod(x,p,mode);

[nx ndum]=size(x);
[np ndum]=size(p);

exp1=exp(p(3)*(x-p(5)));
exp2=exp(p(4)*(x-p(5)));

logist1=1./(1+exp1);
logist2=1./(1+exp2);

y=p(1)*(p(2)+logist1).*logist2;

if mode==0, return, end

% computation of G: the derivatives of the function y with respect to the
% parameter vector p
% used, e.g., while estimating p 

G(:,1)=y/p(1);

G(:,2)=p(1)*logist2;

G(:,3)= -p(1)*(x-p(5)).*(logist1.^2).*exp1.*logist2;

G(:,4)= -p(1)*(p(2)+logist1).*(x-p(5)).*exp2.*(logist2.^2);

G(:,5)=p(1)*( p(3)*exp1.*(logist1.^2).*logist2+p(4)*(p(2)+logist1).*exp2.*(logist2.^2) ); 
