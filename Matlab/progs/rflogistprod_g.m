% rflogistprod_g.m
% function [y,G]=rflogistprod_g(x,p,mode);
% dominant T + Gauss
% models dominant T wave as:
% function: (constant + logist_up)*logist_down +
% y=(p(1)*(p(2)+1/(1+exp(p(3)*(x-p(5)))*1/(1+exp(p(4)*(x-p(5)))
%   + p(6)*exp(p(7)*(x-p(8))^2) (Uwave)
% p(1)  an overall scaling factor of Tdom
% p(2)  the initial value of the Tdom like function; = approx: y(0)
% -p(3) determines the (positive) slope leading up to the apex
% -p(4) determines the (negative) slope following the apex
% p(5)  the timing of the apex of Tdom)
% p(6)  amplitude U wave
% p(7)  width (sigma) of the U wave
% p(8)  timing of the U wave
% TMP, with normalized upstroke, may be computed as 1-cumsum(y)/sum(y)
% if mode==0, just the function is computed, else
% also the derivatives of the function with respect of its parameters 

% funtype 7

% A. van Oosterom; 20130114

function [y,G]=rflogistprod_g(x,p,mode);

nx=size(x,1);
np=size(p,1);

exp1=exp(p(3)*(x-p(5)));
exp2=exp(p(4)*(x-p(5)));

logist1=1./(1+exp1);
logist2=1./(1+exp2);
exp3=exp(-0.5*((x-p(8))/p(7)).^2);

uwave=p(6)*exp3;
yy=(p(2)+logist1).*logist2;
y=p(1)*yy+uwave;

if mode==0, return, end

% computation of G: the derivatives of the function y with respect to the
% parameter vector p
% used, e.g., while estimating p 

G(:,1)=yy/p(1);

G(:,2)=p(1)*logist2;

G(:,3)= -p(1)*(x-p(5)).*(logist1.^2).*exp1.*logist2;

G(:,4)= -p(1)*(p(2)+logist1).*(x-p(5)).*exp2.*(logist2.^2);

G(:,5)=p(1)*( p(3)*exp1.*(logist1.^2).*logist2+p(4)*(p(2)+logist1).*exp2.*(logist2.^2) ); 

G(:,6)=exp3;
G(:,7)=uwave.*(x-p(8)).^2/p(7)^3;
G(:,8)=uwave.*(x-p(8))/p(7)^2;