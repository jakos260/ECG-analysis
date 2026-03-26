% rflogist_prods.m
% function [y,G]=rflogist_prods(x,p,mode);
% models negative derivative of the TMP as:
% function: (constant + logist_up)*logist_down:
% y=p(1)/(p(2)+exp(p(3)*(x-p(5)))*1/(1+exp(p(4)*(x-p(6)))
% x must be a column vector; p a row vector
% p(1)  overall scaling 
% p(2)  ST-like early shift 
% p(3)  determines the slope of the first logist
% p(4)  determines the slope of the second logst
% p(5)  the timing of the first logist
% p(6)  the timing of the second logist

% funtype=13


% A. van Oosterom; 20070606

function [y,G]=rflogist_prods(x,p,mode);

[nx ndum]=size(x);
if nx < ndum,
    x=x';
    nx=ndum;
end
[np ndum]=size(p);
np=max(np,ndum);

exp1=exp(p(3)*(x-p(5)));
exp2=exp(p(4)*(x-p(6)));

logist1=1./(p(2)+exp1);
logist2=1./(1+exp2);
log1log2=logist1.*logist2;

y=p(1)*log1log2;

if mode==0, return, end

% computation of G: the derivatives of the function y with respect to the
% parameter vector p
% used, e.g., while estimating p 
G=zeros(nx,np);

G(:,1)=y/p(1);

G(:,2)=-y.*logist1;

G(:,3)= -y.*logist1.*exp1.*(x-p(5));

G(:,4)= -y.*logist2.*exp2.*(x-p(6));

G(:,5)= p(3)* y.*logist1.*exp1; 

G(:,6)= p(4)* y.*logist2.*exp2;