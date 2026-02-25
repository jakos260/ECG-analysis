% rfgauss_over_sqrt.m
function [y,G]=rfgauss_over_sqrt(x,p,mode);
% models negative derivative of the TMP as:
% function:
% y=p(1)+p(2)*exp(-p(3)*x.^2))./(sqrt(1+p(4)*x.^2);


functype=12;

% A. van Oosterom; 20070606

[nx ndum]=size(x);
if nx < ndum,
    x=x';
    nx=ndum;
    'x changed to the required column vector format'
end
[np ndum]=size(p);
np=max(np,ndum);

pw=0.5;
exp1=exp(-p(3)*x.^2);
denom=(1+p(4)*x.^2);


y=p(1)+p(2)*exp1./denom.^pw;

if mode==0, return, end

% computation of G: the derivatives of the function y with respect to the
% parameter vector p
% used, e.g., while estimating p 
G=zeros(nx,np);
G(:,1)= ones(nx,1);
G(:,2)= exp1./denom.^pw;

G(:,3)=-p(2)*x.^2.*exp1./denom.^pw;

G(:,4)=-p(2)*pw*exp1.*x.^2./denom.^(pw+1);
