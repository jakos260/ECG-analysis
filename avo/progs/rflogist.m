

% rflogist.m
% function [y,G]=rflogist(x,p,mode);
% function: constant + logist_up)
% y=p(1)+p(2)/(1+exp(p(3)*(x-p(4)))
% p(1)  an overall shift
% p(2)  overall scaling factor 
% p(3)  determines the slope 
% p(4)  determines the timing of the slope

% funtype=11;


% A. van Oosterom; 20070613

function [y,G]=rflogist(x,p,mode);

[nx ndum]=size(x);
[np ndum]=size(p);
np=max(np,ndum);

exp1=exp(p(3)*(x-p(4)));


logist=1./(1+exp1);
loglogist=logist.*logist;
y=p(1)+p(2)*logist;

if mode==0, return, end

% computation of G: the derivatives of the function y with respect to the
% parameter vector p
% used, e.g., while estimating p 

G=zeros(nx,np);

G(:,1)=1;

G(:,2)=logist;

G(:,3)= -p(2)*(x-p(4)).* exp1.*loglogist;

G(:,4)=       p(2)*p(3)* exp1.*loglogist;

