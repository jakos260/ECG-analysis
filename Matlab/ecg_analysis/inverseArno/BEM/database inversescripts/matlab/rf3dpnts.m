% rf3dpnts.m
% function(y,G)=rf3dpnts(X,par,mode)
% called by paramest to fit np points X(1:np,1:3) in 3D to points (Y(1:np,1:3) 
% by minimizing their 3D (squared distance) by means by optimizing rotation followed by shift
% par(1:3): the rotation angles indegrees; par(4:6), the shift
% y=Y(:); with size=[np*3, 1)
% mode==0 computes function values only; else: gradients are also computed
% G gradient matrix with size=[np*3 6]

function [y,G]=rf3dpnts(X,par,mode) 
npar=length(par);
% 20041025; A. van Oosterom
	
% compute function
[np jdum]=size(X);;
if npar>3,
  Y=rotash(X,par(1:3),par(4:6));
else,
 Y=rotash(X,par(1:3),[0 0 0]);
end   
y=Y(:);
np3=length(y);
if mode==0, return, end
    
% compute gradients
G=zeros(3*np,npar);
Y=rotash(X,[par(1)+0.5 par(2:3)],zeros(1,3));
G(1:np3,1)=Y(:)*pi/180;
Y=rotash(X,[par(1) par(2)+0.5 par(3)], zeros(1,3));
G(:,2)=Y(:)*pi/180;
Y=rotash(X,[par(1:2) par(3)+0.5], zeros(1,3));
G(:,3)=Y(:)*pi/180;

if npar<=3, return,end
G(:,4)=[ones(np,1);zeros(2*np,1)];
G(:,5)=[zeros(np,1);ones(np,1);zeros(np,1);];
G(:,6)=[zeros(2*np,1);ones(np,1)];
