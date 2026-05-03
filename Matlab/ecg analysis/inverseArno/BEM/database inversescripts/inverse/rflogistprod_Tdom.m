% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

% rflogistprod.m
% function [y,G]=rflogistprod(x,p,mode);
% function: 
% y=p(1)*(phidown+phidown+phinotch)
% phiup=   (p(2) * 1/(1+exp(-p(3)*(x-p(4))))
% phinotch= 1-p(5)*1/(1+exp(-p(6)*(x-p(7)))*1/(1+exp(p(8)*(x-p(7)) 
% phidown= (p(2) + 1/(1+exp(p(10)*(x-p(12)))*1/(1+exp(p(11)*(x-p(12))

% p(1)		% overall scaling factor 
% p(2)		% the initial value of phiup, phidown,phinotch
% p(3)      % determines the (positive) slope leading up to the apex
% p(4)      % determines the (negative) slope following the apex
% p(5)		% the timing of the replolarisation phase

% P.M. van Dam 07-10-08

function [y,G]=rflogistprod_Tdom(x,p,mode)

%%
% p(1)= 1.1276;	% overall scaling factor 
% p(2)= 0.84;	% the initial value of phidown 
exp1=exp( p(3)*(x-p(5)));	exp1(isinf(exp1)==1)=sign(exp1(isinf(exp1)==1))*realmax/2;
exp2=exp( p(4)*(x-p(5)));	exp2(isinf(exp2)==1)=sign(exp2(isinf(exp2)==1))*realmax/2;

%up
pup(1)= 1.;	% the initial value of phiup
pup(2)= -3.46;	% determines the slope of the upstroke 
pup(3)= 10;
exp3=exp( pup(2)*(x-pup(3)));
logist3=1./(1+exp3);
phiup=pup(1)*logist3;


logist1=1./(1+exp1);
logist2=1./(1+exp2);

phidown=(p(2)*logist1).*logist2;
phiplateau=1;
y=p(1)*phidown.*phiplateau.*phiup;
% figure(1);clf;plot(p(1)*phidown,'r'); hold on; plot(phiup,'g'); plot(phiplateau,'m');plot(y,'k')

if mode==0, return, end
%%
% computation of G: the derivatives of the function y with respect to the
% parameter vector p used, e.g., while estimating p 

% down
G(:,1)=y/p(1).*phiplateau.*phiup;
G(:,2)= p(1).*logist2.*phiplateau.*phiup;

% phidown= (p(2) + 1/(1+exp(p(3)*(x-p(5)))*1/(1+exp(p(4)*(x-p(5))
G(:,3)= -p(1)*(x-p(5)).*(logist1.^2)  .*exp1.* logist2.*phiplateau.*phiup;
G(:,4)= -p(1)*(p(2)+logist1).*(x-p(5)).*exp2.*(logist2.^2).*phiplateau.*phiup;
G(:,5)=  p(1)*(p(3)*exp1.*(logist1.^2).*logist2+p(4)*(p(2)+logist1).*exp2.*(logist2.^2) ).*phiplateau.*phiup; 

% G(:,6)= -p(1)*logist3								.*phidown.*phiup;
% G(:,7)= -p(1)*plp(1).*(x-plp(3)).*exp3.*(logist3.^2).*phidown.*phiup;


for i=1:5
	a=G(:,i);
	a(isnan(a)==1)=0;
	G(:,i)=a;
end
