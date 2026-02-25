% legpol.m
% function p=legpol(m,n,x)
% legendre polynomial of degree n and order m
% implemented for m=0 or m=1 only

% METHOD:
% (n-m)*P_n(x)=(2*n-1)*x*P_(n-1)(x)-(n-1+m)*P_(n-2)(x) (ZIE ABRAMOVITCH 8.5.3)
% implemented only for m=0 and  m=1
% for higher order: use standard matlab function: legendre

function  pn=legpol(m,n,x)
    nx=length(x);
	x=min(x,ones(1,nx));
	x=max(x,-1*ones(1,nx));
    n=max(n,0);
	m=max(m,0); m=min(m,1);
if m==0,  
    pn=ones(1,nx);
    pnm1=pn;
    pnm2=zeros(1,nx);           
    k=1; 
else,
    if n==0, pn=zeros(1,nx);k=1; end
    if n>0, pn=sqrt(1-x.^2); pnm1=pn; pnm2=zeros(1,nx); end
    k=2; 
end
while k<=n,
    pn=((2*k-1)*x.*pnm1-(k+m-1)*pnm2)/(k-m);
    pnm2=pnm1;
	pnm1=pn;
    k=k+1;
end   



