% INTDIF.m
% study numerical variants of
% differentiation and integration
n=31;
DIF=zeros(n);
for i=2:n-1,
DIF(i-1,i)=-.5;
DIF(i+1,i)=.5;
end
DIF(1,1)=-1;
DIF(n-1,n)=-1;
DIF(2,1)=1;
DIF(n,n)=1;

INT=zeros(n);
for j=1:n,
% INT(1,j)=0.5;
%for i=1:j,

for i=2:j,
INT(i,j)=1;
end
%INT(j,j)=0.5;
end

DIFinv=pinv(DIF);
INTinv=pinv(INT);
sig=zeros(1,n);
ib=round(n/3);
ie=round(2*n/3);
for i=ib:ie,
sig(i)=2;
%sig(i)=i-ib;
end
for i=1:ib,
sig(i)=1;
end

% eerst diff dan int
clf
plot(sig)
hold on
plot(sig*DIF,'r') 
plot(sig*DIF*INT,'b')
pause
clf
plot(sig)
hold on
plot(sig*DIF,'r') 
plot(sig*DIF*DIFinv,'b')

pause
clf
plot(sig)
hold on
plot(sig*INTinv,'r') 
plot(sig*INTinv*INT,'b')

pause
clf
plot(sig)
hold on
plot(sig*INTinv,'r') 
plot(sig*INTinv*DIFinv,'b')


pause
% eerst int dan diff
clf
plot(sig)
hold on
plot(sig*INT,'r') 
plot(sig*INT*DIF,'b')

pause
clf
plot(sig)
hold on
plot(sig*INT,'r') 
plot(sig*INT*INTinv,'b')

pause
clf
plot(sig)
hold on
plot(sig*DIFinv,'r') 
plot(sig*DIFinv*DIF,'b')

pause
clf
plot(sig)
hold on
plot(sig*DIFinv,'r') 
plot(sig*DIFinv*INTinv,'b')

pause
clf

LAP=zeros(n);
for i=2:n-1,
LAP(i-1,i)=1;
LAP(i,i)=-2;
LAP(i+1,i)=1;
end
LAP(1,1)=-1;
LAP(n-1,n)=-1;
LAP(2,1)=1;
LAP(n,n)=1;

e1=zeros(n,1);
e1(1)=10;

A1=LAP';
B1=zeros(n);

C1=[DIF e1];
D1=[eye(n) e1];

A=A1'*A1;
YA=A1'*B1;

B=C1*C1';
YB=D1*C1';
Y=YA+YB;
[U,SA,VA]=svd(A);
[UB,SB,V]=svd(B);
Ytil=U'*Y*V;
Xtil=zeros(n);
for i=1:n,
for j=1:n,
Xtil(i,j)=Ytil(i,j)/(SA(i,i)+SB(j,j));
end
end
X=U*Xtil*V';

clf
plot(sig)
hold on
plot(sig*X,'r')
plot(sig*X*DIF,'b')

