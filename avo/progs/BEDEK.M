% file bedek
% poging twee driehoeken elkaar volledig
% late bedekken  vanuit een gezichtspunt:
% de oorsprong
% input eerst de driehoek met het kleinste oppervlak
% per driehoek: eerst de grootste zijde, daarna aflopend
echo off
clf
clear

a=[5 4 3];
b=[8 5 5];
% set up A in the plane y=0 with A(:,1)
A=[0    a(1) 
   0    0
   0    0];
A(1,3)=-(a(2)^2-a(3)^2-a(1)^2)/(2*a(1));    
A(2,3)=0;
A(3,3)=sqrt(a(3)^2-A(1,3)^2);  
% set up B in the plane y=a1(1)
B=[0        b(1)    
   a(1)     a(1)
   0         0];
B(1,3)=-(b(2)^2-b(3)^2-b(1)^2)/(2*b(1));    
B(2,3)=a(1);
B(3,3)=sqrt(b(3)^2-B(1,3)^2);  
q=[0;0;0;];
frac=.54;
phi=frac*pi/2;

xbeg=-2.;
xend=2.;
nxsteps=51;
xincr=(xend-xbeg)/(nxsteps-1);
ybeg=-2.1;
yend=2.1;
nysteps=53;
yincr=(yend-ybeg)/(nysteps-1);
opt=1.e+10;
for k=1:nxsteps,
x(k)=xbeg+(k-1)*xincr;
for l=1:nysteps,
p=B(:,1);
p(1)=p(1)+x(k);
y(k)=ybeg+(l-1)*yincr;
p(2)=p(2)+y(k);
psi=atan2(p(2),p(1));
q(1)=p(1)-b(1)*(cos(psi)*cos(phi)-sin(psi)*sin(phi));
q(2)=p(2)-b(1)*(sin(psi)*cos(phi)+cos(psi)*sin(phi));
%pause
v=p(1:2)-A(1:2,1);
w=q(1:2)-A(1:2,2);

M=[v -w];
lambda=pinv(M)*(q(1:2)-p(1:2));
r=p(1:2)+lambda(1)*v;
r(3)=0;
l1=norm(r);
l2=a(3);
l3=norm(A(:,3)-r);
cosa=(A(:,3)-r)'*(-r)/(l1*l3);
l1=norm(p-r);
l2=b(3);
mu=(l1*(cosa+sqrt((l2/l1)^2-1+cosa^2)))/norm(A(:,3)-r);

org=r;
s=org+mu*(A(:,3)-r);

dif1=s-p;
f(1)=(norm(dif1))^2;
res(1)=f(1)-b(3)^2;

dif2=s-q;
f(2)=(norm(dif2))^2;
res(2)=f(2)-b(2)^2;

sqerr=res(1)^2+res(2)^2;
if sqerr< opt,
opt=sqerr;
specs=[p q r s];
end
RESULT(l,k)=sqerr;
end
end
%plot(x,result)
[frac min(min(RESULT))]
[X Y]=meshgrid(xbeg:xincr:xend, ybeg:yincr:yend); 
%contour(X,Y,RESULT)
mesh(X,Y,RESULT)
xlabel('delx')
ylabel('dely')
specs
pause
clf
p=specs(:,1);
q=specs(:,2);
r=specs(:,3);
s=specs(:,4);
axis([-3 10 -5 15 0 3])
xlabel('xaxis')
ylabel('yaxis')
%box
line([A(1,:) A(1,1)],[A(2,:) A(2,1)],[A(3,:) A(3,1)])
set(line([p(1) q(1) s(1) p(1)],[p(2) q(2) s(2) p(2)],[p(3) q(3) s(3) p(3)]),'Color','green');
set(line([r(1) p(1)],[r(2) p(2)],[r(3) p(3)]),'Color','red');
set(line([r(1) q(1)],[r(2) q(2)],[r(3) q(3)]),'Color','red');
set(line([r(1) s(1)],[r(2) s(2)],[r(3) s(3)]),'Color','red');
set(line([r(1) A(1,1)],[r(2) A(2,1)],[r(3) A(3,1)]),'Color','red');
set(line([r(1) A(1,2)],[r(2) A(2,2)],[r(3) A(3,2)]),'Color','red');
set(line([r(1) A(1,3)],[r(2) A(2,3)],[r(3) A(3,3)]),'Color','red');
tekst=sprintf('%s','P')
text(r(1),r(2),r(3)+.1,tekst);