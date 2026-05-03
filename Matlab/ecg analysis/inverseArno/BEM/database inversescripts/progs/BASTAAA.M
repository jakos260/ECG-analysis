% basta.m
clear
clf
TAB=loadasci('constitu.sor');
dim=size(TAB);
xx=TAB(:,3);
yy=TAB(:,4);
%y=sqrt(yy./(xx.^3));
%x=yy.^(1/3)./xx;
x=xx;
y=yy;
npp=dim(1);

% sex specific means
nm=0;
nf=0;
for i=1:npp
if TAB(i,1)==1,
nf=nf+1;
xf(nf,1)=x(i,1);
yf(nf,1)=y(i,1);
else
nm=nm+1;
xm(nm,1)=x(i,1);
ym(nm,1)=y(i,1);
end
end
[mean(xf) mean(xm)]
[mean(yf) mean(ym)]
corrcoef(x,y)

plot(x,y,'+')
scal=axis;
scal(1)=0;
scal(3)=0;
hold on
for i=1:npp,
plot(x(i),y(i),'+')
if TAB(i,1)== 1,
h=get(gca,'children');
set(h(1),'color',[1 0 0]);
else plot(x(i),y(i),'blue+')
end
end
hold on
