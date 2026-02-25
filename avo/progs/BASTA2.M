% basta2.m
clear
clf
TAB=loadasci('constitu.sor');
dim=size(TAB);
x=TAB(:,3);
y=TAB(:,4);
npp=dim(1);
z=y./(x.^3);

% sex specific means
nm=0;
nf=0;
for i=1:npp
if TAB(i,1)==1,
nf=nf+1;
zf(nf,1)=z(i,1);
%xf(nf,1)=x(i,1);
%yf(nf,1)=y(i,1);
else
nm=nm+1;
zm(nm,1)=z(i,1);
end
end

[mean(zf) std(zf) mean(zm) std(zm)]

%bins=1.9:.005:2.1;
hist(zf)
h=get(gca,'children');
set(h(2),'color',[1 0 0]);
hold on
hist(zm)
h=get(gca,'children');
set(h(2),'color',[0 0 1]);
