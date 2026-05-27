% basta.m
clear
clf
TAB=loadasci('constitu.txt');
plot(TAB(:,3),TAB(:,4),'+')
dim=size(TAB);
npp=dim(1)
%M=zeros(npp,2);
x=TAB(:,3);
y=TAB(:,4);
%M(:,1)=x.*TAB(:,2);
%M(:,2)=x.*(1-TAB(:,2));
%M(:,3)=TAB(:,12);
%a=pinv(M)*y
%x=M*a;
corrcoeff(x,y)
plot(x,y,'+')
scal=axis;
scal(1)=0;
scal(3)=0;
hold on
for i=1:npp,
if TAB(i,1)== 1,
plot(x(i),y(i),'red+')
else plot(x(i),y(i),'blue+')
end
end
%COR=corrcoef(TAB)
