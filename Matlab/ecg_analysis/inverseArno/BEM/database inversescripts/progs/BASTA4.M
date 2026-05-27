% basta4.m
clear
clf
TAB=loadmat('basicsta.lst');
dim=size(TAB);
npp=dim(1);
nvar=dim(2);

%select=[3 10 4 5 18 19 22 24 25 26 20 21 17 11 12 14 13 8 7 6 9 16 15];
%DATA=TAB(:,select);
DATA=TAB;

save=0;
if save ~=0,
STAT=zeros(npp,6);
for k=1:nvar,
nm=0;
nf=0;

for i=1:npp,
dat=DATA(i,k);
if DATA(i,2)==0,
nf=nf+1;
yf(nf,1)=dat;
else
nm=nm+1;
ym(nm,1)=dat;
end
end
[t df]=ttest2(yf(:,1),ym(:,1));
STAT(k,:)=[mean(yf) std(yf) mean(ym) std(ym) mean(ym)/mean(yf) t];
end
STAT
saveasci('temp.lst',STAT)
end

%ANGS=loadasci('angles.lst');

x=DATA(:,11);
y=DATA(:,16);

plot(x,y,'+')
size=axis;
size(1)=0;
size(3)=0;
axis(size)
pietsie=(size(2)-size(1))/100
hold on
for i=1:npp,
if DATA(i,2)==0, 
plot(x(i),y(i),'r+')
else
plot(x(i),y(i),'b+')
end
tekst=sprintf('%i',DATA(i,1));
text(x(i)+pietsie,y(i),tekst)
end

statall=[mean(x) std(x) mean(y) std(y)]
rhoall=corrcoef(x,y)
nm=0;
nf=0;
for i=1:npp
if DATA(i,2)==0,
nf=nf+1;
xf(nf,1)=x(i);
yf(nf,1)=y(i);
else
nm=nm+1;
xm(nm,1)=x(i);
ym(nm,1)=y(i);
end
end
'x-stat'
statfemale=[mean(xf) std(xf)];
statmale=[mean(xm) std(xm)];
[t df]=ttest2(xf(:,1),xm(:,1));
[statfemale statmale mean(xm)/mean(xf) t]
'y-stat'
statfemale=[mean(yf) std(yf)];
statmale=[mean(ym) std(ym)];
[t df]=ttest2(yf(:,1),ym(:,1));
[statfemale statmale mean(ym)/mean(yf) t]

af=xf'*yf/(xf'*xf);
sdaf=sqrt((yf-af*xf)'*(yf-af*xf)/(npp-1))/(xf'*xf);
staf=[af sdaf]
am=xm'*ym/(xm'*xm);
sdam=sqrt((ym-am*xm)'*(ym-am*xm)/(npp-1))/(xm'*xm);
stam=[am sdam]

M(:,1)=x.*(1-DATA(:,2));
M(:,2)=x.*DATA(:,2);
%M(:,3)=ones(npp,1);
a=pinv(M)*y
xx=M*a;
'regr after allowing for sex-adapted scaling: rho(zerointercept)'
corrcoef(xx,y)
aall=x'*y/(x'*x);
sda=sqrt((y-aall*x)'*(y-aall*x)/(npp-1))/(x'*x);
stataall=[aall sda]

save=0;
if save~=0,
  hold off
  plot(x,y,'+')
  scal=axis;
  scal(1)=0;
  pietsie=(scal(2)-scal(1))/100
  scal(3)=0;
  hold on

  for i=1:npp,
   if DATA(i,2)== 0,
    plot(x(i),y(i),'red+')
   else plot(x(i),y(i),'blue+')
   end
  tekst=sprintf('%i',DATA(i,1));
  text(x(i)+pietsi,y(i),tekst)
  end
end

%plot([0 scal(2)],[0 aall*scal(2)],'w:')
%plot([0 scal(2)],[0 a(1)*scal(2)],'r:')
%plot([0 scal(2)],[0 a(2)*scal(2)],'b:')

save=0;
if save ~=0, 
  COR=corrcoef(TAB);
  saveasci('unisex.cor',COR)
  SUMMARY=zeros(nvar,4); 
   for l=1:nvar,
   SUMMARY(l,:)=[mean(TAB(:,l)) std(TAB(:,l)) min(TAB(:,l)) max(TAB(:,l))];
   end
  saveasci('unisumar.lst',SUMMARY)
end

sav=0;
if save~=0,
  corextr=zeros(nvar,7);
  for i=1:nvar,
   extrpos=0;
   extrneg=0;
   small=1;
     for j=1:nvar,
     if i~=j,
     if abs(COR(i,j)) < small,
     small=abs(COR(i,j));
     is=j;
     end
  if COR(i,j)> extrpos,
  extrpos=COR(i,j);
   ip=j;
  else
  if COR(i,j)<extrneg,
   extrneg=COR(i,j);
   in=j;
   end
   end
   end
end
corextr(i,1)=i;
corextr(i,2)=extrneg;
corextr(i,3)=in;
corextr(i,4)=small;
corextr(i,5)=is;
corextr(i,6)=extrpos;
corextr(i,7)=ip;
end
corextr
%saveasci('tmp.lst',corextr)
%saveasci('unisexn.cor',COR)
end
end

