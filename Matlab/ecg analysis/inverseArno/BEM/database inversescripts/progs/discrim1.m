% discrim1.m
% testing:
% twoclass distrimination in 2-d space
clf
clear
file1='constitu.man'
man=loadasci(file1);
GR1=man(:,3:4);
file2='constitu.vr'
vrouw=loadasci(file2);
GR2=vrouw(:,3:4);
nn=size(GR1)
n1=nn(1)
nn=size(GR2)
n2=nn(1)
ALL=[GR1;GR2];
na=size(ALL);
ma=mean(ALL)
sall=std(ALL)
ALL=(ALL-ones(na(1),1)*ma)./(ones(na(1),1)*sall);
GR1=(GR1-ones(n1,1)*ma)./(ones(n1,1)*sall);
GR2=(GR2-ones(n2,1)*ma)./(ones(n2,1)*sall);
plot(GR1(:,1),GR1(:,2),'red+')
hold on
plot(GR2(:,1),GR2(:,2),'blue+')
hold off
%pause
% variates normalized by st.dev. of ALL
ab=-.11;
ae=.12;
nasteps=11;
phib=1.0;
phie=1.1;
nphisteps=21;
dela=(ae-ab)/(nasteps-1);
delphi=(phie-phib)/(nphisteps-1);
[X,Y]=meshgrid(ab:dela:ae,phib:delphi:phie);

opt=[ab,phib,0.];
for j=1:nphisteps,
phi=(phib+(j-1)*delphi);
alpha=phi*pi;
u(1)=cos(alpha);
u(2)=sin(alpha);
dd1=GR1*u';
dd2=GR2*u';

for i=1:nasteps,
a=ab+(i-1)*dela;
d1=dd1-a;
d2=dd2-a;

crit=0.;
for k=1:n1,
crit=crit+exp(d1(k))/n1;
end
for k=1:n2,
crit=crit+exp(-d2(k))/n2;
end

RES(j,i)=crit;
if [i j]==[1 1]
opt(3)=crit;
end

if crit < opt(3);
opt=[a,phi,crit]; 
end
end
end

contour(X,Y,RES)
ylabel(' phi in units pi')
xlabel(' distance ')
opt
hold on
plot(opt(1),opt(2),'red+')
pause

alpha=opt(2)*pi;
a=opt(1);
u(1)=cos(alpha);
u(2)=sin(alpha);

hold off
colorbar off
plot(GR1(:,1),GR1(:,2),'red+')
axis([-3 3 -3 3])
scale=axis
hold on
plot(GR2(:,1),GR2(:,2),'blue+')
m1=mean(GR1)
m2=mean(GR2)
plot([m1(1),m2(1)],[m1(2),m2(2)],'green')
%pause

d1=GR1*u';
d2=GR2*u';
scor=[0 0 0];
aextr=a;
scoremax=0;
for aa=-.01: .002: .01,
anow=a+aa
ncorrect1=0;
for i=1:n1,
test=d1(i)-anow;
   if test <= 0 
   ncorrect1=ncorrect1+1;
   end
end
ncorrect1

ncorrect2=0;
for i=1:n2,
test=d2(i)-anow;
   if test >= 0
   ncorrect2=ncorrect2+1;
   end
end
ncorrect2
anow
score=[ncorrect1/n1 ncorrect2/n2 (ncorrect1+ncorrect2)/(n1+n2)];
   if score(3) > scoremax
   scoremax=score(3);
   aextr=anow
   scor=score
   end
end

y1=min(ALL(:,2));
y2=max(ALL(:,2));
x1=(a-y1*u(2))/u(1);
x2=(a-y2*u(2))/u(1);
plot([x1,x2],[y1,y2],'yellow')


ab=scale(1);
ae=scale(2);
bb=scale(3);
be=scale(4);

xtekst=ab-0.15*(ae-ab);
ytekst=be+.05*(be-bb);

tekst=sprintf('%s',file1,' n= ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.25*(ae-ab);
tekst=sprintf('%0.5g',n1);
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab);
tekst=sprintf('%s',' correct: ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab);
tekst=sprintf('%0.5g',scor(1));
text(xtekst,ytekst,tekst);


xtekst=xtekst+0.2*(ae-ab);
tekst=sprintf('%s',file2,'   n= ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.25*(ae-ab);
tekst=sprintf('%0.5g',n2);
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab);
tekst=sprintf('%s',' correct: ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab);
tekst=sprintf('%0.5g',scor(2));
text(xtekst,ytekst,tekst);

xtekst=xtekst-0.24*(ae-ab);
ytekst=bb-.1*(be-bb);
tekst=sprintf('%s','total correct: ');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.24*(ae-ab);
tekst=sprintf('%0.5g',scor(3));
text(xtekst,ytekst,tekst);

xtekst=ab-0.15*(ae-ab);
tekst=sprintf('%s','u;a');
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab);
tekst=sprintf('%0.5g',u(1));
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.17*(ae-ab);
tekst=sprintf('%0.5g',u(2));
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.2*(ae-ab);
tekst=sprintf('%0.5g', aextr);
text(xtekst,ytekst,tekst);
