% discrim.m
% testing:
% twoclass distrimination in 2-d space
clf
clear
man=loadasci('constitu.man');
GR1=man(:,3:4);
vrouw=loadasci('constitu.vr');
GR2=vrouw(:,3:4);
n1=size(GR1)
n2=size(GR2)
ALL=[GR1;GR2];
na=size(ALL);
ma=mean(ALL)
sall=std(ALL)
ALL=(ALL-ones(na(1),1)*ma)./(ones(na(1),1)*sall);
GR1=(GR1-ones(n1(1),1)*ma)./(ones(n1(1),1)*sall);
GR2=(GR2-ones(n2(1),1)*ma)./(ones(n2(1),1)*sall);
plot(GR1(:,1),GR1(:,2),'red+')
hold on
plot(GR2(:,1),GR2(:,2),'blue+')
hold off
pause
% variates normalized by st.dev. of ALL
ab=0.;
ae=1;
nasteps=11;
phib=.0;
phie=2.;
nphisteps=11;
dela=(ae-ab)/(nasteps-1);
delphi=(phie-phib)/(nphisteps-1);
[X,Y]=meshgrid(ab:dela:ae,phib:delphi:phie);
opt=[ab,phib,0.];
for j=1:nphisteps,
phi=(phib+(j-1)*delphi);
alpha=phi*pi;
for i=1:nasteps,
a=ab+(i-1)*dela;
u(1)=cos(alpha);
u(2)=sin(alpha);
d1=GR1*u'-a;
d2=GR2*u'-a;
crit=0.;
%for k=1:n1(1),
%crit=crit+exp(d1(k))/n1(1);
%end
%for k=1:n2(1),
%crit=crit+exp(-d2(k))/n2(1);
%end
crit=d1'*d1+d2'*d2;
RES(j,i)=crit;
if [i j]==[1 1]
opt(3)=crit;
end
if crit < opt(3);
opt=[a,phi,crit]; 
end
end
end
%steps=0.9 :.01:  2;
contour(X,Y,RES)
ylabel(' phi in units pi')
xlabel(' distance ')
colorbar
opt
hold on
plot(opt(1),opt(2),'red+')
pause
%a=sqrt(2)/2;
%alpha=.25*pi;
alpha=opt(2)*pi;
a=opt(1);
u(1)=cos(alpha);
u(2)=sin(alpha);
y1=min(ALL(:,2));
y2=max(ALL(:,2));
x1=(a-y1*u(2))/u(1);
x2=(a-y2*u(2))/u(1);
d1=GR1*u'-a;
d2=GR2*u'-a;
crit=0.;
%for k=1:n1(1),
%crit=crit+exp(d1(k))/n1(1);
%end
%for k=1:n2(1),
%crit=crit+exp(-d2(k))/n2(1);
%end
crit=d1'*d1+d2'*d2;

hold off
colorbar off
plot(GR1(:,1),GR1(:,2),'red+')
axis([-3 3 -3 3])
hold on
plot(GR2(:,1),GR2(:,2),'blue+')
plot([x1,x2],[y1,y2],'yellow')
m1=mean(GR1)
m2=mean(GR2)
plot([m1(1),m2(1)],[m1(2),m2(2)],'green')