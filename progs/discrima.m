% program discrima

% linear discriminant analysis
% for two groups

echo off
clf

filegr1='constitu.man'
filegr2='constitu.vr'

GR1=loadasci(filegr1);
GR2=loadasci(filegr2);
n1=size(GR1);
n2=size(GR2);
ALL=[GR1;GR2];
na=size(ALL);
sall=std(ALL);
ALL=ALL./(ones(na(1),1)*sall);
GR1=GR1./(ones(n1(1),1)*sall);
GR2=GR2./(ones(n2(1),1)*sall);
% variates normalized by SD

plot(GR1(:,3),GR1(:,4),'blue+')
hold on
plot(GR2(:,3),GR2(:,4),'red+')
m1=mean(GR1)
s1=std(GR1)
m2=mean(GR2)
s2=std(GR2)
mall=mean(ALL);
SHALL=ones(na(1),1)*mall(3:4);
ALLsh=ALL(:,3:4)-SHALL;
axis('equal')
hold on
%COVALL=cov(ALLsh(:,1:2));
[U,S,V]=svd(ALLsh(:,1:2)');
PROJ1=ALLsh*U(:,1)*U(:,1)'+SHALL;
plot(PROJ1(:,1),PROJ1(:,2),'white.')
PROJ2=ALLsh*U(:,2)*U(:,2)'+SHALL;
plot(PROJ2(:,1),PROJ2(:,2),'white.')
plot([m1(3) m2(3)],[m1(4) m2(4)],'yellow')
newdir=(m1(3:4)-mall(3:4))/norm(m1(3:4)-mall(3:4))
sep1=(GR1(:,3:4)-ones(n1,1)*mall(3:4))*newdir';
sep2=(GR2(:,3:4)-ones(n2,1)*mall(3:4))*newdir';
[tdirmeans,df]=ttest2(sep1,sep2)
y=-(m2(3)-m1(3))/(m2(4)-m1(4))*(21.-mall(3))+mall(4);
plot([mall(3) 21],[mall(4) y],'yellow')
a=(GR1(:,3:4)-ones(n1,1)*mall(3:4))*U(:,1);
b=(GR2(:,3:4)-ones(n2,1)*mall(3:4))*U(:,1);
thassenall=ttest2(a,b)
tlengte=ttest2(GR1(:,3),GR2(:,3))
[tgew,df]=ttest2(GR1(:,4),GR2(:,4))
%saveasci('survey.lst',survey);
