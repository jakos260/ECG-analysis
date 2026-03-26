% discrim2.m
% twoclass discrimination in 2-d space
clf
clear
file1='constitu.man';
man=loadasci(file1);
GR1=man(:,3:4);
file2='constitu.vr';
vrouw=loadasci(file2);
GR2=vrouw(:,3:4);

n1=size(GR1)
n2=size(GR2)
CV1=cov(GR1)
CV2=cov(GR2)
plot(GR1(:,1),GR1(:,2),'red+')
hold on
plot(GR2(:,1),GR2(:,2),'blue+')
hold off
scal=axis
%pause
m1=mean(GR1)
m2=mean(GR2)
sig11=sqrt(CV1(1,1))
sig12=sqrt(CV1(2,2))
sig21=sqrt(CV2(1,1))
sig22=sqrt(CV2(2,2))
rho1=CV1(1,2)/(sig11*sig12)
rho2=CV2(1,2)/(sig21*sig22)

ab=scal(1);
ae=scal(2);
nasteps=21;
bb=scal(3);
be=scal(4);
nbsteps=21;
dela=(ae-ab)/(nasteps-1);
delb=(be-bb)/(nbsteps-1);
[X,Y]=meshgrid(ab:dela:ae,bb:delb:be);
inv1=inv(CV1);
inv2=inv(CV2);
for i=1:nbsteps,
x2=bb+(i-1)*delb;
for j=1:nasteps,
x1=ab+(j-1)*dela;
RES(i,j)=-(m1-[x1 x2])*inv1*(m1-[x1 x2])'+(m2-[x1 x2])*inv2*(m2-[x1 x2])';
end
end
classif=log(n2(1)^2/n1(1)^2*(1-rho1^2)/(1-rho2^2))
delta=0
contour(X,Y,RES,[classif classif+delta])
ylabel(' gewicht')
xlabel(' lengte ')
hold on
plot(GR1(:,1),GR1(:,2),'red+')
plot(GR2(:,1),GR2(:,2),'blue+')
ncorrect1=0
for i=1:n1(1),
test=-(m1-GR1(i,:))*inv1*(m1-GR1(i,:))'+(m2-GR1(i,:))*inv2*(m2-GR1(i,:))';
if test >= classif+delta
ncorrect1=ncorrect1+1;
end
end
ncorrect2=0
for i=1:n2(1),
test=-(m1-GR2(i,:))*inv1*(m1-GR2(i,:))'+(m2-GR2(i,:))*inv2*(m2-GR2(i,:))';
if test <= classif+delta
ncorrect2=ncorrect2+1;
end
end
[ncorrect1/n1(1) ncorrect2/n2(1) (ncorrect1+ncorrect2)/(n1(1)+n2(1))]


xtekst=ab-0.15*(ae-ab)
ytekst=be+.05*(be-bb)

tekst=sprintf('%s',file1,' n= ')
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.25*(ae-ab)
tekst=sprintf('%0.5g',n1(1))
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab)
tekst=sprintf('%s',' correct: ')
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab)
tekst=sprintf('%0.5g',ncorrect1/n1(1))
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.2*(ae-ab)
tekst=sprintf('%s',file2,'   n= ')
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.25*(ae-ab)
tekst=sprintf('%0.5g',n2(1))
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.07*(ae-ab)
tekst=sprintf('%s',' correct: ')
text(xtekst,ytekst,tekst);

xtekst=xtekst+0.14*(ae-ab)
tekst=sprintf('%0.5g',ncorrect2/n2(1))
text(xtekst,ytekst,tekst);

