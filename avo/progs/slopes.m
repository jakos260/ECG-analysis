% file slopes.m
% uses function diffmat
echo off
clf
psi='psidlss.asc';
PSI=loadasci(psi);
[nrow,ncol]=size(PSI);
rowb=51;
rowe=51;
scal=15;
t=(1:ncol);
subplot(1,2,1)
plot(t,PSI(rowb:rowe,:))
title('measured')
tekst=sprintf('%s','SLOPES');
text(-15,1.1*scal,tekst)
tekst=sprintf('file: %s',psi);
text(40,1.1*scal,tekst)
axis([0 50 -scal scal])
D=diffmat(ncol);
TRUN=1.5*PSI*D;
subplot(1,2,2)
plot(t,TRUN(rowb:rowe,:))
title('differentiated')
axis([0 50 -scal scal])
tekst=sprintf('rows %i',rowb);
text(-25,-1.2*scal,tekst);
tekst=sprintf('to %i',rowe);
text(-5,-1.2*scal,tekst)
[extr,i]=min(TRUN');
j=(1:nrow);
SOL(:,1)=j';
SOL(:,2)=i';
saveasci('psidlss.tim',SOL(1:64,2))
