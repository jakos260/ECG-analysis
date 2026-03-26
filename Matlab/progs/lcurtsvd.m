echo off
clf
phi='tankqrs.zmean'
ama='ama001.asc'
psiref='perqrs.zmean'
PHI=loadasci(phi);
[nrowphi,ncolphi]=size(PHI);
t=(1:ncolphi);
A=loadasci(ama);
[nrowa,ncola]=size(A);
PSI=loadasci(psiref);
subplot(1,2,1)
plot(t,PSI(19:23,:))
title('measured')
tekst=sprintf('%s','LCURVETSVD');
text(-15,16,tekst)
space='     ';
tekst=sprintf('%s',[phi,space,ama,space,psiref]);
text(-15,-17,tekst)
axis([0 50 -15 15])
[nrowpsi,ncolpsi]=size(A);
[U,S,V]=svd(A);
kt=17;
for i=1:kt,
Sinv(i,i)=1./S(i,i);
end
for i=1:kt,
SOL=V(:,1:i)*Sinv(1:i,1:i)*U(:,1:i)'*PHI;
re(i)=norm(SOL-PSI,'fro')/norm(PSI,'fro');
end
re(1:kt)
subplot(1,2,2)
plot(t,SOL(19:23,:))
title('inverse;   kt=')
tekst1=sprintf('%i',kt);
text(33,16,tekst1)
tekst1=sprintf('rd=%0.5g',re(kt));
text(-15,-17,tekst1)
axis([0 50 -15 15])
%saveasci('solr13.asc',SOL)