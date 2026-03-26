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
[nrowpsi,ncolpsi]=size(PSI);
for j=1:ncolphi,
rms(j)=norm(PHI(:,j))/sqrt(nrowphi);
end
ib=1
ie=7;
subplot(1,2,1)
plot(t,PSI(ib:ie,:))
title('measured')
tekst=sprintf('%s','LCURVEDLT');
text(-15,16,tekst)
space='     ';
tekst=sprintf('%s',[phi,space,ama,space,psiref]);
text(-15,-17,tekst)
axis([0 50 -15 15])
[U,S,V]=svd(A);
krank=min(nrowa,ncola)
%lambda=.8 =opt
lambda=.8; 
for j=1:ncolpsi,
ala(j)=lambda*0.003/(0.003+rms(j));
for i=1:krank,
Sinv(i,i)=S(i,i)/(S(i,i)^2+ala(j)^2);
end
SOL(:,j)=V*Sinv*U(:,1:krank)'*PHI(:,j);
end
%SOL=V(:,1:i)*Sinv(1:i,1:i)*U(:,1:i)'*PHI;
%end
re=norm(SOL-PSI,'fro')/norm(PSI,'fro');

subplot(1,2,2)
plot(t,SOL(ib:ie,:))
title('DLS-inverse')
tekst1=sprintf('  lambda %0.5g',lambda);
text(33,16,tekst1)
tekst1=sprintf('rd=%0.5g',re);
text(-15,-17,tekst1)
tekst2=sprintf('rows %d to %d',ib,ie);
text(30,-17,tekst2)
axis([0 50 -15 15])
%saveasci('solr13.asc',SOL)
peak=max(max(SOL))
through=min(min(SOL))
