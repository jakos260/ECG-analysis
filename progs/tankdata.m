echo off
PHI=loadasci('tankqrs.zmean');
[nrowphi,ncolphi]=size(PHI)
t=(1:ncolphi);
[U,S,V]=svd(PHI);
kt=40;
loadasci('ama001tuned.asc');
[nrowa,ncola]=size(A);
PSI=loadasci('perqrs.zmean');
hold off
plot(t,PSI(25:27,:))
hold on
[nrowpsi,ncolpsi]=size(A);
UTA=U(:,1:kt)'*A;
[nr,nc]=size(UTA)
RHS=S(1:kt,1:kt)*V(:,1:kt)';
[nr,nc]=size(RHS)
[U,S,V]=svd(UTA);
%kt=40;
for i=1:kt,
Sinv(i,i)=1./S(i,i);
end
for i=1:kt,
SOL=V(:,1:i)*Sinv(1:i,1:i)*U(:,1:i)'*RHS;
re(i)=norm(SOL-PSI,'fro')/norm(PSI,'fro');
end
plot(t,SOL(25:27,:))
re(1:kt)'
%saveasci('solr13.asc',SOL)