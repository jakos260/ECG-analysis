echo off
clf
PSI=loadasci('perqrs.zmean');
[nrow,ncol]=size(PSI)
t=(1:ncol);
subplot(1,2,1)
plot(t,PSI(19:23,:))
title('measured')
tekst=sprintf('%s','TSVD ')
text(-15,17,tekst)
tekst=sprintf('%s','file: perqrs.zmean')
text(50,17,tekst)
axis([0 50 -15 15])
[U,S,V]=svd(PSI);
kt=4;
TRUN=U(:,1:kt)*S(1:kt,1:kt)*V(:,1:kt)';
re=norm(TRUN-PSI,'fro')/norm(PSI,'fro');
re
subplot(1,2,2)
plot(t,TRUN(19:23,:))
title('truncated at:  kt=')
tekst=sprintf('%i',kt')
text(35,16,tekst)
axis([0 50 -15 15])
%saveasci('solr13.asc',SOL)