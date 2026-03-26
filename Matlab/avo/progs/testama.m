echo off
clf
phi='tankqrs.zme';
ama='ama001.asc';
psiref='perqrs.zme';
PHI=loadasci(phi);
[nrowphi,ncolphi]=size(PHI);
t=(1:ncolphi);
A=loadasci(ama);
[nrowa,ncola]=size(A);
PSI=loadasci(psiref);
[nrowpsi,ncolpsi]=size(PSI);
%subplot(1,2,1)
%plot(t,PSI(19:23,:))
%title('measured')
%tekst=sprintf('%s','TESTAMA');
%text(-15,16,tekst)
%space='     ';
%tekst=sprintf('%s',[phi,space,ama,space,psiref]);
%text(-15,-17,tekst)
%axis([0 50 -15 15])

[U,S,V]=svd(PSI);

% plot singular values of PSI
subplot(1,3,1)
s=sum(S);
plot(t,s(1:ncolpsi))
axis([0 50 0 150])

SOL=A*U;
% SOL expresses which parts of the basic PSI patterns (U)
% seeps through the A filter

SUB=SOL(:,1:50);
% plot SUB=SOL(:,1:50);

for i=1:ncolphi,
no(i)=norm(SUB(:,i));
[i,no(i)];
end

% plot the norms of these patterns
subplot(1,3,2)
plot(t,no(1:ncolphi))
axis([0 50 0 .2])

%COV=SUB'*SUB;
%for i=1:50,
%d(i)=norm(SUB(:,i));
%end
%D=d'*d;
%COR=COV./D;
%[X,Y]=meshgrid(1:1:20,1:1:20);
%subplot(1,3,3)
%steps=linspace(-1. , 1.,11);
%contour(X,Y,COR(1:20,1:20),steps)
%colorbar

% 
kt=17;
SUB=SOL(:,1:kt);
W=SUB\PHI;
REPR=SOL(:,1:kt)*W(1:kt,:);
for i=1:50,
err(i)=norm(REPR(:,i)-PHI(:,i))./norm(PHI(:,i));
[i,err(i)];
end
rd=norm(REPR-PHI,'fro')/norm(PHI,'fro')

% plot representation error based on kt patterns as a function of t
subplot(1,3,3)
plot(t,err)
axis([0 50 0. 1.])
tekst=sprintf('kt=%i',kt);
text(55., -.05,tekst)
tekst=sprintf('rd=%0.5g',rd);
text(55.,0.9,tekst)
USUB=U(:,1:kt);

% compute inverse errorbased on truncating  at kt 
PSIINV=USUB*W;
reler=norm(PSIINV-PSI,'fro')/norm(PSI,'fro')
tekst=sprintf('re=%0.5g',reler);
text(55.,0.5,tekst)
