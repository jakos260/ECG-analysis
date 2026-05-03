% blockinvers.m
% program for interactive inverse problem in which
% the constraint is used to determine a linear
% combination of the zero-space to be added to the 
% unconstrained sulution
clear
AMA=loadasci('ama001.asc');
PSI=loadasci('perqrs.asc');
PHI=loadasci('tankqrs.zme');
L=ones(85,85);
for i=1:85,
L(i,i)=-2;
end
[U,S,V]=svd(AMA);
kt=16;
Sinv=S;
for i=1:kt;
Sinv(i,i)=1/S(i,i);
end
SOL=V(:,1:kt)*Sinv(1:kt,1:kt)*U(:,1:kt)'*PHI;
reldif=norm(AMA*SOL-PHI,'fro')/norm(PHI,'fro')
relerr=norm(SOL-PSI,'fro')/norm(PSI,'fro')
norphi=norm(PHI,'fro')/sqrt(192*50);
norpsi=norm(PSI,'fro')/sqrt(85*50);
norsol=norm(SOL,'fro')/sqrt(85*50);

LX=L*SOL;
LVn=L*V(:,kt+1:85);
BETA=Pinv(LVn)*(-LX);
XX=SOL+V(:,kt+1:85)*BETA;
reldifxx=norm(AMA*XX-PHI,'fro')/norm(PHI,'fro')
relerrxx=norm(XX-PSI,'fro')/norm(PSI,'fro')

