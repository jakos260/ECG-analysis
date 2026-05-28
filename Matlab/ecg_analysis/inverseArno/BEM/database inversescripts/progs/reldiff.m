% reldiff.m
% function rd=reldiff(PHI,PSI)
% 'rms1 rms2 rms1-2 rms1-2/rms2 rms1-2/sqrt(rms1*rms2) rho'  
% RMS-based norms and relative differences between all elements two matrices
function rd=reldiff(PHI,PSI)
M1=PHI;
dim1=size(M1);
M2=PSI;
dim2=size(M2);
ntot=dim1(1)*dim1(2);
rd(1)=norm(M1,'fro')/sqrt(dim1(1)*dim1(2));
rd(2)=norm(M2,'fro')/sqrt(dim2(1)*dim2(2));
rd(3)=norm(M1-M2,'fro')/sqrt(dim1(1)*dim1(2));
rd(4)=rd(3)/rd(2);
rd(5)=rd(3)/sqrt(rd(1)*rd(2));
rho=corrcoef(reshape(PHI,1,ntot),reshape(PSI,1,ntot));
rd(6)=rho(1,2);
