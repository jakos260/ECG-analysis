% reldiffs_mat.m
% function RD=reldiffs_mat(PHI,PSI);
% computes the matrix of rms-differences between all rows of
% matrices PHI and PSI, relative to those of PSI
% 2005-10-20

function RD=reldiffs_mat(PHI,PSI);

[nl nt]=size(PSI);
for i=1:nl,
    denom=norm(PHI(i,:));
    for j=1:nl,
        RD(i,j)=norm(PHI(i,:)-PSI(j,:))/denom;
    end
end





    


