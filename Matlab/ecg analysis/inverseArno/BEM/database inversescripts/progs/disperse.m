% dispersie.m
% determine various measures for the interindividual variation in a collection
% of BSMovies.
clear

allfiles=loadchars('normals.labs');
nfiles=size(allfiles);
for i=1:nfiles(1),
file1=allfiles(i,:)
PHI=loadmat(file1);
for j=1:nfiles(1),
file2=allfiles(j,:);
PSI=loadmat(file2);
rd=reldiff(PSI,PHI);
RD(i,j)=rd(4);
SCD(i,j)=rd(5);
RHO(i,j)=rd(6);
end
end
saveasci('reldiffs.asc',RD)
saveasci('scaldiffs.asc',SCD)
saveasci('rhos.asc',RHO)
