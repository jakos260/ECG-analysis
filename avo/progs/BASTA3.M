% basta3m
clear
DATA=loadasci('angles.lst');
dim=size(DATA);
npp=dim(1);
nvar=dim(2);
TAB=loadasci('basicsta.lst');
STAT=zeros(nvar,12);

for j=1:nvar
meanall=mean(DATA(:,j));
stdall=std(DATA(:,j));
small=min(DATA(:,j));
big=max(DATA(:,j));
STAT(j,1:4)=[meanall stdall small big];

nm=0;
nf=0;
for i=1:npp
if TAB(i,2)==0,
nf=nf+1;
xf(nf)=DATA(i,j);
else
nm=nm+1;
xm(nm)=DATA(i,j);
end
end

fmean=mean(xf);
fstd=std(xf);
small=min(xf);
big=max(xf);
STAT(j,5:8)=[fmean fstd small big];

mmean=mean(xm);
mstd=std(xm);
small=min(xm);
big=max(xm);
STAT(j,9:12)=[mmean mstd small big];
end
STAT
saveasci('anglesst.lst',STAT)
