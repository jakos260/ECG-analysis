% variabil.m
% find measures of (relative) variability for 64 lead positions
clear
allfiles=loadchars('normals.labs');
dim=size(allfiles);
ni=dim(1);
TABEL=loadasci('basicstat.lst');
areas=TABEL(:,19);

pp=allfiles(1,3:4);
file='pp08/pp08.dip';
file(3:4)=pp(1:2);
file(8:9)=pp(1:2);
PHI=loadmat(file);
PHI=PHI';
dim=size(PHI);
nt=dim(2);
nl=dim(1);
%nt=nt/2;
ALLPHI=zeros(nl*ni,nt);

for i=1:ni,
pp=allfiles(i,3:4);
file='pp08/pp08.dip'
file(3:4)=pp(1:2);
file(8:9)=pp(1:2);
nbeg=1+(i-1)*nl;
nend=i*nl;
B=loadmat(file);
B=B';
%MATT=pinv(B')*REF
A=B/areas(i);
%A=B;
%A=A/norm(A,'fro');
ALLPHI(nbeg:nend,1:nt)=A(1:nl,1:nt);
end

for l=1:nl,
for t=1:nt,

for i=1:ni,
n=(i-1)*nl+l;
mn(i)=ALLPHI(n,t);
end

GEM(l,t)=mean(mn);
SD(l,t)=std(mn);

SSQDIF(l,t)=(SD(l,t)*sqrt(ni-1)).^2;
SSQ(l,t)=norm(mn).^2;
RVARLT(l,t)=sqrt(SSQDIF(l,t)./SSQ(l,t));
end
end

%saveasci('rvarltTWAVE.asc',RVARLT)

for l=1:nl,
svarl(l)=sum(SSQDIF(l,:));
nvarl(l)=sum(SSQ(l,:));
rvarl(l)=sqrt(svarl(l)./nvarl(l));

rmsrvl(l)=sqrt(svarl(l))/(nt);
rmsnvl(l)=sqrt(nvarl(l))/(nt);
end

%saveasci('rvarlTWAVE.asc',rvarl)
%saveasci('rmsrvlTWAVE.asc',rmsrvl)
%saveasci('rmsnvlTWAVE.asc',rmsnvl)

for t=1:nt,
svart(t)=sum(SSQDIF(:,t));
nvart(t)=sum(SSQ(:,t));
rvart(t)=sqrt(svart(t)./nvart(t));
end

%saveasci('rvartTWAVE.asc',rvart)

rvar1=sqrt(sum(svarl)/sum(nvarl))
rvar2=sqrt(sum(svart)/sum(nvart))
rms=sqrt(sum(nvart)/(nl*ni*nt))
%saveasci('meansTWAVE.asc',GEM)
%saveasci('sdvalsTWAVE.asc',SD)
