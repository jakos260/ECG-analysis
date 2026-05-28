% getdist.m
clear
allfiles=loadchars('ppees.lst')
EXTR=zeros(25,6);
for p=1:25,
file='pp08/pp08tor.dhk3';
file(1:9)= allfiles(p,:)
[TVER,TTRI]=loadtri(file);
file='pp08/pp08.mels';
file(1:9)= allfiles(p,:)
LAMU=loadmat(file);
LAMU=LAMU';
XYZ=lamu2xyz(LAMU,TTRI,TVER);
file='pp08/pp08har.dhk';
file(1:9)= allfiles(p,:)
[HVER,HTRI]=loadtri(file);
dim=size(TVER);
npt=dim(1);
dim=size(HVER);
nph=dim(1);
dim=size(XYZ);
nelecs=dim(1);
DIST=zeros(npt,nph);
for i=1:npt,
for j=1:nph,
DIST(i,j)=norm(TVER(i,:)-HVER(j,:));
end
end
extr=extremes(DIST)
'tussen hart en torso'
EXTR(p,1:3)=extr(1,1:3);
DIST=zeros(nelecs,nph);
for i=1:nelecs,
for j=1:nph,
DIST(i,j)=norm(XYZ(i,:)-HVER(j,:));
end
end
extr=extremes(DIST)
EXTR(p,4:6)=extr(1,1:3);
'tussen hart en electroden'
end
