% bastemp.m
clear
clf
%TAB=loadmat('basicsta.lst');
DATA=loadmat('basicsta.lst');
dim=size(DATA);
npp=dim(1);
nvar=dim(2);
select=[ 2 5 4 18 19 20 21 22 25 26 17 10 7 8 16 15 9 6 12 11 3];
nvar=21;

TAB=DATA(:,select);
COR=corrcoef(TAB);
for i=1:nvar,
COR(i,i)=1;
end
lev=[-.6 -.3  0.3 0.6] 
%mesh(COR)
contour(COR,lev)
axis('square')
%saveasci('corrcoef.lst',COR)
hold on
plot([1 nvar],[1.5 1.5],'w:')
plot([1.5 1.5],[ 1 nvar],'w:')
plot([1 nvar],[nvar-.5 nvar-.5],'w:')
plot([nvar-.5 nvar-.5],[ 1 nvar],'w:')
plot([1 nvar],[17.5 17.5],'w:')
plot([17.5 17.5],[ 1 nvar],'w:')
plot([1 nvar],[10.5 10.5],'w:')
plot([9.5 9.5],[1 nvar],'w:')
plot([12.5 12.5],[1 nvar],'w:')
plot([1 nvar],[12.5 12.5],'w:')

sav=0;
if sav==1,
COR=corrcoef(DATA);
dim=size(COR);
nvar=dim(1);
for j=1:nvar,
[srt,l]= sort(COR(:,j));
SORT(:,j)=srt(:,1);
L(:,j)=l(:,1);
end
saveasci('rhosorted',SORT')
saveasci('Lsorted',L')
end

clf
hold off
scal=[1 nvar 1 nvar]
axis(scal)
axis('square')
axis off
hold on
for i=1:nvar,
for j=1:nvar,
if i==j,
plot(j,22-i,'rx')
else
if COR(i,j) > 0.507,
plot(j,22-i,'r*')
elseif COR(i,j) >0.337,
plot(j,22-i,'r+')
elseif COR(i,j) < -0.337,
text(j,22-i,'-')
elseif COR(i,j) < -0.507,
plot(j,22-i,'b')
else
end
end
end
end
plot([1 nvar],[1.5 1.5],'w:')
plot([1 nvar],[3.5 3.5],'w:')
plot([1 nvar],[10.5 10.5],'w:')
plot([1 nvar],[11.5 11.5],'w:')
plot([1 nvar],[nvar-.5 nvar-.5],'w:')


plot([1.5 1.5],[ 1 nvar],'w:')
plot([13.5 13.5],[1 nvar],'w:')
plot([nvar-.5 nvar-.5],[ 1 nvar],'w:')
plot([nvar-2.5 nvar-2.5],[ 1 nvar],'w:')

text(-0.2,21,'0')
text(-0.5,12,'10')
text(-0.5,1,'20')

text(1,22,'0')
text(10.5,22,'10')
text(21,22,'20')
