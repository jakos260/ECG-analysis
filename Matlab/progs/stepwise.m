% stepwise.m
clear
TAB=loadmat('basicstat.lst');
dim=size(TAB);
npp=dim(1);
nvar=dim(2);
TAB(:,1)=ones(npp,1);
COR=corrcoef(TAB);
%variable to be explained:
expl=14;
search=[2 3 4 5 9 10 11 12 15 16 18 19 20 22 23 24];
nsearch=size(search);
for k=1:nsearch(2),
index=search(k);
rho(k)=COR(expl,index);
end
[y i]=sort(rho');
[y search(i)']
pause
next=20
nor=0
if nor==1,
% normalize for sex
for j=3:nvar,
nmales=0;
nfemales=0;
summales=0;
sumfemales=0;
for i=1:npp,
if TAB(i,2) ==0,
nfemales=nfemales+1;
sumfemales=sumfemales+TAB(i,j);
else
nmales=nmales+1;
summales=summales+TAB(i,j);
end
end
total=summales+sumfemales;
gem=total/npp;
gemmales=summales/nmales;
gemfemales=sumfemales/nfemales;

for i=1:npp,
 if  TAB(i,2) ==0,
    TAB(i,j)=TAB(i,j)*gem/gemfemales;
 else 
    TAB(i,j)=TAB(i,j)*gem/gemmales;
 end
end
end

end
COR=corrcoef(TAB);
%variable to be explained:
search=[2 3 4 5 9 10 11 15 16 18 19 22 23 24];
nsearch=size(search);

rho=zeros(1,nsearch);
for k=1:nsearch(2),
index=search(k);
rho(k)=COR(expl,index);
end
[y l]=sort(rho');
[y search(l)']
refnorm=norm(TAB(:,expl));
MAT(:,1)=TAB(:,1);
MAT(:,2)=TAB(:,next);
MAT(:,3)=[];
RD=norm(TAB(:,expl)-MAT*pinv(MAT)*TAB(:,expl))/refnorm
search=[2 3 4 5 9 10 11 12 15 16 18 19 22 23 24];
nsearch=size(search);
MAT(:,1)=TAB(:,1);
MAT(:,2)=TAB(:,11);
for k=1:nsearch(2);
MAT(:,3)=TAB(:,search(k));
res(k)=norm(TAB(:,expl)-MAT*pinv(MAT)*TAB(:,expl))/refnorm;
end
[z n]=sort(res');
[z search(n)']
