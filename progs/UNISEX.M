% unisex.m
clear
clf
TAB=loadmat('basicsta.lst');
dim=size(TAB);
npp=dim(1);
nvar=dim(2);

for k=1:nvar,
nm=0;
nf=0;
 for i=1:npp
   if TAB(i,2)==0,
     nf=nf+1;
     yf(nf,1)=TAB(i,k);
   else
     nm=nm+1;
     ym(nm,1)=TAB(i,k);
   end
 end
 
%normalize for equal means for both sexes
mmean=mean(ym);
fmean=mean(yf);
allmean=mean(TAB(:,k));

 for i=1:npp
if k> 2,
   if TAB(i,2)==0,
     TAB(i,k)=TAB(i,k)*allmean/fmean;
   else
     TAB(i,k)=TAB(i,k)*allmean/mmean;
   end
 end
end 
end
saveasci('unisex.lst',TAB)
