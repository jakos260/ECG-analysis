echo off
fin=fopen('perqrs.asc','r');
nrow=fscanf(fin,'%d',[1,1]);
ncol=fscanf(fin,'%d',[1,1]);
PHI=fscanf(fin,'%f',[ncol,nrow]);
fclose(fin);
PHI=PHI';
d=sum(PHI)./nrow;
one=ones(nrow,1);
shift=one*d;
PHIS=PHI-shift;
t=(1:ncol);
hold off
plot(t,PHIS(27:27,:))
hold on
[U,S,V]=svd(PHIS);
sh=(sqrt(nrow)+sqrt(ncol))*0.0042
sl=abs((sqrt(nrow)-sqrt(ncol)))*0.0042
for i=1:ncol,
sn=sl+(ncol+1-i)*(sh-sl)/ncol;
newval=S(i,i)^2-sn^2;
if newval<0,
newval=0;
end
S(i,i)=sqrt(newval);
end
kt=ncol;
PHIT=U(:,1:kt)*S(1:kt,1:kt)*V(:,1:kt)';
plot(t,PHIT(27:27,:))
nd=norm(PHIS-PHIT,'fro');
nr=norm(PHIS,'fro');
reldif=nd/nr
