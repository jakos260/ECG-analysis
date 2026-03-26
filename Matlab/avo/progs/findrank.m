echo off
PHI=loadasci('perqrs.zmean');
[nrow,ncol]=size(PHI);
d=sum(PHI)./nrow;
one=ones(nrow,1);
shift=one*d;
PHIS=PHI-shift;
phinorm=norm(PHIS,'fro');
phisrms=phinorm/sqrt(nrow*ncol)
t=(1:ncol);
hold off
%plot(t,PHIS(25:27,:))
hold off
[U,S,V]=svd(PHIS);
ss=sum(S);
%plot(t,ss)
hold off
kt=15;
PHIT=U(:,1:kt)*S(1:kt,1:kt)*V(:,1:kt)';
%plot(t,PHIT(25:27,:),'+')
nd=norm(PHIS-PHIT,'fro');
reldif=nd/phinorm
phitRMS=norm(PHIT,'fro')/sqrt(nrow*ncol)
%saveasci('perqrs.r10',PHIT)
st=svd(PHIT);
plot(t,st,'red')
hold on
noise=randn(nrow,ncol);
noisefac=.1
noise=noise.*(phitRMS*noisefac);
sn=svd(noise);
plot(t,sn,'magenta')
PHITN=PHIT+noise;
rd0=norm(PHITN-PHIT,'fro')/norm(PHIT,'fro')
[U,S,V]=svd(PHITN);
ktr=15;
RECTSVD=U(:,1:ktr)*S(1:ktr,1:ktr)*V(:,1:ktr)';
rd1=norm(RECTSVD-PHIT,'fro')/norm(PHIT,'fro')
stn=sum(S); 
plot(t,stn,'black')
[U,S,V]=svd(PHITN);
phitnRMS=norm(PHITN,'fro')/sqrt(nrow*ncol)
%nf=.05
sigma=ncol*(S(ncol-30,ncol-30)-S(ncol,ncol))/(60*sqrt(ncol))
%phitnRMS*nf;
%sl=sigma*abs(sqrt(nrow)-sqrt(ncol))
sl=S(ncol,ncol);
sh=(sqrt(nrow)+sqrt(ncol))*sigma;
for i=1:ncol,
snew=sl+(ncol+1-i)*(sh-sl)/ncol;
newval=S(i,i)^2-snew^2;
if newval<0,
newval=0;
end
S(i,i)=sqrt(newval);
end
ssnew=sum(S);
plot(t,ssnew,'green')
RECNEW=U*S*V';
rd2=norm(RECNEW-PHIT,'fro')/norm(PHIT,'fro')
rd3=norm(RECNEW-PHITN,'fro')/norm(PHITN,'fro')
nfest=sigma/phitnRMS
rankest=(nfest-rd3)/nfest*ncol