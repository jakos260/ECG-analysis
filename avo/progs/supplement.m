%supplement.m
%inversqrs.m
% 01/02/08
% inverse calculation of
% timing at Sh
% version automatic lambda evaluation

% epicardial timing is frozen

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings aimed at treating: DEPOLARIZATION
ama='udl2lds.asc';
potentials='../udlinv/qrstank.bczme';
refmat='../udlinv/udl2epi.new';
refpots='../udlinv/qrsepi.bczme';
subsamp=1;
surflapl='../udlinv/surflaplhart.asc';
Q=1;
%initims='../rank20.init';
initims='../udlinv/qrsepi04.tim';
k=0;% identifies rank of initims
win=2;
regpar=0.01;
logfil='logqrs.004';
outfil='qrs004.tim';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load input data
% transfer
A=loadmat(ama);
A=A';
AREF=loadmat(refmat);
AREF=AREF';

% potentials
PSI=loadmat(potentials);
PSI=PSI';
psidim=size(PSI);

REFPOTS=loadmat(refpots);
REFPOTS=REFPOTS';

PHI=PSI(:,1:subsamp:psidim(2));
REFPOTS=REFPOTS(:,1:subsamp:psidim(2));

phidim=size(PHI);
nl=phidim(1);
nt=phidim(2);
adim=size(A);
n=adim(2);
normphi=norm(PHI,'fro');
normref=norm(REFPOTS,'fro');
% surface laplacian
SL=loadmat(surflapl);
SL=SL';

%prepare data
ATA=A'*A;
LTL=SL'*SL;
w1=sum(PHI');
nw1=norm(w1);
w1=-w1';
t=1:nt;
t=t';
w2=-2*PHI*t;
nw2=norm(w2);
% shape of S-curve
span=1000; %size of look-up table for shape action potential 

if Q==1,
% use two connected parabolae to describe source strength
% during activation
map=getmapQ(span,span);
end
mapgrad=zeros(1,span);
ttemp=1:span;
mapgrad(1,1:span/2)=4*ttemp(1:span/2)/span;
mapgrad(span/2+1:span)=2-4*ttemp(1:span/2)/span;

if Q==0,
% this part estimates shape of 
% action potential from mean of signals on the outside
% of the tank/torso ; use for treating T-wave only
   map=1-getmap0(PHI);
   x=1:nt;
   clf
plot(x,map)
   % get tsteep=point of steepest gradient of S-curve
   % operator for time differentiation
   nd=nt;
   mapgrad=zeros(1,nd);
   mapgrad(1)=map(2)-map(1);
   for t=2:nd-1;
   mapgrad(t)=(map(t+1)-map(t-1))/2;
   end
   mapgrad(nd)=map(nd)-map(nd-1);

   [steep tsteep]=max(abs(mapgrad));
   middle=round(nt/2);
   nextra=abs(middle-tsteep)
   nmap=nt+nextra
	if middle > tsteep,
          tmpmap(1,1:nextra)=zeros(1,nextra);
	  tmpmap(1,nextra+1:nmap)=map(1,1:nt);
        end
	if middle < tsteep,
	  tmpmap(1,1:nt)=map(1,1:nt);
          tmpmap(1,nt+1:nmap)=ones(1,nextra);
        end
   % operator for time differentiation
        nd=nmap;
        mapgrad=zeros(1,nd);
        mapgrad(1)=tmpmap(2)/2;
        for t=2:nd-1;
           mapgrad(t)=(tmpmap(t+1)-tmpmap(t-1))/2;
        end
        mapgrad(nd)=(1-tmpmap(nd-1))/2;
   x=1:nmap;
   xx=1:(nmap-1)/span:nmap;
   dim=size(xx)
   span=dim(2)
   map=spline(x,tmpmap,xx);
   mapgrad=(nmap-1)*spline(x,mapgrad,xx);
   clf
end

ik=0;
iter=0;

% load initial estimate
tims=loadasci(initims);
tims=tims/subsamp;
if mean(tims) < 0.1*nt,
tims=tims+ones(n,1)*nt/2;
end
tims=max(tims,win/2*ones(n,1));
tims=min(tims,(nt-win/2)*ones(n,1));

%tims=ones(n,1)*mean(tims);
% use this option only when treating T-wave 
saveasci('listen.tim',tims)

winj=win*ones(n,1);
l=0;
saveasci('listen.tim',tims)
timsnew=tims;
S=getSmap(nt,timsnew,map',winj);
PHIA=A*S;
RES=PHI-PHIA;
res=norm(RES,'fro');
rd=res/normphi;

rdref=norm(REFPOTS-AREF*S,'fro')/normref;

reg=norm(SL*timsnew);
w1inv=A*timsnew;
relres1=norm(w1inv-w1)/nw1;
w2inv=A*timsnew.^2;
relres2=norm(w2inv-w2)/nw2;
ik=ik+1;
if reg <=100, reg=7000/subsamp;end
mu=res/reg*regpar; 
tresiter=sqrt(rd^2+(mu*reg/normphi)^2);
temp=[iter k mu min(timsnew) max(timsnew)];
temp=[temp  relres1 relres2 rd reg/1000 tresiter rdref]
RESNOW(ik,:)=temp;

labopt=0.1;
lambdamax=labopt*1.e+4;
lambdamin=labopt*1.e-4;

stop1=0;
%%%%%%%%%%%%%%%%%%%%
% iterative approach to solving non linear parameter estimate
while stop1 ==0,
'start next iteration'
iter=iter+1
lamb=labopt;
tims=timsnew;

%compute Sprime and Sprime*Sprime'
Sprime=-getSprime(nt,tims,mapgrad',winj);
for t=1:nt,
Sprime(:,t)=Sprime(:,t)./winj(:,1);
end
SST=Sprime*Sprime';
%'compute GTG'
GTG=ATA.*SST;

%'compute GTres'
gtres=sum(A.*(RES*Sprime'));
gtres=gtres';

lab=ones(1,3);
lab(1)=lamb/2;
lab(2)=lamb;
lab(3)=2*lamb;

for i=1:3,
        lamb=lab(i);
        GTGL=GTG+lamb^2*eye(n)+mu^2*LTL;
        deltau=inv(GTGL)*(gtres-mu^2*LTL*tims);
deltau(358:620)=0;
        testtims=tims+deltau;
        reg=norm(SL*testtims);
        %'compute PHIA'
        S=getSmap(nt,testtims,map',winj);
        PHIA=A*S;
        RES=PHI-PHIA;
        rd=norm(RES,'fro')/normphi;
        tres(i)=sqrt(rd^2+(mu*reg/normphi)^2);
end

% search for lambda value that produces smallest total residual (tres)
%%%%%%%%%
stop2=0;
while stop2 ==0,
  labopt=minpar(tres,lab);
  [lab tres labopt];
  labopt=max(lambdamin,labopt);
  labopt=min(lambdamax,labopt);
  %'        lab(1)--lab(3)               labopt'
  %[lab labopt]
    if labopt <= lab(1),
    %   'shift testinterval to the left'
       lab(3)=lab(2);
       tres(3)=tres(2);
       lab(2)=lab(1);
       tres(2)=tres(1);
       lab(1)=lab(1)/2;
       lamb=lab(1);
       GTGL=GTG+lamb^2*eye(n)+mu^2*LTL;
       deltau=inv(GTGL)*(gtres-mu^2*LTL*tims);
deltau(358:620)=0;
       timstemp=tims+deltau;
       reg=norm(SL*timstemp);
       %'compute PHIA'
       S=getSmap(nt,timstemp,map',winj);
       PHIA=A*S;
       RES=PHI-PHIA;
       rd=norm(RES,'fro')/normphi;
       tres(1)=sqrt(rd^2+(mu*reg/normphi)^2);
    end

     if labopt > lab(1) & labopt < lab(3),
      % 'test testinterval'
       lamb=labopt;
       GTGL=GTG+lamb^2*eye(n)+mu^2*LTL;
       deltau=inv(GTGL)*(gtres-mu^2*LTL*tims);
deltau(358:620)=0;
       testtims=tims+deltau;
       reg=norm(SL*testtims);
       %'compute PHIA'
       S=getSmap(nt,testtims,map',winj);
       PHIA=A*S;
       RES=PHI-PHIA;
       rd=norm(RES,'fro')/normphi;
       trestest=sqrt(rd^2+(mu*reg/normphi)^2);
%'      iter     lab(1)-lab(3)               lambtest    trestest   tresiter'
%[iter lab lamb trestest tresiter]

       if trestest >= tresiter, labopt=4*lab(1); end
       if trestest < tresiter,
          labopt=lamb;
          tresiter=trestest; 'break1'; break;
       end
     end

     if labopt >= lab(3),
        rhs=1;
           while rhs==1,
              %'shift testinterval to the right'
              lab(1)=lab(2);
              tres(1)=tres(2);
              lab(2)=lab(3);
              tres(2)=tres(3);
              lab(3)=2*lab(3);
              lamb=lab(3);
              if lamb>=lambdamax, break; end
              GTGL=GTG+lamb^2*eye(n)+mu^2*LTL;
              deltau=inv(GTGL)*(gtres-mu^2*LTL*tims);
deltau(358:620)=0;
              timstemp=tims+deltau;
              reg=norm(SL*timstemp);
              %'compute PHIA'
              S=getSmap(nt,timstemp,map',winj);
              PHIA=A*S;
              RES=PHI-PHIA;
              rd=norm(RES,'fro')/normphi;
              tres(3)=sqrt(rd^2+(mu*reg/normphi)^2);
          %'lamb tres(3) tres(3)-tresiter'
	  %    [lamb tres(3) tres(3)-tresiter]
          if (tresiter-tres(3))/tresiter >=1.e-4, rhs=0;
	  lab(1)=lab(2);
	  lab(2)=lab(3);
	  lab(3)=2*lab(3);
	  end
        end

        if lamb >=lambdamax, break; end
     end
end

%%%%%%%%%
if lamb>lambdamax, 'break2', break; end
timsnew=testtims;
saveasci('listen.tim',timsnew)
w1inv=A*timsnew;
relres1=norm(w1inv-w1)/nw1;
w2inv=A*timsnew.^2;
relres2=norm(w2inv-w2)/nw2;
ik=ik+1;
rdref=norm(REFPOTS-AREF*S,'fro')/normref;
temp=[iter lamb mu min(timsnew) max(timsnew)];
temp=[temp  relres1 relres2 rd reg/1000 tresiter rdref]
RESNOW(ik,:)=temp;
end
%%%%%%%%%%%%%%%%%%%%

saveasci(outfil,timsnew)
saveasci(logfil,RESNOW)

