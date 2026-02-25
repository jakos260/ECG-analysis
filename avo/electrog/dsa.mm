% dsa.m
% function  [SA,index]=solida(VER,ITRI,obs,small)
% compute the distributed solid angles at observation point obs
% subtended by all elements of a  set of triangles
% 2003-02-23

% if index==0: obs not in plane of triangle
   % if any solid angle (sa) is smaller than small, one third of the
   % of the sa is assigned to every vertex of the triangle
   % basic version (IEEE,BME_30, 1983, pp. 125-126),
   % else, a further computation yields
   % the values according to de Munck: IEEE Trans BME_39,1993, pp986-990; eqn 19.

% else SA(*,1:3)=0; obs in plane of the triangle; moreover:
   % index=1: obs external point of triangle
   % index=2: obs on edge of triangle;
   % index=3: obs internal point of triangle
   % index=4: obs is node of triangle

function  [SA,index]=dsa(VER,ITRI,obs,small)
[nver,idum]=size(VER);
VER=VER-ones(nver,1)*obs;
R1=[VER(ITRI(:,1),1),VER(ITRI(:,1),2),VER(ITRI(:,1),3)];
R2=[VER(ITRI(:,2),1),VER(ITRI(:,2),2),VER(ITRI(:,2),3)];
R3=[VER(ITRI(:,3),1),VER(ITRI(:,3),2),VER(ITRI(:,3),3)];
LR=[ norm3d(R1) norm3d(R2) norm3d(R3)];
R2crR3=cross(R2,R3);
[ntri,dum]=size(ITRI);
index=zeros(1,ntri);
% blockproduct
block=dots(R1,R2crR3);
DOTS=[dots(R1,R2) dots(R2,R3) dots(R3,R1)];
denom=LR(:,1).*LR(:,2).*LR(:,3) + ...
      LR(:,1).*DOTS(:,2)+LR(:,2).*DOTS(:,3) + LR(:,3).*DOTS(:,1);
index(abs(block)<=eps & denom >   eps)=1;
index(abs(block)<=eps & denom >= -eps & denom <=eps)=2;
index(abs(block)<=eps & denom <  -eps)=3;
index(LR(:,1)<=eps|LR(:,2)<=eps|LR(:,3)<=eps)=4;
sa=-2*atan2(block,denom)';
%[block denom sa' index']

sa(abs(block)<=eps)=0;
SA=ones(3,1)*sa/3;
TEST=[1:ntri;sa];
k=TEST(1,abs(TEST(2,:))>small);

if isempty(k)==0,
% further computation required for triangles: k
S1=R2(k,:)-R1(k,:);
S2=R3(k,:)-R2(k,:);
S3=R1(k,:)-R3(k,:);
LS=[norm3d(S1) norm3d(S2) norm3d(S3)];
DOTS=DOTS(k,:);

R2crR3=R2crR3(k,:);
R3crR1=cross(R3(k,:),R1(k,:));
R1crR2=cross(R1(k,:),R2(k,:));
NT=R1crR2+R2crR3+R3crR1;
asq=dots(NT,NT);
NT=(sa(k)'*ones(1,3)).*NT;

% compute the gamma vector
m=[2 3 1];
NUM=LR(k,:).*LS+DOTS-LR(k,:).*LR(k,:);
DNOM=LR(k,m).*LS-DOTS+LR(k,m).*LR(k,m);
GAM=((block(k)*ones(1,3)).*log(NUM./DNOM))./LS;

GV(:,1)=S1(:,1).*GAM(:,1)+S2(:,1).*GAM(:,2)+S3(:,1).*GAM(:,3);
GV(:,2)=S1(:,2).*GAM(:,1)+S2(:,2).*GAM(:,2)+S3(:,2).*GAM(:,3);
GV(:,3)=S1(:,3).*GAM(:,1)+S2(:,3).*GAM(:,2)+S3(:,3).*GAM(:,3);

A=[dots(R2crR3,NT)+dots(S2,GV)...  
   dots(R3crR1,NT)+dots(S3,GV)...  
   dots(R1crR2,NT)+dots(S1,GV)];
A=A./(asq*ones(1,3));
SA(:,k)=A';  

end
 