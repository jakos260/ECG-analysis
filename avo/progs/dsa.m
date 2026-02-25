% dsa.m
% function  [SA,index]=dsa(VER,ITRI,obs,small)
% compute the matrix of the distributed solid angles SA at observation point obs (size:  3, ntri)
% subtended by a set of triangles defined by (VER,ITRI)
% A. van Oosterom; 2015-07-12 
% if index==0: obs lies not in the plane of the triangle
   % if any solid angle (sa) is smaller than small, one third of the
   % of the sa is assigned to every vertex of the triangle
   % basic version (IEEE,BME_30, 1983, pp. 125-126),
   % else, a further computation yields
   % the values according to de Munck: IEEE Trans BME_39,1993, pp986-990; eqn 19.
% else SA(*,1:3)=0; obs lies in the plane of the triangle;
% moreover:
     % index=1: obs is an external point of the triangle
     % index=2: obs lies on the edge of the triangle;
     % index=3: obs is an internal point of the triangle
     % index=4: obs is a vertex of the triangle  
    
% assumes the outward normals of the triangles (of a closed surface) to be defined by a clockwize rotation of the vertices ITRI(:,[1 2 3]). 
% For observation points (obs) inside a CLOSED surface the SUM of all SA
% values (in STAR_RADIANS) of all triangles
% is 4 pi; for external obs SUM is zero; 
% 
% Calling: norm3d; cross; dots; norm3d

function  [SA,index]=dsa(VER,ITRI,obs,small)
[nver,~]=size(VER);
VER=VER-ones(nver,1)*obs;

R1=VER(ITRI(:,1),:);
R2=VER(ITRI(:,2),:);
R3=VER(ITRI(:,3),:);
LR=[ norm3d(R1) norm3d(R2) norm3d(R3)];
R2xR3=cross(R2,R3);
ntri=size(ITRI,1);
index=zeros(1,ntri);
% blockproducts
block=dots(R1,R2xR3);
DOTS=[dots(R1,R2) dots(R2,R3) dots(R3,R1)];
denom=LR(:,1).*LR(:,2).*LR(:,3) + ...
      LR(:,1).*DOTS(:,2)+LR(:,2).*DOTS(:,3) + LR(:,3).*DOTS(:,1);
index(abs(block)<=eps & denom >   eps)=1;
index(abs(block)<=eps & denom >= -eps & denom <=eps)=2;
index(abs(block)<=eps & denom <  -eps)=3;
index(LR(:,1)<=eps|LR(:,2)<=eps|LR(:,3)<=eps)=4;
sa=-2*atan2(block,denom)';

sa(abs(block)<=eps)=0; % crude approximation
SA  =ones(3,1)*sa/3;
TEST=[1:ntri;sa];
k   =TEST(1,abs(TEST(2,:))>small);

if isempty(k)==0
   % dsa computation,  required for triangles: k
   S1=R2(k,:)-R1(k,:);
   S2=R3(k,:)-R2(k,:);
   S3=R1(k,:)-R3(k,:);
   LS=[norm3d(S1) norm3d(S2) norm3d(S3)];
   DOTS=DOTS(k,:);

   R2xR3=R2xR3(k,:);
   R3xR1=cross(R3(k,:),R1(k,:));
   R1xR2=cross(R1(k,:),R2(k,:));
   
   NT    =R1xR2+R2xR3+R3xR1;  % dimension m^2
   asq   =dots(NT,NT);  % dimension m^4
   NT    =(sa(k)'*ones(1,3)).*NT;

   % compute the gamma vector
   m=[2 3 1];
   
   NUM = LR(k,:).*LS+DOTS-LR(k,:).*LR(k,:);
   DNOM= LR(k,m).*LS-DOTS+LR(k,m).*LR(k,m);
 
   GAM = ((block(k)*ones(1,3)).*real(log(NUM./DNOM)))./LS;

   GV(:,1)=S1(:,1).*GAM(:,1) + S2(:,1).*GAM(:,2) + S3(:,1).*GAM(:,3);
   GV(:,2)=S1(:,2).*GAM(:,1) + S2(:,2).*GAM(:,2) + S3(:,2).*GAM(:,3);
   GV(:,3)=S1(:,3).*GAM(:,1) + S2(:,3).*GAM(:,2) + S3(:,3).*GAM(:,3);
 
   A =[dots(R2xR3,NT) + dots(S2,GV)...  
       dots(R3xR1,NT) + dots(S3,GV)...  
       dots(R1xR2,NT) + dots(S1,GV)];
   A =A./(asq*ones(1,3));
   SA(:,k)=A';   
end
 

%    [dots(R2xR3,NT); % dimension: m^4
%    dots(R3xR1,NT);  
%    dots(R1xR2,NT)]
 
%  [dots(S2,GV); % dimension: m^4
%   dots(S3,GV); 
%   dots(S1,GV)]
%   dimension asq: m^4
    
