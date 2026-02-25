% solida.m
% function  [sa,index]=solida(VER,ITRI,obs)
% computes the row vector representing the solid angles sa at observation point obs 
% subtended by all triangles (VER,ITRI): vertices: VER(j,[xj yj zj]); ITRI(i,[ki li mi]), the vertex indices
% as well as a row vector index: 
%   if index(i)==0: obs does not lie in the plane of triangle i
%      else obs lies in the plane of the triangle and 
%         if: index(i)=1: obs is an external point of triangle i 
%             index(i)=2: obs lies on an edge      of triangle i
%             index(i)=3: is an internal point     of triangle i
%             index(i)=4: obs is a vertex          of triangle i
%             BASIC version: (IEEE,BME_30, 1983, pp. 125-126); 
% for distributed solid angles: use dsa.m

% A. van Oosterom; 2006-05-02 ;

function  [sa,index]=solida(VER,ITRI,obs)
[nver,~]=size(VER);
VER=VER-ones(nver,1)*obs;
R1=VER(ITRI(:,1),1:3);
R2=VER(ITRI(:,2),1:3);
R3=VER(ITRI(:,3),1:3);
LR=[ norm3d(R1) norm3d(R2) norm3d(R3)];
R2crR3=cross(R2,R3);
ntri=size(ITRI,1);
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
sa(abs(block)<=eps)=0;
