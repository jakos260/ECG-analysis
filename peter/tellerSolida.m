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

function  [sa]=tellerSolida(VER,ITRI,obs)
  
[ntri ~]=size(ITRI);
NORMT=zeros(ntri,3);
MIDTRI=zeros(ntri,3);
dots= zeros(ntri,1);
for i=1:ntri,
    i1=ITRI(i,1);
    i2=ITRI(i,2);
    i3=ITRI(i,3);
    rm=VER(i2,:)-VER(i1,:);
    rp=VER(i3,:)-VER(i1,:);
    NORMT(i,:)=cross(rp,rm)/2;
    MIDTRI(i,:) = mean(VER([i1 i2,i3],:))-obs;
    dots(i)=dot(NORMT(i,:),MIDTRI(i,:));
end

sa = sum(dots);


