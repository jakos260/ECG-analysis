% dis2tri.m (==lamutri.m)
% function [la,mu,dist]=dis2tri(p,VER,ITRI);
function [la,mu,dist]=dis2tri(p,VER,ITRI);

% find  lambda and mu parameters          
%       specifying the projection(s) of the point p (rowvector) on the plane(s)    
%       of the triangle(s) ITRI
% dist is the distance to the plane of the triangle 	
%      its sign is 
%      positive if p lies in the semi-space   
%               into which the normal (R2-R1)cros(R3-R1) is pointing  
%      zero     if p lies in the plane of the triangle; 
%      negative otherwise   
% projection of p on tri is: (1-la-mu)*VER(1,:)+ la*VER(2,:)+mu*VER(3,:) 

% 20080826
[nver,jdum]=size(VER);
[ntri,jdum]=size(ITRI);
la=zeros(ntri,1);
mu=zeros(ntri,1);
dist=zeros(ntri,1);

       S1   =VER(ITRI(:,2),1:3)-VER(ITRI(:,1),1:3); % edge(s) 1 of the triangle(s)
       S2   =VER(ITRI(:,3),1:3)-VER(ITRI(:,1),1:3); % edge(s) 3 of the triangle(s)
       B    =ones(ntri,1)*p(1:3)-VER(ITRI(:,1),1:3);% vector(s) pointing to obs
       N    =cross(S1,S2);
       rn   =max(ones(ntri,1)*eps,norm3d(N));
       deter=det3d(S1,S2,N)+ones(ntri,1)*eps;
       la   =det3d(B,S2,N)./deter;
	   mu   =det3d(S1,B,N)./deter;
       dist =det3d(S1,S2,B)./rn;

