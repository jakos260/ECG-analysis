% lamutri.m
% function [la,mu,dist]=LAMUTRI(p,VER,ITRI);
function [la,mu,dist]=LAMUTRI(p,VER,ITRI);
% find  lambda and mu parameters          
%       specifying the projections of the point p (rowvector) on the plane(s)    
%       of the triangle(s) ITRI
% dist is the distance to the plane of the triangle 	
%      its sign is positive if p lies in the semi-space   
%                           into which the normal (R2-R1)cros(R3-R1) is pointing  
%      zero indicates that p lies in the plane of the    
%      triangle; else dist is negative                  

% 2003-04-18
[nver,jdum]=size(VER);
[ntri,jdum]=size(ITRI);
la=zeros(ntri,1);
mu=zeros(ntri,1);
dist=zeros(ntri,1);

       S1   =VER(ITRI(:,2),1:3)-VER(ITRI(:,1),1:3);
       S2   =VER(ITRI(:,3),1:3)-VER(ITRI(:,1),1:3);
       B    =ones(ntri,1)*p(1:3)-VER(ITRI(:,1),1:3);
       N    =cross(S1,S2);
       rn   =norm3d(N);
       deter=det3d(S1,S2,N);
	   la   =det3d(B,S2,N)./deter;
	   mu   =det3d(S1,B,N)./deter;
	   dist =det3d(S1,S2,B)./rn;
