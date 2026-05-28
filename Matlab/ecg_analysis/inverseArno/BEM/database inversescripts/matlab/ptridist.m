% ptridist.m 
% function [la,mu,dist]=ptridist(p,VER,ITRI);
% find  lambda and mu parameters          
%       specifying the projection(s) of the point p (rowvector) on the triangle(s) ITRI  
% dist is the distance to the triangle 	
% projection of p on tri is: (1-la-mu)*VER(1,:)+ la*VER(2,:)+mu*VER(3,:) 
% ( 1,2 and 3 being ITRI(n,1), ITRI(n,2) and ITRI(n,3) resp. if size(ITRI,1)>1 )
%      positive if p lies in the semi-space   
%               into which the normal is pointing  
%      zero     if p lies within the triangle; 
%      negative otherwise   
% projection of p on tri is: (1-la-mu)*VER(1,:)+ la*VER(2,:)+mu*VER(3,:) 


%20120611 oostep1: adapted from ptridist in vecfun.cc

function [la,mu,dist]=ptridist(p,VER,ITRI);

[la,mu,dist]=dis2tri(p,VER,ITRI); % projection to planes of triangles in ITRI
distsign=sign(dist);
distsign(dist==0)=1; % in case dist is 0

% Handle when projected point is outside of triangles
for i=1:size(ITRI,1)
    if la(i)<0
        la(i)=0;
        [dist(i),mu(i)]=dis2lineseg(p,VER(ITRI(i,1),:),VER(ITRI(i,3),:));
        dist(i)=dist(i)*distsign(i);
    elseif mu(i)<0
        mu(i)=0;
        [dist(i),la(i)]=dis2lineseg(p,VER(ITRI(i,1),:),VER(ITRI(i,2),:));
        dist(i)=dist(i)*distsign(i);
    elseif (la(i)+mu(i))>1
        [dist(i),mu(i)]=dis2lineseg(p,VER(ITRI(i,2),:),VER(ITRI(i,3),:));
        la(i)=1-mu(i);
        dist(i)=dist(i)*distsign(i);
    end
end




