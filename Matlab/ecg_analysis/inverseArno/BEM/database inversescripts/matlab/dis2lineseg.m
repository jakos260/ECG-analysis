% dis2lineseg.m
% function [dist,lambda]=dis2line(NODES,lp1,lp2)
% determine distances (column vector) of nn NODES (size=nn,3) to
% the line segment running from rowvector lp1 to rowvector lp2
% lambda is the fraction: (projection-lp1)/(lp2-lp1)
% 20120611 oostep1. Adapted from dis2line.m


function [dist,lambda]=dis2line(NODES,lp1,lp2)
dim=size(NODES);
nn=dim(1);
a=norm(lp2-lp1);
V1=ones(nn,1)*(lp2-lp1)/a;
V2=NODES-ones(nn,1)*lp1;
dots=V1(:,1).*V2(:,1)+V1(:,2).*V2(:,2)+V1(:,3).*V2(:,3);
D=dots*ones(1,3).*V1-V2;
dist=norm3d(D);
lambda=dots/a;

% check if nearest point is on line segment
% TBD: vectorize

for i=1:size(NODES,1)
    if lambda(i)<0
        lambda(i)=0;
        dist(i)=norm(NODES(i,:)-lp1);
    elseif lambda(i)>1
        lambda(i)=1;
        dist(i)=norm(NODES(i,:)-lp2);
    end
end