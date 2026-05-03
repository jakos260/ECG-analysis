% skipneighbours.m
% SINGLES=skipneighbours(LIST,near)
% skip neighbours (rows) that are too close for comfort

% A. van Oosterom; 2006/01/23
function SINGLES=skipneighbours(LIST,near)
     SINGLES=[];
     while isempty(LIST)==0,
           nlist=size(LIST,1);
           single=LIST(1,:);
           SINGLES=[SINGLES;single];
           d=sum((LIST-ones(nlist,1)*single).^2,2);
           list=find(d<=near);
           LIST(list,:)=[];
     end
     SINGLES=[SINGLES;SINGLES(1,:)];

