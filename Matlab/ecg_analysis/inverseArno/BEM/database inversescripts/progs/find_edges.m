% find_edges.m
% ALLEDGES=find_edges(ITRI,VER)
% ALLEDGES(:,1:2): edge indices;  ALLEDGES(:,3): triangle index;

% date: 2012313

function ALLEDGES=find_edges(ITRI)

   ntri=size(ITRI,1);
   ALLEDGES=[ITRI(:,1:2) (1:ntri)'; ITRI(:,2:3) (1:ntri)'; ITRI(:,[3 1]) (1:ntri)'];
   ALLEDGES=sortrows(ALLEDGES,1);
 
%    EDGES=[];
%    while isempty(ALLEDGES)==0,
%        nall=size(ALLEDGES,1);
%        k=find(ALLEDGES(:,1)==ones(nall,1)*ALLEDGES(1,2) & ALLEDGES(:,2)==ones(nall,1)*ALLEDGES(1,1));
%        EDGES=[EDGES;[min(ALLEDGES(1,:)) max(ALLEDGES(1,:))] ];
%        ALLEDGES([1; k],:)=[];    
%    end
%  
% if nargin>1,
%    D=VER(EDGES(:,1),:)-VER(EDGES(:,2),:);
%    el=[ALLEDGES norm3d(D)];
% end
       
   
        