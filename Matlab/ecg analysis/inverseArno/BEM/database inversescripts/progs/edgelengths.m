% edgelengths.m
% EDGES=edgelengths(VER,ITRI)
% edgelengths finds all edges contained in ITRI and their length
% EL(:,1:2) edge indices;  EL(:,3) lengths

% date: 20110427

function EDGES=edgelengths(VER,ITRI)
   ntri=size(ITRI,1);
   ALLTRIS=[ITRI(:,1:2); ITRI(:,2:3); ITRI(:,[3 1])];
   ALLTRIS=sortrows(ALLTRIS,1);
   
   
   EDGES=[];
   while isempty(ALLTRIS)==0,
       nall=size(ALLTRIS,1);
       k=find(ALLTRIS(:,1)==ones(nall,1)*ALLTRIS(1,2) & ALLTRIS(:,2)==ones(nall,1)*ALLTRIS(1,1));
       EDGES=[EDGES;[min(ALLTRIS(1,:)) max(ALLTRIS(1,:))]];
       ALLTRIS([1; k],:)=[];    
   end
 
   D=VER(EDGES(:,1),:)-VER(EDGES(:,2),:);
   EDGES=[EDGES norm3d(D)];
       
   
        