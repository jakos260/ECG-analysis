% edgelengths.m
% EDGES=edgelengths(VER,ITRI)
% edgelengths finds all edges contained in ITRI and their length
% EL(:,1:2) edge node indices;  EL(:,3) lengths

% date: 20140422   

function EDGES=edgelengths(VER,ITRI)
   ntri=size(ITRI,1);
   ALLEDGES=[ITRI(:,1:2) ; ITRI(:,2:3);ITRI(:,[3 1]) ]; 
   EDGES=[];
   while isempty(ALLEDGES)==0,
       n=size(ALLEDGES,1);
       k=find(ALLEDGES(2:n,2)==ones(n-1,1)*ALLEDGES(1,1) & ALLEDGES(2:n,1)==ones(n-1,1)*ALLEDGES(1,2));

       EDGES=[EDGES;ALLEDGES(1,1:2)];
       ALLEDGES([1 k+1],:)=[];     
   end
   
   D=norm3d(VER(EDGES(:,1),:)-VER(EDGES(:,2),:));
   EDGES=[EDGES D];
 
       
   
    