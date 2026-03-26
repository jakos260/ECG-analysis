% edgedist.m
% dist=edgedist(nver,ITRI,node)
% edgdist finds de distances between node in a triangulated
% surface and all other nodes; 
% the distances are taken to be the minimum number of edges that have to
%   be passed before reaching the other nodes starting from %node: node.
%   the surface need not be closed; infinite values denote %%non-connected  
% date: 050416; older version !!
%!! use the more recent edge_dist(ITRI,node).m 
% graphdist should be upgraded to 

function dist=edgedist(nver,ITRI,node)
   ntri=size(ITRI,1);
   dist(1:nver)=inf;
   dist(node)=0;
   oldn=node;
   steps=0;
   TRI=ITRI;
   while isempty(TRI)==0,
       % find any triangles carrying oldn
       newn=[]; index=[];
       nold=length(oldn);
       for i=1:nold,
           BTRI=[];buur=[];
           [buur,BTRI]=buren(TRI,oldn(i));
           newn=[newn buur'];
           index=[index; BTRI(:,1)];
       end
       index=unique(index);
       newn=unique(newn);
       state=dist(newn);
       
       
       if isempty(newn)==1, break, end
       newn(sum([newn; state])<inf)=[];
       
       if isempty(newn)==1, break, end
       index=skipreplica(index);
       steps=steps+1;
       dist(newn)=steps;
       oldn=newn;
       TRI(index,:)=[];
    end
   
   
       
        