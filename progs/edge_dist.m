% edge_dist.m
% D=edge_dist(ITRI,node)
% edgdist finds de distances D(:,2) between any node: node of a triangulated
% surface ITRI and all other nodes
% the distances are taken to be the minimum number of edges that have to
% be passed before reaching the other nodes starting from node: node.
% the surface need not be closed
% date: 090313; version; permitting ITRI labels that are a subset of the vertices.
% % see edgedist(nver,ITRI,node) for an older version, still used in
% some other scripts
 

function D=edge_dist(ITRI,node)

   nodes=unique(ITRI);
   nnodes=size(nodes,1);
   
   D=[nodes ones(nnodes,1)*inf];    % NB: non-connecting nodes will be assigned the value inf.
   D(nodes==ones(nnodes,1)*node,2)=0;
   
   oldn=node;
   nsteps=0;
   TRI=ITRI;
   
   while isempty(TRI)==0,
       % find all ITRIs carrying an oldn
       newn=[]; itris=[];
       nold=length(oldn);
       for i=1:nold,
           BTRI=[];buur=[];
           [buur,BTRI]=buren(TRI,oldn(i));
           if isempty(buur)==0,
              newn=[newn buur'];
              itris=[itris; BTRI(:,1)]; % indices of triangles carying the node
           end
       end
       
       if isempty(newn)==1, break, end
       itris=unique(itris);
       newn=unique(newn);
       newn(ismember(newn,oldn))=[];
       nsteps=nsteps+1;
       D(ismember(D(:,1),newn'),2)= nsteps;
       oldn=newn;
       TRI(itris,:)=[];
    end
   
   
       
        