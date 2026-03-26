% nodedist.m
% DIST=nodedist(nver,ITRI,node)
% nodedist finds de distances between node in a triangulated
% surface and all other nodes of ITRI; the distances are taken to be the minimum
% number of edges that have to
% be passed before reaching the other nodes starting from node: node.
% the surface need not be closed
% see also edgedist
% date: 050412

%%%%         under construction; see edgedist



function DIST=nodedist(nver,ITRI,node)
   
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
       index=skipreplica(index);
       newn=skipreplica(newn);
       state=dist(newn);
       newn(sum([newn; state])<inf)=[];
       
       if isempty(newn)==1, break, end
       index=skipreplica(index);
       steps=steps+1;
       dist(newn)=steps;
       TRI(index,:)=[];
    end
   
   
       
        