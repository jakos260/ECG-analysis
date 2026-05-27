% edge_usage.m
% function EDGES=edge_usage(ITRI)
% 20090415
% for all edges of ITRI: list the edges, their triangle and their match, the latter zero if absent; 
% EGES(: edge itri matching_itri); in a (correctly meshed) closed geometry:  nedges=3/2*ntri
% see also tricheck;  trisect


function EDGES=edge_usage(ITRI)
ntri=size(ITRI,1);

ind=[1:ntri 1:ntri 1:ntri]; 
EDGES=[ITRI(:,1) ITRI(:,2); ITRI(:,2) ITRI(:,3); ITRI(:,3) ITRI(:,1)];
EDGES=[EDGES ind' zeros(3*ntri,1)];
nedg=3*ntri;
  
for i=1:nedg,
    if EDGES(i,4)==0,
       % find its paired version
       k=find(EDGES(:,1)==ones(nedg,1)*EDGES(i,2) & EDGES(:,2)==ones(nedg,1)*EDGES(i,1));
       if max(size(k))>1,
          'warning: edge usage greater than 2'
          i
          k  
          EDGES([i;k],:)
          %pause
       end
       
       if isempty(k)==0;
          EDGES(i,4)=EDGES(k(1),3);
          EDGES(k(1),4)=EDGES(i,3);
       end
       
    end
end

    
 