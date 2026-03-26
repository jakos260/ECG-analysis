% tri_re_arrange.m
% 20141015
% [VER,ITRI,oldlabs]=tri_relabel(VERIN,ITRIN,seed);
% find new labels (1:size(unique(ITRI),1) for the nodes based on their proximity to node: "seed"
% VER and ITRI are sorted accordingly
% oldlabs are the vertex indices of VERIN;
% NB: size(VER,1) will be size(unique(ITRI),1)
% previously called: tri_relabel

% brute force version; preferably to be perfected later by spiral type

function [VER,ITRI,oldlabs]=tri_relabel(VERIN,ITRIN,seed);

nnodes=size(unique(ITRIN(:)),1)

sortednodes(1)=seed;
k=1;

while size(sortednodes,1) < nnodes,
    
    bver=buren(ITRIN,seed);
    bver(ismember(bver,sortednodes))=[];
    if isempty(bver)==0,
        sortednodes=[sortednodes;bver];
    end
    k=k+1;
    seed=sortednodes(k);
    %sortednodes
    %pause
end

LABS=[(1:nnodes)' sortednodes];
LABS=sortrows(LABS,2);

ntri=size(ITRIN,1);
ITRI=zeros(ntri,3);
%pause

for i=1:ntri,
    ITRI(i,1:3)=[LABS(LABS(:,2)==ones(nnodes,1)*ITRIN(i,1),1) ...
                 LABS(LABS(:,2)==ones(nnodes,1)*ITRIN(i,2),1) ...
                 LABS(LABS(:,2)==ones(nnodes,1)*ITRIN(i,3),1) ];   
end

% now sort triangles
sortedtris=[];
for i=1:nnodes
    [bver,BTRI]=buren(ITRI,i);
    newtris=BTRI(:,1);
    newtris(ismember(newtris,sortedtris))=[];
    if isempty(newtris)==0,
        sortedtris=[sortedtris;newtris];
    end
end
ITRI=ITRI(sortedtris,:);    
VER=VERIN(sortednodes,:);
oldlabs=sortednodes;




