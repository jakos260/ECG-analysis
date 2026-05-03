% find_clusters.m
% function clus=find_clusers(A)
% find clusters of connected nodes in a graph from its adjacency matrix A
% A. van Oosterom; 20050412

function CLUS=find_clusters(A)

CLUS=[];
nn=size(A,1);
nclus=0;
placed=zeros(1,nn);
clus=[];

while sum(placed)<nn,
    nclus=nclus+1;
    k=find(placed==0);
    testrow=k(1);
    
    while isempty(testrow)~=1,
        row=testrow(1);
        A(row,:);
        clus=[clus row];
        testrow(1)=[];
        add=find(A(row,:)~=0);
        if isempty(add)==0, 
            add(ismember(add,clus))=[];
            if isempty(add)==0, 
                testrow=[testrow add];
                testrow=unique(testrow);
            end
        end
    end
    
    placed(clus)=1;
    clus=unique(clus);
    [ni nj]=size(CLUS);
    nc=length(clus);
    if nc>ni, CLUS=[CLUS; zeros(nc-ni,nj)]; end
    if nc<ni, clus=[clus zeros(1,ni-nc)];  end  
    nc=length(clus);
    CLUS=[CLUS clus'];
    clus=[];
end


    
