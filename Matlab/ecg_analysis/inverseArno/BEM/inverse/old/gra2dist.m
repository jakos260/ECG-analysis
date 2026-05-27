% function D = gra2dist(G)
%
% Compute distances (shortest path) between all nodes of a triangular 
% surface based on G, the weighted adjacency matrix of a the involved 
% graph (computed, e.g. by tri2gra).
%
% 2004-oct-27
%
% Version 1: 01-apr-2015

function D = gra2dist(G)

dim = size(G);
n   = dim(1);
D   = zeros(n,n);

for node = 1:n,
    fix         = zeros(1,n);
    fix(node)   = 1;
    nfix        = 1;
    
    dist        = G(node,:);
    % dist(j) contains the distances from node:node to nodes(j) j=1:n
    
    while nfix < n,
        % find smallest dist value of a free node
        ind     = find(dist & ~fix);
        [y,id]  = min(dist(ind));
        
        % fix it
        newfix      = ind(id(1));
        fix(newfix) = 1;
        nfix        = nfix+1;
        
        if nfix >= n, break, end
        
        % set new values at free neighbours of newfix
        row     = G(newfix,:);
        ina     = find(row & ~fix);
        distnow = dist(ina);
        
        % suggested new distance values
        newdist = row(ina)+y(1);
        
        % accept newdist(ina) only if this leads to smalled dist values
        dist(ina) = min(dist(ina)+~dist(ina).*newdist,newdist);
    end
    
    D(node,:) = dist;
end