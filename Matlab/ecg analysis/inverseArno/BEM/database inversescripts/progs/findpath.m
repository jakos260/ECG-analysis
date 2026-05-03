% findpath.m
% function [path,dist]=findpath(G,i,j)
% computes the route derived from shortest path between nodes i and j of a
% triangulated surface
% based on the precomputed G, the (weighted) adjacency matrix of the related graph
% computed by means of , e.g., tri2gra.m;
% path: the route;   dist: the distances along the route
% 20050110; A. van Oosterom

function [path,dist]=findpath(G,i,j)

[d p]=distgra1(G,i);
if i==j, p(i)=i; p(j)=i; end
pj=p(j);
path=[p(j) j];

while pj~=i,
  pj=p(pj);
  path=[pj path];
end
dist=d(path);
