% distgra1.m
% function called by findroute
% function [d p]=distgra1(G,node);
% vector d  contains shortest path length from node to
% to all nodes along edges of a graph
% p are the parents of the nodes along the shortest route 
% G should contain pre-computed (weighted) direct internode distances  
%          G(i,i)=0; G(i,j)=0 if nodes are not connected
% graph G computed, e.g. by tri2gra

function [d, parents]=distgra1(G,node);
[n n]=size(G);
nnodes=1;
D=zeros(1,n);

fix=zeros(1,n);
parents=zeros(1,n);
newfix=node;
fix(newfix)=1;
parents(node)=0;
nfix=1;
dist=G(node,:);
parents(find(dist))=node;

% in each of the next iterations dist(j), j=1:n, is the 
% current estimate of the smallest distance from node to node j. 
% nb: dist==0 if, in fact,  dist=inf

while nfix <= n,
  % find smallest dist value of a free node
  [dist;fix];
  ind=find(dist & ~fix);
  [y,id]=min(dist(ind));

  % determine newfix; then fix it 
  newfix=ind(id(1));
  fix(newfix)=1;
  
  nfix=nfix+1;
  if nfix>= n, break, end

  % set new values at free neighbours of newfix
  
   row=G(newfix,:);

  ina=find(row & ~fix);
  distnow=dist(ina);

  % suggested new distance values
  newdist=row(ina)+y(1);
  % accept those newdist(ina) values that lead to smaller dist values
  dist(ina)=min(dist(ina)+~dist(ina).*newdist,newdist);
  newpar=ina(dist(ina) < distnow | distnow==0);
  if length(newpar) > 0, parents(newpar)=newfix; end
end
d=dist;
