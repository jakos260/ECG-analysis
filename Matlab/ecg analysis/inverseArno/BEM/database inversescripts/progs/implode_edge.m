% implode_edge.m
% function [VER,ITRI]=implode_edge(VER,ITRI,edge)
% remove an edge running from nodes: edge(1) to edge(2)
% adapt mesh; introduces new vertex, halfway in between n1 and n2 
% old nodes may be removed later by means of tri_clean
% removes 2 triangles; NB; usefull mainly for removing very short edges
% A. van Oosterom; 20090329


function [VER,ITRI]=implode_edge(VER,ITRI,edge)
VER=[VER; mean(VER(edge,:))];
nver=size(VER,1);
ITRI(sum(ismember(ITRI,edge),2)>1,:)=[];
k=find(sum(ismember(ITRI,edge),2)==1);
ITRI(ismember(ITRI,edge))=nver;


  
  




