% buurtris.m
% function [btri,BVER]=buurtris(ITRI,tri)
% for triangle  tri of a triangulate surface ITRI: find 
%               btri: the indexes of triangles sharing an edge with tri
%               BVER: labels the vertex indices of tri shared with btri
%               example: BVER(i,1:4)=[3 1 2 1] indicates that 
%                                       vertex 3 of tri = vertex 2 of btri(i)
%                                and    vertex 1 of tri = vertex 1 of btri(i)

% for correctly triangulated surfaces the cyclic order of BVER(i,1:2) and
%                                                 that of BVER(i,3:4) 
%                                                 should be opposite
% see also: buren,voisin, findedge, findloop

% A. van Oosterom; 20050201

function [btri, BVER]=buurtris(ITRI,tri)

  btri=[]; BVER=[];
  ntri=size(ITRI,1);
  itri=ITRI(tri,:);
  ITRIS=[(1:ntri)' ITRI];
  index=sum(ITRIS(:,2:4)==itri(1)|ITRIS(:,2:4)==itri(2)|ITRIS(:,2:4)==itri(3),2)==2;
  if isempty('index')==1; return, end
  
  btri=ITRIS(index,1);
  TEST=ITRI(btri,:);
  TEST1=(TEST==itri(1));
  TEST2=(TEST==itri(2))*2;
  TEST3=(TEST==itri(3))*3;
  ALL=TEST1+TEST2+TEST3;
  nall=size(ALL,1);
  
  test=1:3;
  for i=1:nall,
      BVER(i,3:4)=[test(ALL(i,:)==1) test(ALL(i,:)==2) test(ALL(i,:)==3)];
      tmp=ALL(i,ALL(i,:)~=0);
      BVER(i,1:2)=tmp(tmp~=0);
  end
  
  
