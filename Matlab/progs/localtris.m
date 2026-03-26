% localtris.m
% find all triangles that have at least one node within a radius rho from node
% function tris=localtris(VER,ITRI,node,rho);
% 20050314

function tris=localtris(VER,ITRI,node,rho);
  nver=size(VER,1);
  dsq=[norm3d(VER-ones(nver,1)*VER(node,:)) (1:nver)'];
  nodes=dsq(dsq(:,1)<=rho,2);
  nnodes=length(nodes);
  tris=[];
  for i=1:nnodes, 
     [bver,BTRI]=buren(ITRI,nodes(i));
     tris=[tris;BTRI(:,1)];
  end
tris=skipreplica(tris);
