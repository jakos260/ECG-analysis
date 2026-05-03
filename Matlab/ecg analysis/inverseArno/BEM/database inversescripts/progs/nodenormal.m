% nodenormal.m
% function normal=nodenormal(VER,ITRI,node)
% computes the unit_outward normal of a node: node on a triangulated  surface
% from the mean of the norms of the triangles on which node is situated
% ITRI specifies three vertexes for each triangle
% VER  are the vertex coordinates
% see also: nodearea; triarea; trinormals
% 20090409; a van oosterom


function normal=nodenormal(VER,ITRI,node)
[bver,BTRI]=buren(ITRI,node);
nbtri=size(BTRI,1);
nbver=length(bver);
normal=[0 0 0];
if nbtri~=nbver;
    'incomplete triangulation around node';
    return
end
for i=1:nbtri,
    i1=BTRI(i,2);
    i2=BTRI(i,3);
    i3=BTRI(i,4);
    rm=VER(i2,:)-VER(i1,:);
    rp=VER(i3,:)-VER(i1,:);
    normal=normal+cross(rm,rp);
end
normal=-normal/norm(normal);