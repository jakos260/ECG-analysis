% zoomin_tris.m
% function zoomin_tris(VERINP,ITRIINP,node,rho);
% zoom in on the node of a triangulated surface
% all triangles having all nodes within distance: rho from nodes are shown

function list=zoomin_tris(VERINP,ITRIINP,node,rho);

[ntri jdum]=size(ITRIINP);
[nver jdum]=size(VERINP);

dist=norm3d(VERINP-ones(nver,1)*VERINP(node,:));
% list of nodes in view:
list=find(dist<=rho)
ITRI=[];
for i=1:ntri,
    if any(list==ITRIINP(i,1)) & any(list==ITRIINP(i,2)) & any(list==ITRIINP(i,3)),
        ITRI=[ITRI; ITRIINP(i,:)];
    end
end
ITRI

VER=VERINP;

tri_nodes=unique(ITRI(:));
meanv=mean(VER(tri_nodes,:))

VER=VER-ones(size(VER,1),1)*meanv;

figure(1)
clf
VALS=VER;
triplot


% figure(2)
% clf
% crossec
% 
% figure(1)