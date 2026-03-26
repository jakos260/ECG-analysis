function normal = VertexNormal(VER,ITRI,vertexIndex)

%  Computes the normal at a given vertex by avarging the neighboring face
%  normals
[ti,l]=find(ITRI == vertexIndex);
normal=cross(VER(ITRI(ti(:),3),:)-VER(ITRI(ti(:),1),:),VER(ITRI(ti(:),2),:)-VER(ITRI(ti(:),1),:));
normal = sum(normal./(norm3d(normal)*ones(1,3)));
% normal = mean(b); 
normal = normal./norm(normal); 