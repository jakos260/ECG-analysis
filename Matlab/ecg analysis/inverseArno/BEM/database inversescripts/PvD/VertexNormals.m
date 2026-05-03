function normals = VertexNormals(VER,ITRI)

%  Computes the normal at a given vertex by avarging the neighboring face
%  normals

normals = VER;
for i=1:length(VER)
    [ti,l]=find(ITRI == i);
    normal = cross(VER(ITRI(ti, 3),:)-VER(ITRI(ti, 1),:),...
                   VER(ITRI(ti, 2),:)-VER(ITRI(ti, 1),:));
    normal = sum(normal./(norm3d(normal)*ones(1,3)));
    normals(i,:) = normal./norm(normal); 
    if isnan(normals(i,:))
        stip=1;
    end
end