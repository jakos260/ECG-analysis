% nodearea.m
% function area=nodearea(VER,ITRI)
% 000327
% a van oosterom
% compute the local area represented by  the nodes of
%         a triangulated surface
% ITRI specifies three vertexes for each triangle
% VER  are the vertex coordinates
% see also: triareas; trinormals

function area=nodearea(VER,ITRI)
dim=size(VER);
nver=dim(1);
ntri=size(ITRI,1);

% compute, for all vertices, 1/3 of the surface area of the triangles around it
area=zeros(nver,1);
for i=1:ntri,
           i1=ITRI(i,1);
           i2=ITRI(i,2);
           i3=ITRI(i,3);
           
           rm=VER(i2,:)-VER(i1,:);
           rp=VER(i3,:)-VER(i1,:);
           ov=cross(rm,rp);
           opp=norm(ov)/6;
           area(i1)=area(i1)+opp;
           area(i2)=area(i2)+opp;
           area(i3)=area(i3)+opp;
end
