% trinormals.m
% function [NORMV, NORMT]=trinormals(VER,ITRI)
% NORMT: the outward normals (computed from the cross products of two edges)
%                of a set of triangles
% NORMV: the normals at VER are computed as the sum of 1/3 the normals of the
%        triangles that carry VER; the norm of element i is the
%        area represented by node i 
% ITRI specifies three vertices for each triangle
% VER  specifies the vertex coordinates
% NB: to compute areas around each node and/or areas of individual triangles
%     use: triareas

% 20050523 a. van oosterom


function [NORMV,NORMT]=trinormals(VER,ITRI)
   [ntri, ~]=size(ITRI);
   [nver, ~]=size(VER);
   NORMT=zeros(ntri,3);
   NORMV=zeros(nver,3);
   
for i=1:ntri
    i1=ITRI(i,1);
    i2=ITRI(i,2);
    i3=ITRI(i,3);
    rm=VER(i2,:)-VER(i1,:);
    rp=VER(i3,:)-VER(i1,:);
    
    NORMT(i,:)=cross(rp,rm)/2;
    NORMV([i1 i2 i3],:) = NORMV([i1 i2 i3],:) + ones(3,1)*NORMT(i,:)/3;
end
NORMV = bsxfun(@times,NORMV, 1./ norm3d(NORMV));
NORMT = bsxfun(@times,NORMT, 1./ norm3d(NORMT));
