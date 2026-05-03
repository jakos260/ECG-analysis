% triangles.m
% function TRIANG=triangles(VER,ITRI)
% compute the angles at the vertices of the planar triangles 
% ITRI specifies three vertices for each triangle
% VER specifies the vertex coordinates
% see also trinormals

% 050201 a van oosterom

function TRIANG=triangles(VER,ITRI)
 
   R1=VER(ITRI(:,2),:)-VER(ITRI(:,1),:);
   R2=VER(ITRI(:,3),:)-VER(ITRI(:,2),:);
   R3=VER(ITRI(:,3),:)-VER(ITRI(:,1),:);
   SHARP=sign([sum(R1.*R3,2) sum(-R2.*R1,2) sum(R2.*R3,2)]) 
   LR=[ norm3d(R1) norm3d(R2) norm3d(R3)];
   R1crR3=norm3d(cross(R1,R3));
   SINS=[R1crR3./(LR(:,1).*LR(:,3)) R1crR3./(LR(:,2).*LR(:,1)) R1crR3./(LR(:,2).*LR(:,3))];
   SINS(SINS>1)=1; SINS(SINS<-1)=-1;
   TRIANG=asin(SINS);
   TRIANG(SHARP<0)=pi-TRIANG(SHARP<0);


