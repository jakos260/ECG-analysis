% trisense.m
% function [sense]=trisense(ITRI,VER)
% determine sense of the triangles of a surface
% 010714
function [sense]=trisense(ITRI,VER)
   R1(:,1:3)=VER(ITRI(:,2),1:3)-VER(ITRI(:,1),1:3);
   R2(:,1:3)=VER(ITRI(:,3),1:3)-VER(ITRI(:,1),1:3);
   R3=cross(R1,R2);
   sense=sign(R3(:,1));

