% file trivolume.m
% function vol=trivolume(VER,ITRI)
% compute the volume within a closed, triangulated surface;
% demands consistent orientation of triangle normals !!!!
% sign of output depends on the orientation (outward or inward)
% A. van Oosterom; 202_12_20
function vol=trivolume(VER,ITRI)
[ntri,jdum]=size(ITRI);
d=det3d(VER(ITRI(:,1),:),VER(ITRI(:,2),:),VER(ITRI(:,3),:));
vol=-sum(d)/6;


