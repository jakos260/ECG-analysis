% file getvolume.m
% function vol=getbolume(VER,ITRI)
function vol=getvolume(VER,ITRI)
[ntri,jdum]=size(ITRI);
d=det3d(VER(ITRI(:,1),:),VER(ITRI(:,2),:),VER(ITRI(:,3),:));
vol=-sum(d)/6;


