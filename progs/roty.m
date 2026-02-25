% roty.m
% function variant of the script rot_y
% rotate VER about its center of gravity by roty DEGREES around positive y-axis; positive: anti-clockwise
% total of successive rotations may be computed by ROTS=ROT*ROTS
% 20090311

function [VEROUT,ROMA]=roty(VER,roty)

roty=roty/180*pi;
nver=size(VER,1);
vermean=mean(VER);
SHIFT=ones(nver,1)*vermean;
VER_zm=VER-SHIFT;
ROMA=[ cos(roty) 0 -sin(roty); 0 1 0; sin(roty) 0 cos(roty)];
VEROUT=(ROMA*VER_zm')'+SHIFT;
