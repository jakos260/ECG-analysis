% rotx.m
% function variant of the script rot_x
% rotate VER about its center of gravity by rotx DEGREES around positive x-axis; positive: anti-clockwise
% total of successive rotations may be computed by ROTS=ROT*ROTS
% NB: mean position is retained
% 20090311

function [VEROUT,ROMA]=rotx(VER,rotx)

rotx=rotx/180*pi;
nver=size(VER,1);
vermean=mean(VER);
SHIFT=ones(nver,1)*vermean;
VER_zm=VER-SHIFT;
ROMA=[1 0 0; 0 cos(rotx) sin(rotx); 0 -sin(rotx) cos(rotx)];
VEROUT=(ROMA*VER_zm')'+SHIFT;
