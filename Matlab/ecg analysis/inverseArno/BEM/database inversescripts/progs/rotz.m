% rotz.m
% function variant of the script rot_z
% rotate VER about its center of gravity by rotz DEGREES around positive z-axis; positive: anti-clockwise
% total of successive rotations may be computed by ROTS=ROT*ROTS
% 20090311

function [VEROUT,ROMA]=rotz(VER,rotz)

rotz=rotz/180*pi;
nver=size(VER,1);
vermean=mean(VER);
SHIFT=ones(nver,1)*vermean;
VER_zm=VER-SHIFT;
ROMA=[cos(rotz)  sin(rotz) 0;  -sin(rotz) cos(rotz) 0; 0 0 1];
VEROUT=(ROMA*VER_zm')'+SHIFT;
