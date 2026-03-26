% rot_z.m
% callback of fitgeometries
% rotate VER1 by rotz degrees around positive z-axis; positive: anti-clockwise
% total of successive rotations stored in ROTS
% 20060219

rotz=get(slrotz,'Val');
rotz=rotz/180*pi;

nver=size(VERA,1);
vermean=mean(VERA);
SHIFT=ones(nver,1)*vermean;
VERA_zm=VERA-SHIFT;
ROMA=[cos(rotz)  sin(rotz) 0;  -sin(rotz) cos(rotz) 0; 0 0 1];
ROTS=ROMA*ROTS;
VERA=(ROMA*VERA_zm')'+SHIFT;

VER(1:nver,:)=VERA;
set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);

set(slrotz,'val',0)
figure(2)
crossub
figure(1)