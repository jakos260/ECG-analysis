% rot_x.m
% callback of fitgeometries
% rotate VER1 about its center of gravity by rotx degrees around positive x-axis; positive: anti-clockwise
% total of successive rotations are stored in ROTS
% 20090119

rotx=get(slrotx,'Val');

rotx=rotx/180*pi;

nver=size(VERA,1);
vermean=mean(VERA);
SHIFT=ones(nver,1)*vermean;
VERA_zm=VERA-SHIFT;
ROMA=[1 0 0; 0 cos(rotx) sin(rotx); 0 -sin(rotx) cos(rotx)];
ROTS=ROMA*ROTS;
VERA=(ROMA*VERA_zm')'+SHIFT;

VER(1:nver,:)=VERA;
set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);

set(slrotx,'val',0)
figure(2)
crossub
figure(1)