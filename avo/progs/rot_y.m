% rot_y.m
% callback of fitgeometries
% rotate VER1 by roty degrees around positive y-axis; positive: clockwise
% total of successive rotations stored in ROTS
% 20090119

roty=get(slroty,'Val');
roty=roty/180*pi;

nver=size(VERA,1);
vermean=mean(VERA);
SHIFT=ones(nver,1)*vermean;
VERA_zm=VERA-SHIFT;

ROMA=[ cos(roty) 0 -sin(roty); 0 1 0; sin(roty) 0 cos(roty)];
ROTS=ROMA*ROTS;
VERA=(ROMA*VERA_zm')'+SHIFT;

VER(1:nver,:)=VERA;
set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);

set(slroty,'val',0)
figure(2)
crossub
figure(1)