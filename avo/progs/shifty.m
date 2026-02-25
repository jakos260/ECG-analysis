% shifty.m
% script of fitgeometries
% 20090119

newshift=get(sly,'Val');
shift(2)=shift(2)+newshift;
VERA(:,2)=VERA(:,2)+newshift;
nver=size(VERA,1);
VER(1:nver,:)=VERA;

set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);
set(uit5,'string',num2str(shift(2)));
set(sly,'val',0);
figure(2)
crossub



