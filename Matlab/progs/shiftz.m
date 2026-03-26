% shiftz.m
% script of fitgeometries
% 20090119

newshift=get(slz,'Val');
shift(3)=shift(3)+newshift;

VERA(:,3)=VERA(:,3)+newshift;

nver=size(VERA,1);
VER(1:nver,:)=VERA;

set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);
set(uit6,'string',num2str(shift(3)));
set(slz,'val',0)

figure(2)
crossub
