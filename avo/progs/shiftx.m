% shiftx.m
% script of fitgeometries

newshift=get(slx,'Val');

shift(1)=shift(1)+newshift;


VERA(:,1)=VERA(:,1)+newshift;

nver=size(VERA,1);

VER(1:nver,:)=VERA;

set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);

set(slx,'val',0)

set(uit4,'string',num2str(shift(1)));

figure(2)
crossub
