% setheta.m
% callback of fitgeometries

theta=get(sltheta,'Val');
VER1=rotash(VERINP,[phi/180 theta/180 gamma/180],shift);
nver1=size(VER1,1);
VER(1:nver1,:)=VER1;
set(hs,'Vertices',VER)
set(ht,'Pos',[1.03*VER(node,1) 1.03*VER(node,2) 1.03*VER(node,3)]);
set(uit2,'string',num2str(theta));
figure(2);
crossnew
figure(1)