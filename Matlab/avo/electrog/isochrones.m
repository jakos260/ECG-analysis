% isochrones.m
% plot isochrones within annulus
% 20051105
clear
a=40; b=50;

% epicardial focus
pnt1=[0 b 0];
pnt2=[0 a 0];

nr=15;
nphi=61;

[VER,ITRI]=make_annulus(a,b,nr,nphi,0);

nver=size(VER,1);

for i=1:nver;
    VALS(i,1)=distance_annulus(pnt1,VER(i,:),a,b);
    VALS(i,2)=distance_annulus(pnt2,VER(i,:),a,b);
end

VALS(:,3)=min(VALS')';

% ms if v=1m/s
VALS=VALS/1;
 
grsw=0;

figure(1)
clf
cmap='tims.mcm';
zebra=-5;
lsw=0;
iview=9;

showPatch(VER,ITRI,VALS(:,2));view(0,90)
