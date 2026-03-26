% dirstat.m
% statistics of directional data
clear
DATA=loadmat('angles.lst');
dim=size(DATA);
npp=dim(1);
nvar=dim(2);
males=[1 2 3 4 7 8 9 10 11 17 19 20 21 24 25];
females=[5 6 12 13 14 15 16 18 22 23];
MALES=DATA(males,:);
FEMALES=DATA(females,:);
nmales=15;
nfemales=10;
col1=12;
col2=13;
azi=DATA(:,col1)*pi/180;
ele=DATA(:,col2)*pi/180;
x=-cos(ele).*sin(azi);
y=cos(ele).*cos(azi);
z=sin(ele);
xm=mean(x);
ym=mean(y);
zm=mean(z);
precision=norm([xm ym zm])
azim=atan2(-xm,ym)*180/pi
elem=asin(zm)*180/pi

azi=FEMALES(:,col1)*pi/180;
ele=FEMALES(:,col2)*pi/180;
x=-cos(ele).*sin(azi);
y=cos(ele).*cos(azi);
z=sin(ele);
xm=mean(x);
ym=mean(y);
zm=mean(z);
precision=norm([xm ym zm])
azim=atan2(-xm,ym)*180/pi
elem=asin(zm)*180/pi

azi=MALES(:,col1)*pi/180;
ele=MALES(:,col2)*pi/180;
x=-cos(ele).*sin(azi);
y=cos(ele).*cos(azi);
z=sin(ele);
xm=mean(x);
ym=mean(y);
zm=mean(z);
precision=norm([xm ym zm])
azim=atan2(-xm,ym)*180/pi
elem=asin(zm)*180/pi
