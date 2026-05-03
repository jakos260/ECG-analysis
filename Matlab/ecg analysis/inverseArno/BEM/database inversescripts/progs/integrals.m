% integrals.m
% version: 20030107
% a. van Oosterom
% test various integrqals

clear
clf
% 1/R related integrals
% int(log(a+sqrt(b^2+x^2)); |a| < |b|
a=1; c=.1;
b=sqrt(a^2+c^2);

x=0:.01:1;
r=sqrt(b^2+x.^2);
y=log(a+r);
plot(x,y)
pause
z1=introws(y)/100;
plot(x,z1)
hold on
z2=x.*log(a+r)+a*log(x+r)+2*c*atan(c*x./((a+b)*(b+r)))-x;
z2=z2-z2(1)+z1(1);
plot(x,z2,'r')
z3=cumsum(y)/100;
z3=z3-z3(1);
plot(x,z3,'g')
pause
clear
clf
% int(log(a+sqrt(b^2+x^2)); |a| > |b|

b=1; c=3;
a=sqrt(b^2+c^2);

x=0:.01:1;
r=sqrt(b^2+x.^2);
y=log(a+r);
plot(x,y)
pause
z1=introws(y)/100;
plot(x,z1)
hold on
arg=abs( ((a-b)*x+(b+r)*c)./((a-b)*x-(b+r)*c) );
z2=x.*log(a+r)+a*log(x+r)-c*log(arg)-x;
z2=z2-z2(1)+z1(1);
plot(x,z2,'r')


