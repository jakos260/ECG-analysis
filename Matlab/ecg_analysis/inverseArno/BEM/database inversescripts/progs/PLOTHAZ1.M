echo off
clear
clf
%TAB=loadasci('hm9498.txt');
h=loadasci('hv9498.dat');
t=loadasci('h9498.tvals');
logh=log(h);
plot(t,logh)
scal=axis;
scal(4)=0.
scal(2)=120;
axis(scal);
hold on
size(t)
size(logh)
pause
alpha=3.13;
beta =-14.85;
gamma=-7.383;
tmax=120;
tt=t(1):tmax;
tts=(tt-tmax)/tmax;
T(1:58,1)=(t-tmax)/tmax;
T(1:58,2)=T(1:58,1).^2;
T(1:58,3)=T(1:58,1).^3;
sol=T\logh;
alpha=sol(1)
beta=sol(2)
gamma=sol(3)

fun=alpha*tts+beta*tts.^2+gamma*tts.^3;
plot(tt,fun,'red')
%res=norm(fun-logh)
%relres=norm(fun-logh)/norm(logh)
