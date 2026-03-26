% plothaz.m
% plotje van sterfte tabellen
echo off
clf
tmax=120;
TABM=loadasci('hned81m.dat');
t=loadasci('tned81.dat');
semilogy(t,TABM(:,1),'b')
scal=axis;
scal(1)=30;
scal(2)=tmax;
axis(scal);
hold on
TABV=loadasci('hned81v.dat');
semilogy(t,TABV(:,1),'r')
tmax=110
est=30:tmax;
t=30:tmax;
est=exp(0.097*(t-tmax));
semilogy(t,est)
tmax=120;
t=30:tmax;
est=t;
est=exp(0.037*(t-tmax)-.00053*(t-tmax).^2);
semilogy(t,est,'g')


