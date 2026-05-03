echo off
clear
clf
TAB=loadasci('sterfte.tab');
tvals=TAB(1:77,1);

nsteps=151;
amin=0.4;
amax=0.46;

cmin=24.;
cmax=28.;
adel=(amax-amin)/(nsteps-1);
cdel=(cmax-cmin)/(nsteps-1);
steps=0.: .005: .1;
[A,C]=meshgrid(amin:adel:amax,cmin:cdel:cmax);
for j=1:nsteps,
a=amin+(j-1)*adel;
for i=1:nsteps,
c=cmin+(i-1)*cdel;
est= (1./(1.+a/100*(106.5-tvals))).^c;
FUN(i,j)=norm(est-TAB(1:77,2));
end
end
contour(A,C,FUN,steps)
%mesh(A,C,FUN)
xlabel('alpha')
ylabel('p')
%colorbar

