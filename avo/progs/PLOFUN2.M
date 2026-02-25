% file plofun2.m
% plot van norm van residual tussen
% exp hazard data en negatief eponentieel model
echo off
clear
clf
TAB=loadasci('sterfte.tab');
tvals=TAB(1:77,1);

nsteps=41;
amin=0.08;
amax=0.12;
tmin=106.;
tmax=109.;
adel=(amax-amin)/(nsteps-1);
tdel=(tmax-tmin)/(nsteps-1);
incr=.05;
steps=0.: incr : 1.;
[A,T]=meshgrid(amin:adel:amax,tmin:tdel:tmax);
small=1.e+6;
asmall=amin;
tsmall=tmin;
for j=1:nsteps,
a=amin+(j-1)*adel;
for i=1:nsteps,
t=tmin+(i-1)*tdel;
est= exp(-a.*(t-tvals));
res=est-TAB(1:77,2);
nr=norm(res);
FUN(i,j)=nr;
if nr < small;
small=nr;
asmall=a;
tsmall=t;
end
end
end
contour(A,T,FUN,steps)
[small,asmall,tsmall]
%mesh(A,T,FUN)
xlabel('alpha');
ylabel('tmax');
ytekst=tmax+0.04*(tmax-tmin);
xtekst=amin+0.05*(amax-amin);
tekst=sprintf('min: %0.5g',small);
text(xtekst,ytekst,tekst);
tekst=sprintf('at  alpha: %0.5g',asmall);
text(xtekst+0.35*(amax-amin),ytekst,tekst);
tekst=sprintf('    and  tmax: %0.5g',tsmall);
text(xtekst+0.6*(amax-amin),ytekst,tekst);
%colorbar
ytekst=tmin-0.1*(tmax-tmin);
tekst=sprintf('step: %0.5g',incr);
text(amin,ytekst,tekst);


