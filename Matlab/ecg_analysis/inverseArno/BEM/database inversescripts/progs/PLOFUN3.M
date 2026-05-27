% file plofun3.m
% plot van Rosenbrock's functie
echo off
clear
clf

nsteps=41;
xmin=-1.5;
xmax=1.5;
ymin=-0.5;
ymax=1.5;
xdel=(xmax-xmin)/(nsteps-1);
ydel=(ymax-ymin)/(nsteps-1);
incr=4.;
steps=0.: incr : 100.;
[X,Y]=meshgrid(xmin:xdel:xmax,ymin:ydel:ymax);
small=1.e+6;
xsmall=xmin;
ysmall=ymin;
for j=1:nsteps,
x=xmin+(j-1)*xdel;
for i=1:nsteps,
y=ymin+(i-1)*ydel;
FUN(i,j)=100*(y-x^2)^2+(1-x)^2;
if FUN(i,j) < small;
small=FUN(i,j);
xsmall=x;
ysmall=y;
end
end
end
contour(X,Y,FUN,steps)
[small,xsmall,ysmall]
%mesh(X,Y,FUN)
xlabel('x-vals');
ylabel('y-vals');
ytekst=ymax+0.04*(ymax-ymin);
xtekst=xmin+0.05*(xmax-xmin);
tekst=sprintf('min: %0.5g',small);
text(xtekst,ytekst,tekst);
tekst=sprintf('at  x: %0.5g',xsmall);
text(xtekst+0.35*(xmax-xmin),ytekst,tekst);
tekst=sprintf('    and  y: %0.5g',ysmall);
text(xtekst+0.6*(xmax-xmin),ytekst,tekst);
%colorbar
ytekst=ymin-0.1*(ymax-ymin);
tekst=sprintf('step: %0.5g',incr);
text(xmin,ytekst,tekst);
pause
% now compute and plot norm of the gradient
incr=.1;
steps=0.: incr : 5.;
small=1.e+6;
xsmall=xmin;
ysmall=ymin;
for j=1:nsteps,
x=xmin+(j-1)*xdel;
for i=1:nsteps,
y=ymin+(i-1)*ydel;
valx=-400*(y-x^2)-(1-x);
valy=200*(y-x^2);
FUN(i,j)=sqrt(valx^2+valy^2);
if FUN(i,j) < small;
small=FUN(i,j);
xsmall=x;
ysmall=y;
end
end
end
contour(X,Y,FUN,steps)
[small,xsmall,ysmall]
%mesh(X,Y,FUN)
xlabel('x-vals');
ylabel('y-vals');
ytekst=ymax+0.04*(ymax-ymin);
xtekst=xmin+0.05*(xmax-xmin);
tekst=sprintf('min: %0.5g',small);
text(xtekst,ytekst,tekst);
tekst=sprintf('at  x: %0.5g',xsmall);
text(xtekst+0.35*(xmax-xmin),ytekst,tekst);
tekst=sprintf('    and  y: %0.5g',ysmall);
text(xtekst+0.6*(xmax-xmin),ytekst,tekst);
%colorbar
ytekst=ymin-0.1*(ymax-ymin);
tekst=sprintf('step: %0.5g',incr);
text(xmin,ytekst,tekst);


