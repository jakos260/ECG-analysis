% subscripy fourplots
% used by demogr

clf
subplot(2,2,1)
plot(t,ptot(1:101,2),'white')
scal=axis;
hold on
plot(t,pmi(1:101,2),'blue')
plot(t,pfi(1:101,2),'red')
tekst=sprintf('%i',year);
text(105,1.1*scal(4),tekst)

subplot(2,2,2)
plot(t,born(1:101),'white')
scal=axis;
hold on
plot([0,0],[0,scal(4)],'white')
plot(t,2.5*scal(4)*vrt(1:101,2),'g')
plot([98,100],[0.5*scal(4),0.5*scal(4)],'g')
tekst=sprintf('%s','0.2');
text(105,0.5* scal(4),tekst);
tekst=sprintf('%s','%');
text(105,0.9*scal(4),tekst);

subplot(2,2,3)
plot(t,mdeath(1:101),'blue')
scal=axis;
hold on
plot([0,0],[0,scal(4)],'blue')
plot(t,scal(4)*stm(1:101,2),'g')
tekst=sprintf('%s','0.5');
text(105,0.5*scal(4),tekst)

subplot(2,2,4)
plot(t,fdeath(1:101),'red')
axis(scal);
hold on
plot([0,0],[0,scal(4)],'red')
plot(t,scal(4)*stf(1:101,2),'g')
tekst=sprintf('%s','0.5');
text(105,0.5*scal(4),tekst)
