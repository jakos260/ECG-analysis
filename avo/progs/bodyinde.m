% bodyinde.m
clear
clf
lbeg=1.5;
lend=2.;
jsteps=41;
lincr=(lend-lbeg)/(jsteps-1);
mbeg=50;
mend=100;
isteps=41;
mincr=(mend-mbeg)/(isteps-1);
[L,M]=meshgrid([lbeg:lincr:lend],[mbeg:mincr:mend]);
INDEX=zeros(jsteps,isteps);
for i=1:isteps,
m=mbeg+(i-1)*mincr;
for j=1:jsteps,
l=lbeg+(j-1)*lincr;
INDEX(i,j)=(m/10-l*10+9.65)/1.1736;
end
end
contour(L,M,INDEX,[-3: .5 :3])
hold on
plot(1.66,69.67,'red+')
plot(1.775,81.14,'blue+')
text(1.75,102,'(m/10-10*l+9.65)/1,1736:  Linear Index')
text(1.75,45,'steps  -3 : .5 : 3 (unit:SD-population)')
TAB=loadasci('constitu.sor');
disp=(TAB(:,4)/10-TAB(:,3)*10+9.65)/1.1736;
mean(disp)
std(disp)

