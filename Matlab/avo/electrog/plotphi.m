% plotphi.m
% 20030205
% script of electrogr
clf
% draw wall and signals at nn observation points
xb=b*sin(2*phi);
yb=b*cos(2*phi);
xa=a*sin(2*phi);
ya=a*cos(2*phi);
plot(xb-.5,yb)
axis([-1.5 1.5 -1.5 1.5])
axis square
hold on
plot(xa-.5,ya,'r')
plot(xobs-.5,zobs,'*k')

% plot potentials at nn nodes 
for i=1:nn,
   curve=PHIOBS(i,:);
   shift=1.2-(i-1)*.3;
   obs=[xobs(i) yobs(i) zobs(i)];
   plot(tijd,shift+scal*curve,'r');
   plot(tijd,shift+0*curve,'k:')
   plot([xobs(i)-.5 0.55],[zobs(i) shift],':')
end

cal=-1.4+.5*scal;
calm=-1.4+.46*scal;
calib=plot([-1.25 -1.25 -1.23 -1.25 -1.27],[-1.4 cal calm cal calm]);
if inho==0, text(-1.2, -1.25 ,'20 mV (inf medium)'),
else,text(-1.2, -1.25 ,['20 mV; sigma: ' num2str(sigm)]),
end

if stim=='epi ', title('epicardial stimulus'),
   plot(-.5 , b,'*r')
else,
   title('endocardial stimulus'),
   plot(-.5 , a,'*r')
end

