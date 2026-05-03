% plotphi.m
% 20030205
% script of electrogr
figure(1)
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

% plot observation points; include NEEDLE
% plot(xobs-.5,zobs,'*k')

% plot observation points; EXclude NEEDLE
plot(xobs(1:nn+1) -.5,zobs(1:nn+1),'*k')

% plot potentials at nn nodes 
for i=1:nn,
   curve=PHIOBS(i,:);
   shift=1.2-(i-1)*.3;
   plot(tijd,shift+scal*curve,'r');
   plot(tijd,shift+0*curve,'k:')
   plot([xobs(i)-.5 0.55],[zobs(i) shift],':')
end
   % potential centre sphere
   curve=PHIOBS(nn+1,:);
   shift=-1.5;
   plot(tijd+shift,scal*curve-0.2,'k');
   plot(tijd+shift,0*curve-0.2,'k:')
   
% plot potential profiles along needle
figure(2)
clf
for proft=1:9,
    tnow=(proft-1)*2+6;
    plot(OBS(nn+1:nobs,3), proft*20+20*PHIOBS(nn+1:nobs,tnow),'r')
    hold on
    text(1.3,proft*20-5,['t= ' num2str(tnow)])
    plot([OBS(nn+1,3) OBS(nobs,3)], proft*20+[0 0],'k:')
    if inho==1,
       plot(OBS(nn+1:nobs,3), proft*20+20*PHIOBS_sec(nn+1:nobs,tnow),'g');
       plot(OBS(nn+1:nobs,3), proft*20+20*PHIOBS_inf(nn+1:nobs,tnow),'b');
    end
 end
 assen=axis;
 assen(1:4)=[0 OBS(nobs,3) 0 190];
 axis(assen)
 title(' potential profiles along needle')
 
 %figure(3)
 %clf
 %plot(40*PHIOBS([1 110  nn+1],:)')
 %title('potentials at nodes 1 110 nn+1')

figure(1)

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


