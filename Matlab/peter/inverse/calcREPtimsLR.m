function rep=calcREPtimsLR(GEOM,depfinal,shape_mu)
% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

locMinArea=30;
locmin=[];
for i=1:length(depfinal)
	if depfinal(i)<max(depfinal)	
		if all(depfinal(GEOM.DIST(:,i)<locMinArea)-depfinal(i)>=0)		
			locmin=[locmin i];
		end
	end
end
sourcelocmin=locmin;

locmin=[];
for i=1:length(depfinal)
	if all(depfinal(GEOM.DIST(:,i)<locMinArea)-depfinal(i)<=0) 
		locmin=[locmin i];
	end
end
sinklocmin=locmin;
%%
x=10*ones(length(depfinal),length(sinklocmin)*length(sourcelocmin));
mu=-x;
depK=x;
D=x;
Dk=x;
Dc=x;
dep=depfinal;
for i=1:length(dep)
	k=1;
	for ic=1:length(sourcelocmin)
		for ik=1:length(sinklocmin)
			sinkVeri=sinklocmin(ik);
			sourceVeri=sourcelocmin(ic);
            DIST = GEOM.DIST;
            
			v=(DIST(sinkVeri,sourceVeri))./((dep(sinkVeri)-dep(sourceVeri)));
			L=DIST(sinkVeri,sourceVeri);
			depK(i,k)=(dep(sinkVeri)-dep(sourceVeri));
			depL=L/v;
			mu(i,k)=depL;
			x(i,k)=(dep(i)-dep(sourceVeri))/depL;
			D(i,k)=DIST(sinkVeri,i)+GEOM.DIST(sourceVeri,i);
			Dk(i,k)=DIST(sinkVeri,i);
			Dc(i,k)=DIST(sourceVeri,i);			
			k=k+1;
		end
	end
end
tmpx=x;
tmpmu=mu;
tmpdepK=depK;
%%
x=10*ones(size(dep));
mu=-x;
depK=ones(size(depfinal));
for i=1:length(dep)
	% find all source and sinks within a distance of 40 mm
% 	sel=find(Dk(i,:)<=min(Dk(i,:))+40 & Dc(i,:)<=min(Dc(i,:))+40);	
    sel=find(Dk(i,:)<=min(Dk(i,:))*2 & Dc(i,:)<=min(Dc(i,:))*2);	
	depK(i)=mean(tmpdepK(i,sel));
    depK(i)=mean(tmpdepK(i,sel)) - diff(range(tmpdepK(i,sel)))/4;
	mu(i)=mean(tmpmu(i,sel));
	x(i)=mean(tmpx(i,sel));
end
mu=mu.*shape_mu;
alpha=(exp(-x./mu)-exp((x-max(x))./(mu)));
%limit extremes
maxmin=3*std(alpha)+mean(alpha);	
alpha(alpha>maxmin)=maxmin;	alpha(alpha<-maxmin)=-maxmin;

alpha=alpha.*depK;	
alpha(GEOM.Rfreewallver==0) = alpha(GEOM.Rfreewallver==0) + 5;
alpha(GEOM.Rfreewallver==1) = alpha(GEOM.Rfreewallver==1) - 5;

alpha = alpha / max(abs(alpha)) * max(depfinal) /2;
alpha = ones(size(depfinal)) *  mean(depfinal) - depfinal ;
alpha(GEOM.Rfreewallver==0) = 0.9*alpha(GEOM.Rfreewallver==0) + 10;
alpha(GEOM.Rfreewallver==1) = 0.3*alpha(GEOM.Rfreewallver==1) - 10;

%%
if isempty(strfind(GEOM.type,'atria'))
    PSIREF=baselinecor(GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5)));
    t=0:size(PSIREF,2)-1;T=ones(size(GEOM.DIST,1),1)*t;
    [tdom,ttmref]=gettdom(PSIREF,floor(max(depfinal)));	
    rep=ttmref+(depfinal-mean(depfinal))+alpha;
    S=getSmode(T,depfinal,rep,GEOM.pS,[],4);
    PHI=lowpassma(GEOM.AMA*S,5);	
    % plot(rms(PSIREF)); hold on; plot(rms(PHI),'r');
    [tdom,ttm]=gettdom(PHI,floor(max(depfinal)));
    if diff(abs(ttmref-ttm))<50
        rep=rep+(ttmref-ttm);
    end
    figure(1)
    subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
    subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight
else
    figure(1);
    
    rep=GEOM.specs(5)-GEOM.specs(2)+(depfinal-mean(depfinal))+ 0.2*alpha;
    subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
    subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight

end
f=polyfit(depfinal-min(depfinal),rep-depfinal,1);
fl=polyfit(depfinal(GEOM.Rfreewallver==0)-min(depfinal(GEOM.Rfreewallver==0)),rep(GEOM.Rfreewallver==0)-depfinal(GEOM.Rfreewallver==0),1);
fr=polyfit(depfinal(GEOM.Rfreewallver==1)-min(depfinal(GEOM.Rfreewallver==1)),rep(GEOM.Rfreewallver==1)-depfinal(GEOM.Rfreewallver==1),1);
disp(['slope ARI=f(dep): ' num2str([f(1) fl(1) fr(1)])])

