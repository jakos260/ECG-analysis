function rep=calcREPtims(GEOM,depfinal,shape_mu)
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
if ~isempty(strfind(GEOM.type,'atria'))
	sourcelocmin=sourcelocmin(1);
	sinklocmin=find(GEOM.DIST(:,sourcelocmin)==max(GEOM.DIST(:,sourcelocmin)));
end
%%
x=10*ones(length(depfinal),length(sinklocmin)*length(sourcelocmin));
mu=-x;
depK=x;
D=x;Dk=x;Dc=x;
dep=depfinal;
for i=1:length(dep)
	k=1;
	for ic=1:length(sourcelocmin)
		for ik=1:length(sinklocmin)
			sinkVeri=sinklocmin(ik);
			sourceVeri=sourcelocmin(ic);
			v=(GEOM.DIST(sinkVeri,sourceVeri))./((dep(sinkVeri)-dep(sourceVeri)));
			L=GEOM.DIST(sinkVeri,sourceVeri);
			depK(i,k)=(dep(sinkVeri)-dep(sourceVeri));
			depL=L/v;
			mu(i,k)=depL;
			x(i,k)=(dep(i)-dep(sourceVeri))/depL;
			D(i,k)=GEOM.DIST(sinkVeri,i)+GEOM.DIST(sourceVeri,i);
			Dk(i,k)=GEOM.DIST(sinkVeri,i);
			Dc(i,k)=GEOM.DIST(sourceVeri,i);			
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
	sel=find(Dk(i,:)<=min(Dk(i,:))+40 & Dc(i,:)<=min(Dc(i,:))+40);	
	depK(i)=mean(tmpdepK(i,sel));
	mu(i)=mean(tmpmu(i,sel));
	x(i)=mean(tmpx(i,sel));
end
mu=mu*shape_mu;
alpha=(exp(-x./mu)-exp((x-max(x))./(mu)));
%limit extremes
maxmin=1.5;	
alpha(alpha>maxmin)=maxmin;	alpha(alpha<-maxmin)=-maxmin;
alpha=alpha.*depK;	

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
    f=polyfit(depfinal-min(depfinal),rep-depfinal,1);
    disp(['slope ARI=f(dep): ' num2str(f(1))])

else
%     figure(1);
%     
    rep= GEOM.SPECS.onsetqrs - GEOM.SPECS.onsetP + (depfinal-mean(depfinal))+ 0.2*alpha;
%     subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
%     subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight

end

