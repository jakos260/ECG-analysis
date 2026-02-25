function rep = initRep(GEOM,dep)



rho=zeros(size(dep));

for i=1:length(dep)
    a=find(GEOM.DIST(i,:) < 30 & GEOM.DIST(i,:) > 0);
    rho(i) = sum( (dep(a) - dep(i) ) ./ (GEOM.DIST(a,i).^1.8));
   
end
PSIREF = GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave);
t=0:size(PSIREF,2)+100;
T=ones(length(GEOM.VER),1)*t;
[~,ttmref]=gettdom(PSIREF,floor(max(dep)));	
rrms = rms(PSIREF);
endt = find(rrms(ttmref:end) < rrms(ttmref)*0.2);endt=endt(1);
rangeTwave = endt;
fact = rangeTwave/diff(range(rho));

ari = GEOM.SPECS.time_apexT - mean(dep) + rho * fact - GEOM.SPECS.repCorrection;

ari = GEOM.SPECS.time_apexT - mean(dep) + dep*-0.5 - GEOM.SPECS.repCorrection;
rep =ari + dep;
if 1
    figure(50)
    clf
    plot(dep,rep-dep,'.')
    linear_regression(dep,ari,gcf);
end
if isempty(strfind(GEOM.type,'atria'))    
        S=getSmode(T,dep,rep, GEOM.SPECS ,4,GEOM.ADJsurf);
    [~,ttmref]=gettdom(PSIREF,floor(max(dep)));	
    PHI=lowpassma(GEOM.AMA*S,5);	
    [~,ttm]=gettdom(PHI,floor(max(dep)));

    rep=rep+(ttmref-ttm);
    S=getSmode(T,dep,rep, GEOM.SPECS ,4,GEOM.ADJsurf);
    [~,ttmref]=gettdom(PSIREF,floor(max(dep)));	
    PHI=lowpassma(GEOM.AMA*S,5);	
    [~,ttm]=gettdom(PHI,floor(max(dep)));
    rep=rep+(ttmref-ttm);
    S=getSmode(T,dep,rep, GEOM.SPECS ,4,GEOM.ADJsurf);
    [~,ttmref]=gettdom(PSIREF,floor(max(dep)));	
    PHI=lowpassma(GEOM.AMA*S,5);	
    [~,ttm]=gettdom(PHI,floor(max(dep)));
    
    figure(23);clf
    subplot(3,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');
    ylabel('ARI [ms]'); axis tight
    subplot(3,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');
    ylabel('repolarization time [ms]');axis tight
    subplot(3,1,3);plot(rms(PHI));hold on; plot(rms(PSIREF))
else
%     figure(23);   
%     subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
%     subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight

end


% figure(21);showpatch(GEOM.VER,GEOM.ITRI,rep);
% figure(22);showpatch(GEOM.VER,GEOM.ITRI,ari);

