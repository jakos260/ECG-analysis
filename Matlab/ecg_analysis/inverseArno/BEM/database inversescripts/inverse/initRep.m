function rep = initRep(GEOM,dep)



rho=zeros(size(dep));

for i=1:length(dep)
    a=find(GEOM.DIST(i,:) < 30 & GEOM.DIST(i,:) > 0);
%     if GEOM.typ(i) == 2
        rho(i) = sum( (dep(a) - dep(i) ) ./ (GEOM.DIST(a,i).^2));
   
end
ari = GEOM.SPECS.time_apexT - mean(dep) + rho * 5 + GEOM.SPECS.repCorrection;

rep =ari + dep;
if 1
    figure(50)
    clf
    plot(dep,rep-dep,'.')
    linear_regression(dep,ari,gcf) 
end
PSIREF = GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave);
t=0:max(rep)+100;
T=ones(length(GEOM.VER),1)*t;
if isempty(strfind(GEOM.type,'atria'))    
    [~,ttmref]=gettdom(PSIREF,floor(max(dep)));	
    S=getSmode(T,dep,rep, GEOM.SPECS ,4);
    PHI=lowpassma(GEOM.AMA*S,5);	
    [~,ttm]=gettdom(PHI,floor(max(dep)));

    rep=rep+(ttmref-ttm);

    figure(23)
    subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');
    ylabel('ARI [ms]'); axis tight
    subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');
    ylabel('repolarization time [ms]');axis tight
else
%     figure(23);   
%     subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
%     subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight

end


% figure(21);showpatch(GEOM.VER,GEOM.ITRI,rep);
% figure(22);showpatch(GEOM.VER,GEOM.ITRI,ari);

