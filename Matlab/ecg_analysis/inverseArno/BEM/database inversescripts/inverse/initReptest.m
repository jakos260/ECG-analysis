function rep = initReptest(GEOM,dep,foci)

 

belongsTo= zeros(size(dep));
for i=1:length(dep)
    minD=100000;
    for j=1:length(foci)    
        if GEOM.DIST2W(foci(j),i) < minD
            minD = GEOM.DIST2W(foci(j),i);
            ifocus = j;
        end
    end        
    belongsTo(i) = ifocus;
    
end

for i =1:length(foci)
    focus = foci(i);
    mu = 500/diff(range(dep(belongsTo==i)));
    dist = GEOM.DIST(focus,belongsTo==i);
    alpha=(exp(-dist./mu)-exp((dist-max(dist))./(mu)));
    plot(dep(belongsTo==i),alpha,'.')
    if GEOM.typ(focus) == 2
        ari(belongsTo==i) = GEOM.SPECS.time_apexT + 20 - mean(dep(belongsTo==i)) + alpha * 0.7 * diff(range(dep(belongsTo==i)));
    else
        ari(belongsTo==i) = GEOM.SPECS.time_apexT - 20- mean(dep(belongsTo==i)) + alpha * 0.4 * diff(range(dep(belongsTo==i)));
    end
end
ari = ari';
% rho=zeros(size(dep));
% 
% for i=1:length(dep)
%     a=find(GEOM.DIST(i,:) < 10 & GEOM.DIST(i,:) > 0);
%     rho(i) = sum((dep(a) - dep(i)) ./ (GEOM.DIST2W(a,i).^2));
% end
% ari = GEOM.SPECS.time_apexT - mean(dep) + rho*7 + GEOM.SPECS.repCorrection;

rep =ari + dep;

if 1
    figure(50)
    clf
    plot(dep,rep-dep,'.')
    linear_regression(dep,rep- dep,gcf) 
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

