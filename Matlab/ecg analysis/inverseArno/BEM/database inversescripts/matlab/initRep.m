function rep = initRep(GEOM,dep,showplot)

if ~exist('showplot','var'), showplot = 1; end

rho = zeros(size(dep));
for i = 1:length(dep),
    a       = find(GEOM.DIST(i,:) < 50 & GEOM.DIST(i,:) > 0);
    rho(i)  = sum((dep(a) - dep(i)) ./ (GEOM.DIST(a,i).^2));
end

ari = GEOM.SPECS.time_apexT - mean(dep) + rho*7 + GEOM.SPECS.repCorrection;
rep = ari + dep;

PSIREF  = GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave);
t       = 0:max(rep)+100;
T       = ones(length(GEOM.VER),1)*t;

if isempty(strfind(GEOM.type,'atria')),
    [xdum,ttmref]   = gettdom(PSIREF,floor(max(dep)));
    S               = getSmode(T,dep,rep, GEOM.SPECS ,4);
    PHI             = lowpassma(GEOM.AMA*S,5);
    [xdum,ttm]      = gettdom(PHI,floor(max(dep)));
    
    rep = rep + (ttmref-ttm);
    
    if showplot,
        figure(23)
        subplot(2,1,1); plot(dep,rep-dep,'.'); xlabel('depolarization time [ms]');
        ylabel('ARI: rep-dep [ms]'); axis tight
        subplot(2,1,2); plot(dep,rep,'.'); xlabel('depolarization time [ms]');
        ylabel('repolarization time [ms]'); axis tight
    end    
end