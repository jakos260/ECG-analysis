function rep = initRep(GEOM,dep,useleads)

lpass= 5;
rho=zeros(size(dep));
for i=1:length(dep)
    a=find(GEOM.DIST2W(i,:) < 50 & GEOM.DIST2W(i,:) > 0);
    rho(i) = sum((dep(a) - dep(i)) ./ (GEOM.DIST2W(a,i).^2));
end
qrsduration = GEOM.specs(3)-GEOM.specs(2)+1;
% ari = qrsduration + GEOM.specs(4) - mean(dep) + rho.*(80 /diff(range(rho))) + GEOM.specs(8);
ari = qrsduration + GEOM.specs(4) - mean(dep) + rho*7 + GEOM.specs(8);

rep =ari + dep;


% PSIREF = GEOM.BSM(useleads,GEOM.specs(2):GEOM.specs(5));
PSIREF = GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5)); % oostep1: are only usedleads
t=0:size(PSIREF,2)-1;
T=ones(length(GEOM.VER),1)*t;
if isempty(strfind(GEOM.type,'atria'))    
    [tdom,ttmref]=gettdom(PSIREF,floor(max(dep)));	
    S=getSmode(T,dep,rep,GEOM.pS,[],4);
    PHI=lowpassma(GEOM.AMA*S,5);	
    [tdom,ttm]=gettdom(PHI,floor(max(dep)));
%     if diff(abs(ttmref-ttm))<50
        rep=rep+(ttmref-ttm);
%     end
    figure(23)
    subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
    subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight
else
    figure(23);
   
    subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
    subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight

end


figure(21);showPatch(GEOM.VER,GEOM.ITRI,rep);
figure(22);showPatch(GEOM.VER,GEOM.ITRI,ari);

% PSIA =lowpassma((GEOM.AMA(useleads,:))*getSmode(T,dep,rep,GEOM.pS,[],4),lpass);
PSIA =lowpassma(GEOM.AMA*getSmode(T,dep,rep,GEOM.pS,[],4),lpass); % oostep1: AMA already only used BSM leads
figure(23);clf
normphi=norm(PSIREF,'fro');
sigplot(PSIREF,'',GEOM.LAY,1,'b',1,0);
hold on
sigplot(PSIA,'',GEOM.LAY,1,'r',1,0);
figure(24);clf;plot(rms(PSIREF),'b');hold on;plot(rms(PSIA),'r')
COR=corrcoef(PSIA,PSIREF);
disp(['mean/std rep ' num2str([mean(rep) std(rep)],3) '    cor/rd ' num2str([COR(2,1) norm(PSIREF - PSIA,'fro')/normphi],3)])

