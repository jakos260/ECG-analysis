function rep = initRep(GEOM,dep)

lpass= 5;
rho=zeros(size(dep));
for i=1:length(dep)
    a=find(GEOM.DIST2W(i,:) < 50 & GEOM.DIST2W(i,:) > 0);
    rho(i) = sum((dep(a) - dep(i)) ./ (GEOM.DIST2W(a,i).^2)); % waarom ^2, en niet delen door aantal nodes?
end
qrsduration = GEOM.specs(3)-GEOM.specs(2)+1;
% ari = qrsduration + GEOM.specs(4) - mean(dep) + rho.*(80 /diff(range(rho))) + GEOM.specs(8);
ari = qrsduration + GEOM.specs(4) - mean(dep) + rho*7 + GEOM.specs(8);

rep =ari + dep;


PSIREF = GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
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
    fh=figure(23);set(fh,'Name','Init: ARI/Rep vs Dep');
    subplot(2,1,1);plot(dep,rep-dep,'.');drawnow;xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight; 
    subplot(2,1,2);plot(dep,rep,'.');drawnow;xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight
else
       fh=figure(23);set(fh,'Name','Init: ARI/Rep vs Dep');
   
    subplot(2,1,1);plot(dep,rep-dep,'.');xlabel('depolarization time [ms]');ylabel('ARI [ms]'); axis tight
    subplot(2,1,2);plot(dep,rep,'.');xlabel('depolarization time [ms]');ylabel('repolarization time [ms]');axis tight

end


fh=figure(21);showPatch(GEOM.VER,GEOM.ITRI,rep);set(fh,'Name','Init:Rep on venricles');
fh=figure(22);showPatch(GEOM.VER,GEOM.ITRI,ari);set(fh,'Name','Init:ARI on ventricles')

PSIA =lowpassma((GEOM.AMA)*getSmode(T,dep,rep,GEOM.pS,[],4),lpass);
fh=figure(23);clf;set(fh,'Name','Init: BSM');
normphi=norm(PSIREF,'fro');
sigplot(PSIREF,'',GEOM.LAY,1,'b',1,0);
hold on
sigplot(PSIA,'',GEOM.LAY,1,'r',1,0);
drawnow;pause(1);
fh=figure(24);clf;set(fh,'Name','Init:RMS');
plot(rms(PSIREF),'b');hold on;plot(rms(PSIA),'r')
COR=corrcoef(PSIA,PSIREF);
disp(['mean/std rep ' num2str([mean(rep) std(rep)],3) '    cor/rd ' num2str([COR(2,1) norm(PSIREF - PSIA,'fro')/normphi],3)])

