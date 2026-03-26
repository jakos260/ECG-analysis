function [bestfoci,bestdep,notchPot,outp]=multifociscan02(GEOM,leads,scanmode,RVrelVelocity)

% use all foci found each round

% date:30-12-2010
% identification of one or more focal points of depolarization,
% Peter van Dam; 2010 november. All rights reserved Peacs
% global lpass
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones

useAmpl=1;
SCAN.lpass  = 10;
SCAN.VER    = GEOM.VER;
SCAN.ITRI   = GEOM.ITRI;
SCAN.ADJ    = GEOM.ADJ;
SCAN.DIST   = GEOM.DIST;
SCAN.ADJ2W  = GEOM.ADJ2W;
SCAN.DIST2W = GEOM.DIST2W;
SCAN.ADJsurf= GEOM.ADJsurf;
SCAN.DISTsurf = GEOM.DISTsurf;
SCAN.purkinjever  = GEOM.purkinjever;
SCAN.Rpurkinjever = GEOM.Rpurkinjever;
SCAN.Lpurkinjever = GEOM.Lpurkinjever;
SCAN.Rfreewallver = GEOM.Rfreewallver;
SCAN.RendoVER = GEOM.RendoVER;
SCAN.endoVER = GEOM.endoVER;

SCAN.MAX_PURK_SHIFT=0.4;
SCAN.MAX_MYO_SHIFT=0.8;
anisotropicV    = 10 / 33.3333; % wallthickness ventricle / transmural time
SCAN.MAX_MYO_VS = anisotropicV * GEOM.anisotropyRatio;

tic;SCAN.ACT = compAct(SCAN,RVrelVelocity);toc


% prepare
% ECG signals are only used between start QRS and end Twave
SCAN.qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
SCAN.PSIREF=GEOM.BSM(leads,GEOM.specs(2):GEOM.specs(3));
SCAN.rep=100*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1
t=0:size(SCAN.PSIREF,2)-1;
SCAN.T=ones(length(GEOM.VER),1)*t;
SCAN.AMA=zeromean(GEOM.AMA(leads,:));
SCAN.pS = GEOM.pS;

SCAN.normphi=norm(SCAN.PSIREF,'fro');




%% needs to be implemented

if scanmode > 1 &&  scanmode ~= 4
    if useAmpl
        notchPot = sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
    else
        notchPot = 1 - sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
    end
else
    notchPot=zeros(size(GEOM.DIST2W(:,1)));
end

%% Scan for foci
qrsduration=GEOM.specs(3) - GEOM.specs(2) + 1;
dep=[];
bestfoci=[];
outp=[];
bestdep=dep;
startTime=clock;
k=0;
maxshift = 10;
prevrd = 10;
prevCor = -10;
useCor = 1;
while length(unique(bestfoci)) < 15 %isempty(bestdep) || max(bestdep) > SCAN.qrsduration%
    k=k+1;
    [nofocus,dep,foci,bestcor,bestrd,bestshift]=fociscan(SCAN,dep,bestfoci,scanmode,notchPot,RVrelVelocity, useCor);
%     if bestcor > 0.8 && useCor
%         useCor=0;
%     end
    if bestcor - prevCor < 0.005%bestrd > prevrd
        break;
    end
    prevCor = bestcor;
    prevrd = bestrd;
    if ( nofocus ) ||...
            ( ~isempty(outp)&& round(10000*bestcor) <= round(10000*max(outp(:,1))) )
        if maxshift > qrsduration - 80
            break;
        end
        maxshift = maxshift + 15;
    end
    bestdep=dep;
    bestfoci=[bestfoci foci];
    outp=[outp;[bestcor bestrd bestshift max(dep)]];
    nrClust=nrClusters(GEOM,unique(bestfoci));
    disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
        '   QRS duration /sim: ' num2str(qrsduration,3) '  ' num2str(max(dep),3)...
        '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2)...
        '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    figure(13);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci);drawnow
end

%% =======================================================================

function [nofocus,bestdep,foci,bestcor,bestrd,bestshift]=fociscan(SCAN,initdep,initFoci,scanmode,notchPot,maxshift,useCorrelation)

% shifts
if isempty(initdep)
    shift=0*ones(size(SCAN.DIST2W(:,1)));
else
    % initial shifts are MAX_PURK_SHIFT% of the previous depolarization time
    % for purkinje nodes and MAX_MYO_SHIFT % for myocardial nodes. In case
    % the previous depolarization duration is longer than the QRS complex
    % also the purkinje nodeshift is allowed on myocardial nodes. Most
    % probably there is more than one focus in this situation. The reason
    % for making a diffrence in myocardial nodes and purkinje nodes is that
    % under almost all circumcantces the activation wave ends on myocardial
    % nodes. Optimizing the depolarization sequence without changing this
    % latest site of activation causes only local small minima, resulting
    % in an alternative way of inverse timing estimation, without teh
    % posibility to delay activation!
    %     shift(SCAN.purkinjever==1) = min( [ones(sum(SCAN.purkinjever),1) *(SCAN.qrsduration - 85), initdep(SCAN.purkinjever==1) * SCAN.MAX_PURK_SHIFT],[],2);
    shift(SCAN.purkinjever==1) = min( [ones(sum(SCAN.purkinjever),1) *(maxshift), initdep(SCAN.purkinjever==1) * SCAN.MAX_PURK_SHIFT],[],2);
    shift(SCAN.purkinjever==0) = initdep(SCAN.purkinjever==0) * SCAN.MAX_MYO_SHIFT;
    %     shift(SCAN.Rfreewallver==1) = initdep(SCAN.Rfreewallver==1) * SCAN.MAX_MYO_SHIFT;
    shift(shift<1)=0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
% init
foci=0;
nofocus=0;
bestdep=initdep;
bestcor=-1;
bestrd=10;
bestshift=-1;
cors = -1*ones(size(SCAN.VER,1),1);
rds  = 10*ones(size(SCAN.VER,1),1);
deps=zeros(size(SCAN.ADJ));		% all depolarization sequences

donotuse = zeros(length(SCAN.VER),1);
useCorrelation=1;
if ~isempty(initdep)
    % this results in suboptimalisation
    donotuse(shift < 1) = 1;
%     if length(initFoci) <= 3 
%         donotuse(SCAN.RendoVER==1) = 1;
%     elseif rem(length(initFoci),2) == 0
%         donotuse(SCAN.RendoVER==1) = 1;
%         donotuse(SCAN.endoVER==0) = 1;
%     else
%         donotuse(SCAN.endoVER==0) = 1;
%         donotuse(SCAN.RendoVER==0) = 1;
%         donotuse(SCAN.Rfreewallver==1) = 0;
%         useCorrelation=0;
%     end
    
    %     fociDist= min(SCAN.DIST(:,initFoci),[],2);
    % next focus should be at least 35 mm from any focus
    %     donotuse(fociDist < 35) = 1;
    %     % expect for neighbors of focal nodes to enable fine tuning
    %     donotuse(max(SCAN.ADJsurf(:,initFoci),[],2)> 0)=0;
elseif SCAN.qrsduration < 110 % initial focus during sinus rhythm should be in the left chamber
    donotuse(SCAN.Lpurkinjever==0) = 1;
end

if useCorrelation
    for inode=1:length(SCAN.VER)
        % these nodes are already focus and cannot be activated earlier
        if  donotuse(inode) || donotuse(jnode)
            continue;
        end
        dep = min( SCAN.ACT(:,inode),SCAN.ACT(:,jnode));
        if isempty(initdep) % no previous focus determined
            if SCAN.purkinjever(inode)
                dep = dep / (max(dep))* SCAN.qrsduration;
            end
        else
            dep=min(shift(inode) + dep,initdep);
            dep = dep / (max(dep))* SCAN.qrsduration;
        end
        PHI = lowpassma(SCAN.AMA*getSmode(SCAN.T,dep,SCAN.rep,SCAN.pS,notchPot,scanmode),SCAN.lpass);
        COR = corrcoef(PHI,SCAN.PSIREF);
        if cors(inode) < COR(2,1)
            cors(inode)=COR(2,1);
            deps(:,inode) = dep;
            if COR(2,1) == max(cors)
            iinode = inode;
            jjnode = jnode;
            end
            rds(inode) = norm(SCAN.PSIREF - PHI,'fro')/SCAN.normphi;
        end
    end
    % detemine the node with the highest correlation
    % cors(SCAN.purkinjever==0)=cors(SCAN.purkinjever==0)-0.1;
    A=[(1:length(cors))' cors rds];
    if isempty(initdep) && qrsduration < 110
        A(SCAN.Lpurkinjever==0,:)=[];
    elseif isempty(initdep) && qrsduration < 130
        A(SCAN.endoVER==0,:)=[];
    end

%     if isempty(initdep) &&  SCAN.qrsduration < 120
%         cors(SCAN.Rpurkinjever==1 | SCAN.Rfreewallver==1 ) = cors(SCAN.Rpurkinjever==1 |SCAN.Rfreewallver==1 ) - 0.2;
%     end
    A=sortrows(A,2);
    A=A(end:-1:1,:);
    A(A(:,2)<=0,:)=[];
    if isempty(A)
        return;
    else
%         if max(cors) > 0.8
%             A=A(A(:,2)>0.8,:);
%             A=sortrows(A,3);
%         end
        A=A(1,:);
        bestdep=deps(:,A(1,1));
        foci=A(1);
        bestcor=A(2);
        bestshift=shift(A(1));
        PHI =lowpassma(SCAN.AMA*getSmode(SCAN.T,bestdep,SCAN.rep,SCAN.pS,notchPot,scanmode),SCAN.lpass);
        bestrd = norm(SCAN.PSIREF - PHI,'fro')/SCAN.normphi;
%         figure(1);leadv16(SCAN.PSIREF,PHI,'leadsys','nim','max',[-3 3],'paperspeed',200,'do9',1);

        figure(1);clf%leadv16(SCAN.PSIREF,PHI,'leadsys','nim','max',[-3 3],'paperspeed',200,'do9',1);
        sigplot(SCAN.PSIREF,'',GEOM.LAY,1,'b',1,0);
        hold on
        sigplot(PHI,'',GEOM.LAY,1,'r',1,0);
    end
else
    for inode=1:length(SCAN.VER)
        % these nodes are already focus and cannot be activated earlier
        if  donotuse(inode)
            continue;
        end
        dep = SCAN.ACT(:,inode);
        if isempty(initdep) % no previous focus determined
            if SCAN.purkinjever(inode)
                dep = dep / (max(dep))* SCAN.qrsduration;
            end
        else
            dep=min(shift(inode) + dep,initdep);
            dep = dep / (max(dep))* SCAN.qrsduration;
        end
        PHI = lowpassma(SCAN.AMA*getSmode(SCAN.T,dep,SCAN.rep,SCAN.pS,notchPot,scanmode),SCAN.lpass);
        RD = norm(SCAN.PSIREF - PHI,'fro')/SCAN.normphi;
        if rds(inode) > RD
            rds(inode) = RD;
            deps(:,inode) = dep;
        end
    end
    % detemine the node with the highest correlation
    % cors(SCAN.purkinjever==0)=cors(SCAN.purkinjever==0)-0.1;
    A=[(1:length(rds))' rds ];
    A=sortrows(A,2);
    A(A(:,2) >= 10,:)=[];
    if isempty(A)
        return;
    else
        A=A(1,:);
        bestdep=deps(:,A(1,1));
        foci=A(1);
        bestrd=A(2);
        bestshift=shift(A(1));
        PHI =lowpassma(SCAN.AMA*getSmode(SCAN.T,bestdep,SCAN.rep,SCAN.pS,notchPot,scanmode),SCAN.lpass);
        COR = corrcoef(PHI,SCAN.PSIREF);
        bestcor = COR(2,1);
        figure(1);%leadv16(SCAN.PSIREF,PHI,'leadsys','nim','max',[-3 3],'paperspeed',200,'do9',1);
        sigplot(PSIREF,'',GEOM.LAY,1,'b',1,0);
        hold on
        sigplot(PHI,'',GEOM.LAY,1,'r',1,0);
    end
end
    

%% delay scan-----------------------------------------------
% experiment, probably not usefull.
function [bestdep,bestcor,bestrd]=delayscan(GEOM,initdep,leads,scanmode,notchPot)

global lpass;

bestdep=initdep;
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
PSIREF=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
PSIREF=baselinecor(zeromean(PSIREF));
PSIREF=PSIREF(leads,1:qrsduration);
normphi=norm(PSIREF,'fro');

rep=100*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1
t=0:size(PSIREF,2)-1;
T=ones(length(GEOM.VER),1)*t;
AMA=zeromean(GEOM.AMA(leads,:));

PSIA=baselinecor(lowpassma(AMA*getSmode(T,initdep,rep,GEOM.pS,notchPot,scanmode),lpass));
COR=corrcoef(PSIA,PSIREF);
bcor=COR(2,1);
bestrd=norm(PSIREF-PSIA,'fro')/normphi;
%%
keepdep=initdep;
DELAYFACT=0.005;
keepdelay=zeros(size(keepdep));
for i=1:10
    for inode=1:length(GEOM.VER)
        dep=keepdep+keepdelay;
        delay=keepdelay;
        nodes=[inode find(GEOM.DISTsurf(inode,:)<30)];
        delay(nodes)=((-1)^i)*dep(nodes)*DELAYFACT;
        dep=dep+delay;
        dep(dep>qrsduration)=qrsduration;
        dep(dep<0)=0;
        dep(nodes)=max(dep(nodes),0.1*keepdep(nodes));
        dep(nodes)=min(dep(nodes),1.1*keepdep(nodes));
        PSIA=baselinecor(lowpassma(AMA*getSmode(T,dep,rep,GEOM.pS,notchPot,scanmode),lpass));
        COR=corrcoef(PSIA,PSIREF);COR=COR(2,1);
        RD=norm(PSIREF-PSIA,'fro')/normphi;
        if RD<bestrd && round(100*COR)/100>=round(100*bcor)/100
            bestrd=RD;bcor=COR;	bestdep=dep;
            keepdelay(nodes)=dep(nodes)-keepdep(nodes);
        end
    end
end
bestcor=bcor;

%% compute activation from a single node
function ACT = compAct(SCAN,purkVelocity)

area = 50;
ACT= SCAN.DIST2W;
% return;
ADJP = SCAN.ADJ2W;%ADJ / SCAN.MAX_MYO_VS;
ADJPorg= ADJP;
ADJPorg(SCAN.ADJsurf>0) =1000;
ADJPorg(ADJPorg==0) =1000;
ADJP(SCAN.ADJsurf==0)= 0;
for inode =1:length(SCAN.VER)
    if SCAN.purkinjever(inode) %node in the purkinje system area
        mini =find(ADJPorg(inode,:) == min(ADJPorg(inode,:)));
        ADJP(inode,mini) = ADJPorg(inode,mini);
        ADJP(mini,inode) = ADJPorg(inode,mini);
    end
end
ADJP = ADJP / SCAN.MAX_MYO_VS;
relPurkVelocity = purkVelocity / SCAN.MAX_MYO_VS;
% ADJ2P(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1) = SCAN.ADJ(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1)/ SCAN.MAX_MYO_VS;
% ADJP(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1) = SCAN.ADJ2W(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1)/1.5;%RVrelVelocity;

% ADJP(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1) = SCAN.ADJ(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1)/1.75;

% ADJ2W(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1) = SCAN.ADJ(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1)/ SCAN.MAX_MYO_VS;
% ADJ2W(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1) = SCAN.ADJ(SCAN.Rfreewallver==1,SCAN.Rfreewallver==1)/1.1;%RVrelVelocity;
% ADJ2W(SCAN.Rpurkinjever==1,SCAN.Rpurkinjever==1) = ADJ(SCAN.Rpurkinjever==1,SCAN.Rpurkinjever==1) / RVrelVelocity;
% ADJ2W(SCAN.Lpurkinjever==1,SCAN.Lpurkinjever==1) = ADJ(SCAN.Lpurkinjever==1,SCAN.Lpurkinjever==1) / RVrelVelocity;
% 
% ACT = graphdist(ADJ2W);
% 
DISTP= graphdist(ADJP);
hw=waitbar(0,'computing anisotropic purkinje distance matrix');
for inode =1:length(SCAN.VER)
    if SCAN.Lpurkinjever(inode) %node in the purkinje system area
        adj = ADJP;
        indx = find(SCAN.DISTsurf(:,inode) < area & SCAN.Lpurkinjever' ==1);
        adj(inode,indx) = ADJP(inode, indx) / relPurkVelocity;
        adj(indx ,inode) =adj(inode, indx );
        act = graphdistone(adj,inode);
        DISTP(:,inode) = act;
    elseif SCAN.Rpurkinjever(inode) &&~SCAN.Rfreewallver(inode)%node in the purkinje system area
        adj = ADJP;
        indx = find(SCAN.DISTsurf(:,inode) < area & SCAN.Rpurkinjever' ==1);
        adj(inode,indx) = ADJP(inode, indx)/ relPurkVelocity;
        adj(indx ,inode) =adj(inode, indx );
        act = graphdistone(adj,inode);
        DISTP(:,inode) = act;
    end
    waitbar(inode/length(SCAN.VER),hw);
end
close(hw)
% ADJP(SCAN.Rpurkinjever==1,SCAN.Rpurkinjever==1) = ADJP(SCAN.Rpurkinjever==1,SCAN.Rpurkinjever==1) / relPurkVelocity ;
% ADJP(SCAN.Lpurkinjever==1,SCAN.Lpurkinjever==1) = ADJP(SCAN.Lpurkinjever==1,SCAN.Lpurkinjever==1) / relPurkVelocity ;
% DISTP= graphdist(ADJP);
ACT = min(SCAN.DIST2W,DISTP);
