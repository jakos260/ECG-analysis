function [bestfoci,bestdep,outp]=multifociscanAtria(GEOM,usecor,sinkScan)

% use all foci found each round

% date:15032011
% identification of one or more focal points of depolarization,
% Peter van Dam
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones 

SCAN.scanmode= 4;
SCAN.usecor  = 0;
SCAN.sinkScan = 0;
if nargin >= 2
    SCAN.usecor = usecor;
end
if nargin >= 3
    SCAN.sinkScan = sinkScan;
end

SCAN.lpass   = 10;

SCAN.MAX_MYO_VS = 0.8;
SCAN.DELTA_SHIFT=0.2;
SCAN.MAX_MYO_SHIFT = SCAN.DELTA_SHIFT;

if strfind(GEOM.type,'atria')
	SCAN.MAX_MYO_VS=.9;
	SCAN.MAX_MYO_SHIFT=0.8;
end
SCAN.regionR = 30;

SCAN.qrsduration = GEOM.specs(3)-GEOM.specs(2)+1;
%% prepare 
% ECG signals are only used between start QRS and end Twave
if max(rms(GEOM.BSM-lowpassma(GEOM.BSM,20))) < 0.5 % check for pacing spike
    SCAN.PSIREF = baselinecor(GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5)));
else
    difSig=abs(diffrows(rms(GEOM.BSM-lowpassma(GEOM.BSM,20))));
    tonset=find( difSig > 0.5*std(difSig));
    tonset = tonset(1);
    SCAN.qrsduration = tonset-GEOM.specs(2)+1;

    SCAN.PSIREF = baselinecor(GEOM.BSM(:,GEOM.specs(2):tonset));
end
% SCAN.PSIREF=baselinecor(SCAN.PSIREF,1,GEOM.specs(3)-GEOM.specs(2));
SCAN.PSIREForg = SCAN.PSIREF;

SCAN.t=0:(GEOM.specs(3) - GEOM.specs(2))-1;%size(SCAN.PSIREF,2)-1;
SCAN.T = ones(length(GEOM.VER),1)*SCAN.t;


SCAN.AMA = GEOM.AMA;
SCAN.AMAorg = SCAN.AMA;

SCAN.normphi=norm(SCAN.PSIREF,'fro');
SCAN.rep=300*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1

SCAN.L = surflapl(GEOM.VER/1000,GEOM.ITRI,0);
% SCAN.L = GEOM.DIST.^-2;
% SCAN.L(isnan(SCAN.L)) = 0;
% SCAN.L(isinf(SCAN.L)) = 0;
% SCAN.L(SCAN.L < (1/(40^2))) = 0;
% SCAN.L =SCAN.L - diag(sum(SCAN.L));
% SCAN.L = 60 * SCAN.L;
SCAN.DIST2W = GEOM.DIST2W / SCAN.MAX_MYO_VS;
SCAN.ADJ2W  = GEOM.ADJ2W;
% SCAN.ADJ2W = calcAnisADJ(GEOM,GEOM.anisotropyRatio);	

% SCAN.ADJ2W(SCAN.ADJ2W>45)=0;
SCAN.ADJ2W = SCAN.ADJ2W / SCAN.MAX_MYO_VS;

%% Scan for foci
dep       = [];
bestdep   = [];
bestfoci  = [];
bestsinks = [];
outp      = [];
bestcor   = -1;
bestrd    = 10;

startTime=clock;

nrClust = 0;
sinks=[];
while  nrClust < 3 
	[focusfound,dep,foci,bestshift,SCAN]=fociscan(GEOM,SCAN,dep,bestfoci);   
    if SCAN.sinkScan
        [sinkfound,dep,sinks,SCAN]=sinkscan(GEOM,SCAN,dep,[bestfoci foci]);
    end
   
    delay =0;
    [dep,delay] = delayDep(GEOM,SCAN,dep);
    if ~focusfound 
        break;
    end
    maxt = min(SCAN.qrsduration, round(max(dep)) );   
    t= 1:maxt;
    T=ones(length(GEOM.VER),1)*t;   
    PSIA =lowpassma(SCAN.AMAorg*getSmode(T,dep,SCAN.rep,[],[],1),SCAN.lpass);
    COR=corrcoef(PSIA,SCAN.PSIREForg(:,1:size(PSIA,2)));
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREForg(:,1:size(PSIA,2)) - PSIA,'fro')/SCAN.normphi;
    if ( SCAN.usecor && cor <= bestcor + 0.01) ||...
       (~SCAN.usecor && rd >= bestrd - 0.03 ) || ...
            (~isempty(bestdep) && max(dep) > max(bestdep) + 5)
%        ( rd > bestrd && length(GEOM.VER) > 400 ) % coarse geometry makes the rd less trustworthy
        break;	
    end	 	
    bestrd = rd;
    bestcor = cor;
    bestdep=dep;
	outp=[outp;[bestcor bestrd bestshift max(dep)]];
    bestfoci=[bestfoci foci];
    bestsinks=[bestsinks; sinks];
	nrClust=nrClusters(GEOM,unique(bestfoci));
	disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
		  '   QRS duration /sim: ' num2str(SCAN.qrsduration,3) '  ' num2str(diff(range(dep)),3)...
		  '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2) ' ' num2str(delay,2) ...
		  '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    if size(PSIA,1) == size(GEOM.LAY,1)-1
        maxAmpl = (max(max(abs(PSIA))));
        figure(100);clf; sigplot(SCAN.PSIREForg(:,1:size(PSIA,2)),'',GEOM.LAY,1.3/maxAmpl,'b',1,0); 
        hold on; 
        sigplot(PSIA,'',GEOM.LAY,1.3/maxAmpl,'r',1,0);	
    end
    figure(101);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.qrsduration]);drawnow		

    SCAN.leads = 1:size(PSIA,1);
%     if bestrd < 1
%         leadCor=zeros(size(PSIA,1),1);
%         leadRd=zeros(size(PSIA,1),1);
%         for i=1:size(PSIA,1)
%             COR=corrcoef(PSIA(i,:),SCAN.PSIREForg(i,:));
%             leadCor(i) = COR(1,2);
%             leadRd(i) = norm(SCAN.PSIREForg(i,:) - PSIA(i,:),'fro')/norm(SCAN.PSIREForg(i,:),'fro');
%         end
%     %     SCAN.leads(leadRd < bestrd) = []; 
%         SCAN.leads(leadRd > mean(leadRd)) = []; 
%     end
    SCAN.PSIREF = SCAN.PSIREForg(SCAN.leads,:);
    SCAN.AMA = SCAN.AMAorg(SCAN.leads,:);
end


%% =======================================================================
function [focusfound,bestdep,foci,bestshift,SCAN] = fociscan(GEOM,SCAN,initdep,initfoci)

% init
focusfound  = 0;
bestdep     = initdep;
bestshift   = -1;
cors        = -1 * ones(size(GEOM.DIST,1),1);
rds         = 10 * ones(size(GEOM.DIST,1),1);
deps        = zeros(size(GEOM.DIST));		% all depolarization sequences
% shifts
shift=3*ones(size(GEOM.DIST(:,1)));
if ~isempty(initdep)
	% initial shifts are DELTA_SHIFT% of the previous depolarization time 
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
    shift = initdep  * SCAN.DELTA_SHIFT;  
	shift(shift < 1) = 0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
for inode=1:length(GEOM.VER)
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if ~isempty(initdep)
        if shift(inode)==0
            continue;
        end
        dep = shift(inode) + SCAN.DIST2W(:,inode);
        dep=min([dep,initdep],[],2);
    else
        dep = shift(inode) + SCAN.DIST2W(:,inode);
        dep=min([dep,initdep],[],2);
        dep = dep * ( SCAN.qrsduration / max(dep) );         
    end
    rep = dep + 130;
    deps(:,inode) = dep;
    maxt = min(SCAN.qrsduration, round(max(dep)) );   
    t= 1:maxt;
    T=ones(length(GEOM.VER),1)*t;
    PSIA =lowpassma(SCAN.AMA*getSmode(T,dep,rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
    

        COR=corrcoef(PSIA,SCAN.PSIREF(:,1:maxt));
        cors(inode) = COR(2,1);
%     else
        rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA(:,1:maxt),'fro')/SCAN.normphi;
%     end
end
% select focus
A=[(1:length(cors))' cors rds];


if SCAN.usecor
    A(A(:,2)< 0,:)=[];  
    A=sortrows(A,2);
    A=A(end:-1:1,:);
else
    A(A(:,3)>=10,:)=[];  
    A=sortrows(A,3);
end
if ~isempty(A)
    focusfound  = 1;
    A           = A(1,:);
    bestdep     = deps(:,A(1));
    foci        = A(1);
    bestshift   = shift(A(1)); 
end

% if isempty(initdep)
%     SCAN.use = 0;
% end
   
%% =======================================================================
function [sinkfound,bestdep,sink,SCAN]=sinkscan(GEOM,SCAN,initdep,foci)

sinkfound   = 0;
sink       = [];
cors        = -1 * ones(size(GEOM.VER,1),1);
rds         = 10 * ones(size(GEOM.VER,1),1);
deps        = zeros(size(GEOM.ADJ));		% all depolarization sequences
useSurf     = zeros(size(GEOM.DIST,1),1);


% Keep the foci at t=0 otherwise the whole dep sequence starts to shift.
% Leave the first minTFact% (10%) of the activation sequence as is.
T=ones(length(GEOM.VER),1)*(1:max(initdep));
 % surf
PSIA = lowpassma( SCAN.AMA * getSmode(T, initdep, SCAN.rep, [],[],SCAN.scanmode),SCAN.lpass);

if SCAN.usecor
    COR0 = corrcoef(PSIA, SCAN.PSIREF(:,1:size(PSIA,2))); 
    COR0=COR0(2,1);
else
    rd0 = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/SCAN.normphi; 
end


%% scan for sink -----------------------------------------------
tic;
cors = ones(size(initdep))*-1;
for inode=1:length(GEOM.VER)
    if cors(inode) == -1
        ADJ = SCAN.ADJ2W;
        neighs = find(GEOM.DIST(inode,:) < 25 );
        ADJ(neighs,:) = ADJ(neighs,:) * SCAN.prolongDistFact;
        ADJ(:,neighs) = ADJ(:,neighs) * SCAN.prolongDistFact;
        deps(:,inode) = focipaths(ADJ,foci,initdep(foci));
%         if max(deps(:,inode)) > 1.05 * SCAN.qrsduration 
%             continue
%         end
    %     dep = dep - min(dep);
    %     dep  = dep  *( SCAN.qrsduration / max(dep) );
        T=ones(length(GEOM.VER),1)*(1:max(deps(:,inode)));
        PSIA = lowpassma( SCAN.AMA*getSmode(T, deps(:,inode), SCAN.rep, [],[],SCAN.scanmode),SCAN.lpass);

        if SCAN.usecor
            COR = corrcoef(PSIA, SCAN.PSIREF(:,1:size(PSIA,2))); 
    %         cors(inode) = max(cors(neighs),ones(length(neighs),1)*COR(2,1));      
            cors(inode) = COR(2,1);      
            if cors(inode) < COR0  && ...
               min(GEOM.DISTsurf(inode,foci)) > 35
               candidates = find(cors > COR0);
               if isempty(candidates) || ...
                  min(GEOM.DISTsurf(inode,candidates)) > 35
                    cors( GEOM.DISTsurf(inode,:) < 35 ) = cors(inode);              
               end
            end
        else
            rd = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/SCAN.normphi; 
%             if rds(inode) < rd0 && ...
%                min(GEOM.DISTsurf(inode,foci)) > 35
%                candidates = find(rds > COR0);
%                if isempty(candidates) || ...
%                   min(GEOM.DISTsurf(inode,candidates)) > 35
%                     cors( GEOM.DISTsurf(inode,:) < 35 ) = cors(inode);              
%                end
%             end
    %         rds(neighs) = min(rds(neighs),ones(length(neighs),1)*rd);
            rds(inode) = rd;
        end
    end
end
toc
%%
if SCAN.usecor
    sinkfound = any(cors > COR0);
    if sinkfound 
        sink = find(cors==max(cors));
    else
        sink =[];
    end
%     sel = find((cors-COR0) > 0.8*max(cors-COR0));
%     selection =[];
%     for i=1:length(sel)
%         selection = [selection find(GEOM.DISTsurf(sel(i),:) < 25)];
%     end
%     selection = unique(selection);
%     selection = find(cors-COR0 > 0.8*max(cors-COR0) );
else
    sinkfound = any(rds < rd0);
    if sinkfound 
        sink = find(rds == min(rds));
    else
        sink =[];
    end
end
if sinkfound
    bestdep  = deps(:,sink) ;% *( SCAN.qrsduration / max(dep) );    
    neighs = find(GEOM.DISTsurf(sink,:) < 25 );
    ADJ = SCAN.ADJ2W;
    ADJ(neighs,:) = ADJ(neighs,:) * SCAN.prolongDistFact;
    ADJ(:,neighs) = ADJ(:,neighs) * SCAN.prolongDistFact;
    SCAN.ADJ2W = ADJ;
else
    bestdep = initdep;
end
return
% sink = selection;
% sinkfound = ~isempty(sink);
% ADJ(selection,:) = ADJ(selection,:) * SCAN.prolongDistFact;
% ADJ(:,selection) = ADJ(:,selection) * SCAN.prolongDistFact;
% dep = focipaths(ADJ,foci,initdep(foci));
% dep = dep - min(dep);


% PSIA = lowpassma( SCAN.AMA*getSmode(T, bestdep, SCAN.rep, GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
%%
if SCAN.usecor
    COR = corrcoef(PSIA, SCAN.PSIREF(:,mint:SCAN.qrsduration)); 
    if COR(2,1) > COR0
        SCAN.ADJ2W = ADJ;
    else
        bestdep=initdep;
    end
else
    rd = norm(SCAN.PSIREF(:,mint:SCAN.qrsduration) - PSIA,'fro')/SCAN.normphi; 
    if rd < rd0
        SCAN.ADJ2W = ADJ;
    else
        bestdep=initdep;
    end

end

return


for inode=1:length(GEOM.VER)
    fixed=[];
    % only delay parts that are already later
    if initdep(inode) > SCAN.qrsduration * minTFact
        
        fixedNodesSurf = find(GEOM.DISTsurf(:,inode)> SCAN.regionR  | initdep < SCAN.qrsduration * minTFact);
        depSurf = adaptdep(SCAN,initdep,fixedNodesSurf);
       
        minDist = min(GEOM.DIST(inode,GEOM.ADJsurf(inode,:)==0 & GEOM.ADJ(inode,:)~= 0 ));
        if ~isempty(minDist)
            fixedNodesTrans = find( ( GEOM.DIST(:,inode)   > SCAN.regionR ) | ...
                                      initdep < SCAN.qrsduration * minTFact);                                  
            depTrans = adaptdep(SCAN,initdep,fixedNodesTrans);
        else
            fixedNodesTrans = fixedNodesSurf;
        end
        depSurf = depSurf * (SCAN.qrsduration / max(depSurf) );
        depTrans= depTrans * (SCAN.qrsduration / max(depTrans) );
        depTrans(depTrans < initdep) = initdep(depTrans < initdep) ;
        depSurf(depSurf < initdep) = initdep(depSurf < initdep) ;

        % surf
        PSIA = lowpassma( SCAN.AMA*getSmode(T, depSurf, SCAN.rep, GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
        
        if SCAN.usecor
            CORSu = corrcoef(PSIA, SCAN.PSIREF(:,mint:SCAN.qrsduration)); 
        else
            rdSu = norm(SCAN.PSIREF(:,mint:SCAN.qrsduration) - PSIA,'fro')/SCAN.normphi; 
        end
        
        %trans
        PSIA = lowpassma( SCAN.AMA*getSmode(T, depTrans, SCAN.rep, GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
        if SCAN.usecor
            CORTr = corrcoef(PSIA, SCAN.PSIREF(:,mint:SCAN.qrsduration)); 
            cors(inode) = max(CORSu(2,1),CORTr(2,1));            
            if CORSu(2,1) > CORTr(2,1)
                useSurf(inode) = 1;
                fixed = fixedNodesSurf;
            else
                fixed = fixedNodesTrans;
            end
        else
            rdTr = norm(SCAN.PSIREF(:,mint:SCAN.qrsduration) - PSIA,'fro')/SCAN.normphi; 
            rds(inode) = min(rdSu,rdTr);
            if rdSu < rdTr
                useSurf(inode) = 1;
                fixed = fixedNodesSurf;
            else
                fixed = fixedNodesTrans;
            end
        end

        if ( SCAN.usecor && cors(inode) == CORSu(2,1)) || ...
            ~SCAN.usecor && rds(inode) == rdSu
            deps(:,inode) = depSurf;                     
        else
            deps(:,inode) = depTrans;
        end
      
    end
end
%% select the best sink
A=[(1:length(cors))' cors rds ]; 
if SCAN.usecor
    A(A(:,2) <= COR0,:)=[];	
    A=sortrows(A,2);
    A=A(end:-1:1,:);
else
    A(A(:,3)>=rd0,:)=[];  
    A=sortrows(A,3);
end
if ~isempty(A) && length(cors) - length(fixed) > 10 % at least 10 nodes are affected
    sinkfound = 1;
    A=A(1,:);
    bestdep=deps(:,A(1));
    sink = A(1);
    fixed = find(GEOM.DISTsurf(:,sink)> SCAN.regionR  | initdep < SCAN.qrsduration * minTFact);
    if ~useSurf(sink)
        minDist = min(GEOM.DIST(inode,GEOM.ADJsurf(sink,:)==0 & GEOM.ADJ(sink,:)~= 0 ));
        if ~isempty(minDist)
            fixed = find( ( GEOM.DIST(:,sink)   > SCAN.regionR ) | initdep < SCAN.qrsduration * minTFact);                                  
        end
    end
    SCAN = adaptScanDist(GEOM,SCAN,fixed);   
else
    bestdep = initdep;
end
   
%% =======================================================================
function [bestdep,delay] = delayDep(GEOM,SCAN,initdep)


delays= -2:0.5:2;
cors = ones(size(delays)) *-1;
rds  = ones(size(delays)) *10;

for i= 1:length(delays)
    delay= delays(i);
    maxt = min(SCAN.qrsduration, round(max(initdep) +delay) );   
    t= 1:maxt;    T=ones(length(GEOM.VER),1)*t;
    PSIA =lowpassma(SCAN.AMA*getSmode(T,initdep+delay,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
    
    COR=corrcoef(PSIA,SCAN.PSIREF(:,1:maxt));    
    cors(i) = COR(2,1);
    rds(i) = norm(SCAN.PSIREF(:,1:maxt) - PSIA(:,1:maxt),'fro')/SCAN.normphi;
end

if SCAN.usecor    
    delay = delays(find(cors==max(cors)));     
else
    delay = delays(find(rds == min(rds)));  
end
bestdep = initdep + delay;
