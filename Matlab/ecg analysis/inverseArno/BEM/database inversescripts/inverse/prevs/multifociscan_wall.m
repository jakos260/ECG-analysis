function [bestfoci,bestdep,outp]=multifociscan_wall(varargin)

% use all foci found each round
% date:05102012
% identification of one or more focal points of depolarization,
% Peter van Dam
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones 


GEOM = varargin{1};
usecor = 1;
sinkscan = 0;
velocity = Inf;
isSinus = 0;
clusters =5;
pp=2;
if length(varargin) < 1
	error('This routine needs at least one parameter');
else    
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'usecor'
                    usecor=varargin{pp+1};pp=pp+2;
                case 'sinkscan'
                    sinkscan=varargin{pp+1};pp=pp+2;
                case 'velocity'
                    velocity=varargin{pp+1};pp=pp+2;								
                case 'issinus'
                    isSinus=varargin{pp+1};pp=pp+2;	
                case 'clusters'
                    clusters = max(varargin{pp+1},1);pp=pp+2;	
                otherwise
                    error('unknown parameter');
            end
        end
    end
end



SCAN.scanmode= 1;
SCAN.usecor  = usecor;
SCAN.sinkScan = sinkscan;

SCAN.lpass   = 10;
SCAN.qrsduration = GEOM.specs(3)-GEOM.specs(2)+1;
SCAN.MAX_MYO_VS = velocity;
SCAN.isSinus = isSinus;
if isSinus
    SCAN.MAX_MYO_VS = 0.8;
end

SCAN.DELTA_SHIFT=0.35;
SCAN.MAX_MYO_SHIFT = 0.5;
if SCAN.qrsduration > 120    
    SCAN.prolongDistFact = 1.75; % might have ischemica areas (dead tissue expected)
else
    SCAN.prolongDistFact = 1.25; % only mild adapations, no dead tissue expected
end
if strfind(GEOM.type,'atria')
	SCAN.MAX_MYO_SHIFT=0.8;
end
SCAN.regionR = 30;


%% prepare 
% ECG signals are only used between start QRS and end Twave

SCAN.PSIREF = GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
SCAN.PSIREF = baselinecor(SCAN.PSIREF,1,GEOM.specs(3)-GEOM.specs(2));
SCAN.PSIREF = [zeros(size(SCAN.PSIREF,1),10) SCAN.PSIREF];
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
SCAN.ADJ2W  = GEOM.ADJ2W;
SCAN.DIST2W = GEOM.DIST2W;
if ~isinf(SCAN.MAX_MYO_VS)
	SCAN.DIST2W = GEOM.DIST2W / SCAN.MAX_MYO_VS;
    SCAN.ADJ2W = SCAN.ADJ2W / SCAN.MAX_MYO_VS;
end


%% Scan for foci
dep       = [];
bestdep   = [];
bestfoci  = [];
bestsinks = [];
outp      = [];
bestcor  = -1;
bestrd = 10;
prevCor = -1;
prevrd  = 10;
startTime=clock;

nrClust = 0;
sinks=[];

while  nrClust < clusters 
    [focusfound,dep,foci,bestshift,SCAN]=fociscan(GEOM,SCAN,dep,bestfoci);   
    if SCAN.sinkScan
        [sinkfound,dep,sinks,SCAN]=sinkscan(GEOM,SCAN,dep,[bestfoci foci]);
    end
    delay =0;
%     [dep,delay] = delayDep(GEOM,SCAN,dep);
    if ~focusfound
        break;
    end
    PSIA =lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * (1:max(SCAN.qrsduration, round(max(dep)) )),dep,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);    
    COR=corrcoef(PSIA,SCAN.PSIREF(:,1:size(PSIA,2)));
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
    if ( SCAN.usecor && cor <= bestcor + 0.01) ||...
       (~SCAN.usecor && rd >= bestrd - 0.03 )
        break;
    end
    
    prevCor = bestcor;
    prevrd = bestrd;
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
        maxAmpl = round(max(max(abs(PSIA))));
        figure(100);clf; sigplot(SCAN.PSIREF(:,1:size(PSIA,2)),'',GEOM.LAY,1.3/maxAmpl,'b',1,0); 
        hold on; 
        sigplot(PSIA,'',GEOM.LAY,1.3/maxAmpl,'r',1,0);	
    end
    figure(101);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.qrsduration]);drawnow		
end

%% =======================================================================
function [focusfound,bestdep,foci,bestshift,SCAN] = fociscan(GEOM,SCAN,initdepIn,initfoci)

% init
focusfound  = 0;
bestdep     = initdepIn;
bestshift   = -1;
cors        = -1 * ones(size(GEOM.DIST,1),1);
rds         = 10 * ones(size(GEOM.DIST,1),1);
deps        = zeros(size(GEOM.DIST));		% all depolarization sequences

% shifts
shift= 10*ones(size(GEOM.DIST(:,1)));
initdep = initdepIn;
if  ~isempty(initdep)
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
    shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) * SCAN.DELTA_SHIFT;  
    shift(GEOM.purkinjever==0) = initdep(GEOM.purkinjever==0) * SCAN.MAX_MYO_SHIFT;
    if SCAN.isSinus
        shift(GEOM.purkinjever==1) = min(30, shift(GEOM.purkinjever==1) );
    end  
	shift(shift < 1) = 0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
for inode=1:length(GEOM.VER)
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if SCAN.isSinus && ~GEOM.purkinjever(inode)
        continue
    end
    if any(initfoci==inode) 
        continue;
    end
    dep = min([shift(inode) + SCAN.DIST2W(:,inode),initdep],[],2);
	% if the velocity is undefined scale depolarization by the QRS width
    if isinf(SCAN.MAX_MYO_VS)
        dep = dep * (SCAN.qrsduration / max(dep) );
    end

%     PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:max(SCAN.qrsduration, round(max(dep)) )),dep,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
% 	% PETER OOSTERHOFF  PETER OOSTERHOFF PETER OOSTERHOFF PETER OOSTERHOFF PETER OOSTERHOFF
% 	% compute the delay/shit in time based on the xcorr signals. Maybe you should use the 
% 	% all twelve lead signals to allign (GEOM.v12 contain teh 12 lead indices).  
% 	shitXcorr = ????
%     dep = dep + shitXcorr;
    PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:max(SCAN.qrsduration, round(max(dep)) )),dep,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
    COR=corrcoef(PSIA,SCAN.PSIREF(:,1:size(PSIA,2)));
    cors(inode) = COR(2,1);
    rds(inode) = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
    deps(:,inode) = dep;
end

disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)] )])

% below a possible refinement to plac the focus somewhat with in the wall
% currently switche off, because is it requires too much computation time
if 1
    minCor = max(cors) - 0.25 * std(cors);
    for inode=1:length(GEOM.VER)
        if cors(inode) < minCor
            continue
        end
        if SCAN.isSinus && ~GEOM.purkinjever(inode)
            continue
        end
        if any(initfoci==inode)  
            continue;
        end

        % foci should be found in the first part of the QRS complex, if shift
        % is 0 it cannot be adpated anymore
        if cors(inode) > minCor 
            cor0 = cors(inode);
            rd0 = rds(inode);
        
            ADJ = SCAN.ADJ2W;
            b = find(GEOM.DIST(inode,:) < 20);
            ADJ(inode,b) = ADJ(inode,b) * 0.75;
            ADJ(b,inode) = ADJ(b,inode) * 0.75;
            dep1 = focipaths(ADJ,[initfoci inode],[initdep(initfoci); shift(inode)]);
            if isinf(SCAN.MAX_MYO_VS)
                dep1 = dep1 * (SCAN.qrsduration/ max(dep1) );
            end          
            maxt = max(SCAN.qrsduration, round(max(dep1)) );
            if maxt >300
                %TODO ????
                continue
            end
            t = 1:maxt;   
            T=ones(length(GEOM.VER),1) * t;
            PSIA =lowpassma(SCAN.AMA*getSmode(T,dep1,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);

            COR=corrcoef(PSIA,SCAN.PSIREF(:,1:maxt));
            cor1 = COR(1,2);
            rd1 = norm(SCAN.PSIREF(:,1:maxt) - PSIA(:,1:maxt),'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
            cor2 = cor1;
            rd2 = rd1;
%             if cor1 > cor0 ||  ( isempty(initdep) && min(rds) < 0.7 )
%                 ADJ = SCAN.ADJ2W;
%                 ADJ(inode,b) = ADJ(inode,b) * 0.5;
%                 ADJ(b,inode) = ADJ(b,inode) * 0.5;
%                 dep2 = focipaths(ADJ,[initfoci inode],[initdep(initfoci); shift(inode)]);
%                 if isinf(SCAN.MAX_MYO_VS)
%                     dep2 = dep2 * ((SCAN.qrsduration+10) / max(dep2) );
%                 end
% 
%                 maxt = max(SCAN.qrsduration, round(max(dep2)) );
%                 t= 1:maxt;
%                 T=ones(length(GEOM.VER),1)*t;
%                 PSIA =lowpassma(SCAN.AMA*getSmode(T,dep2,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
% 
%                 COR=corrcoef(PSIA,SCAN.PSIREF(:,1:maxt));
%                 cor2 = COR(1,2);
%                 rd2 = norm(SCAN.PSIREF(:,1:maxt) - PSIA(:,1:maxt),'fro')/SCAN.normphi;
%             end
            if SCAN.usecor
                if cor0 > cor1 && cor0 > cor2

                elseif cor1 >= cor2
                    cors(inode) = cor1;
                    deps(:,inode) = dep1;
                else
                    cors(inode) = cor2;
                    deps(:,inode) = dep2;
                end
            else
                if rd0 < rd1 && rd0 < rd2
                    
                elseif rd1 <= rd2
                    rds(inode)= rd1;
                    deps(:,inode) = dep1;
                else
                    rds(inode)= rd2;
                   deps(:,inode) = dep2;
                end        
            end
        end
    end
end
disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)] )])
if isfield(GEOM,'typ')
    cors(GEOM.typ>3)=0;
    rds(GEOM.typ>3)=-10;
end
% select focus
A=[(1:length(cors))' cors rds];
if ~SCAN.usecor 
    A(A(:,3)>=10,:)=[];  
    A=sortrows(A,3);
else
    A(A(:,2)< 0,:)=[];
    A=sortrows(A,2);
    A=A(end:-1:1,:);
    
end
if ~isempty(A)
    focusfound  = 1;
    A           = A(1,:);
    bestdep     = deps(:,A(1));
    foci        = A(1);
    bestshift   = shift(A(1)); 
end
   
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
PSIA = lowpassma( SCAN.AMA * getSmode(T, initdep, SCAN.rep, GEOM.pS,[],SCAN.scanmode),SCAN.lpass);

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
        PSIA = lowpassma( SCAN.AMA*getSmode(T, deps(:,inode), SCAN.rep, GEOM.pS,[],SCAN.scanmode),SCAN.lpass);

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
    maxt = max(SCAN.qrsduration, round(max(initdep) +delay) );   
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


%% 
function dep = adaptdep(SCAN,initdep,fixed)

% from=getfrom(GEOM.ADJ,initdep);
% ADJ = SCAN.ADJ2W;
ADJSLOW = SCAN.ADJ2W * (SCAN.prolongDistFact); 
done=zeros(size(initdep));
% normals= done;
% normals(fixed) = 1;
% done(initdep < min(initdep(~normals))) = 1;
done(fixed)=1;

dep=initdep;
volgorde = [(1 : length(done) )' dep];
volgorde = volgorde(done==0,:);
volgorde = sortrows(volgorde,2);
for j=1:length(volgorde)
    i=volgorde(j,1);
%     neighNorm = find(ADJSLOW(:,i) > 0 & done==1 & normals==1);
%     neighSlow = find(ADJSLOW(:,i) > 0 & done==1 & normals==0);
%     dep(i) = min([dep(neighNorm) + ADJ(neighNorm,i); dep(neighSlow) + ADJSLOW(neighSlow,i)]);

    neighSlow = find(ADJSLOW(:,i) & done );
    if ~isempty(neighSlow)
        dep(i) = min(dep(neighSlow) + ADJSLOW(neighSlow,i));
    else
        a=1;
    end
    done(i)=1;
end
%% This is a cruede appoximation, because teh velocity outside the adapted
% zone is normal. For the addition of a focus this is reletive
% unimportant, because once started in teh ischemic zone it will not affect
% much outside this area. For teh initial estimate this is good enough

function SCAN = adaptScanDist(GEOM,SCAN,fixed)

% TODO make more accurate save an interpolation in the foci search
a=[1:length(SCAN.DIST2W)];
a(fixed)=[];

SCAN.DIST2W(a,:) = SCAN.DIST2W(a,:) *(SCAN.prolongDistFact);
SCAN.DIST2W(:,a) = SCAN.DIST2W(:,a) *(SCAN.prolongDistFact);
SCAN.ADJ2W(a,:) = SCAN.ADJ2W(a,:) *(SCAN.prolongDistFact);
SCAN.ADJ2W(:,a) = SCAN.ADJ2W(:,a) *(SCAN.prolongDistFact);
inode = 1;
dep = focipaths(SCAN.ADJ2W,162,0);
% showPatch(GEOM.VER,GEOM.ITRI,dep,'nodes',inode)

%%
function interDep = iterpolatedDep(SCAN,dep,nodes2)

nn=length(dep);
n2=length(nodes2);
n1=nn-n2;
nodes1=zeros(1,n1);
index=zeros(1,nn);
index(nodes2)=1;

% create submatrices
L1=zeros(nn,nn-n2);
L2=zeros(nn,n2);
% form L1 and L2
k=0;
l=0;

for j=1:nn,
  if index(j) == 0
    k=k+1;
    L1(:,k) = SCAN.L(:,j);
    nodes1(k)=j;
  else
    l=l+1;
    L2(:,l) = SCAN.L(:,nodes2(l));
  end
end

%compute the interpolating matrix
INT=-inv(L1'*L1)*L1'*L2;

interDep = dep;
interDep(nodes1) = INT * dep(nodes2);
