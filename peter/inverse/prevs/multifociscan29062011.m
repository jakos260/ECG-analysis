function [bestfoci,bestdep,outp]=multifociscan(GEOM,leads)

% use all foci found each round

% date:15032009
% identification of one or more focal points of depolarization,
% Peter van Dam
% global lpass
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones 

SCAN.scanmode= 1;
SCAN.usecor  = 1;
SCAN.leads   = leads;
SCAN.lpass   = 5;
%% anisotropic matrix calculation
% GEOM.ADJ2W=calcAnisADJ(GEOM,anisotropyRatio);	
SCAN.MAX_MYO_VS=0.7;
SCAN.DELTA_SHIFT=0.3;
SCAN.MAX_MYO_SHIFT=0.9;
SCAN.reduceVeloFact = 0.75;
if strfind(GEOM.type,'atria')
	SCAN.MAX_MYO_VS=.9;
	SCAN.MAX_MYO_SHIFT=0.8;
end
SCAN.DIST2W = GEOM.DIST2W / SCAN.MAX_MYO_VS;

%% prepare 
% ECG signals are only used between start QRS and end Twave
SCAN.qrsduration = GEOM.specs(3)-GEOM.specs(2)+1;
SCAN.PSIREF = GEOM.BSM(leads,GEOM.specs(2):GEOM.specs(5));
SCAN.PSIREF = SCAN.PSIREF(:,1:SCAN.qrsduration+1);
SCAN.PSIREF=zeromean(baselinecor(SCAN.PSIREF));
SCAN.PSIREForg = SCAN.PSIREF;

SCAN.t=0:size(SCAN.PSIREF,2)-1;
SCAN.T=ones(length(GEOM.VER),1)*SCAN.t;


SCAN.AMA = zeromean(GEOM.AMA(leads,:));
SCAN.AMAorg = SCAN.AMA;

SCAN.normphi=norm(SCAN.PSIREF,'fro');
SCAN.rep=300*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1

% SCAN.L = surflapl(GEOM.VER/1000,GEOM.ITRI,0);
% SCAN.L = GEOM.DIST.^-2;
% SCAN.L(isnan(SCAN.L)) = 0;
% SCAN.L(isinf(SCAN.L)) = 0;
% SCAN.L(SCAN.L < (1/(40^2))) = 0;
% SCAN.L =SCAN.L - diag(sum(SCAN.L));
% SCAN.L = 60 * SCAN.L;

%% Scan for foci
dep       = [];
bestfoci  = [];
bestsinks = [];
outp      = [];
bestcor   = -1;
bestrd    = 10;

startTime=clock;

nrClust = 0;
sinks=[];

while  nrClust < 6
	[focusfound,dep,foci,bestshift,SCAN]=fociscan(GEOM,SCAN,dep);
    % find first 3 foci for sinus rythm
    if SCAN.qrsduration < 110 && length(bestfoci) < 2 
%         if length(bestfoci) == 2 
%             [sinkfound,dep,sinks,SCAN]=sinkscan(GEOM,SCAN,dep);    
%         end
    else
        [sinkfound,dep,sinks,SCAN]=sinkscan(GEOM,SCAN,dep);
    end
    if ~focusfound
        break;
    end
    PSIA =lowpassma(SCAN.AMAorg*getSmode(SCAN.T,dep,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);
    COR=corrcoef(PSIA,SCAN.PSIREForg);
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREForg - PSIA,'fro')/SCAN.normphi;
    if ( SCAN.usecor && cor <= bestcor + 0.005) ||...
       (~SCAN.usecor && rd >= bestrd - 0.005 ) || ...
       ( rd > bestrd && length(GEOM.VER) > 400 ) % coarse geometry makes the rd less trustworthy
        break;	
    end	 	
    bestrd = rd;
    bestcor = cor;
    bestdep=dep;
	outp=[outp;[bestcor bestrd bestshift max(dep)]];
    bestfoci=[bestfoci foci];
    bestsinks=[bestsinks sinks];
	nrClust=nrClusters(GEOM,unique(bestfoci));
	disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
		  '   QRS duration /sim: ' num2str(SCAN.qrsduration,3) '  ' num2str(max(dep),3)...
		  '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2)...
		  '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    
    figure(100);clf; sigplot(SCAN.PSIREForg,'',GEOM.LAY,1,'b',1,0); hold on; sigplot(PSIA,'',GEOM.LAY,1,'r',1,0);
	figure(101);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks);drawnow		

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
function [focusfound,bestdep,foci,bestshift,SCAN] = fociscan(GEOM,SCAN,initdep)

% init
focusfound  = 0;
bestdep     = initdep;
bestshift   = -1;
cors        = -1 * ones(size(GEOM.DIST,1),1);
rds         = 10 * ones(size(GEOM.DIST,1),1);
deps        = zeros(size(GEOM.DIST));		% all depolarization sequences
VS          = zeros(size(GEOM.DIST));
% shifts
shift=zeros(size(GEOM.DIST(:,1)));
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
% 	if abs(max(initdep)-qrsduration)<1e-10
%     shift = initdep;
%     shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) - 40;
%     shift(GEOM.purkinjever==0) = initdep(GEOM.purkinjever==0) - 10;
%     shift(shift < 0) = 0;
    shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) * SCAN.DELTA_SHIFT;
    shift(GEOM.purkinjever==0) = initdep(GEOM.purkinjever==0) * SCAN.MAX_MYO_SHIFT;
%         shift=initdep*DELTA_SHIFT;
% 	else
%  		shift=initdep*DELTA_SHIFT;
% 	end
	shift(shift<1)=0; % prevent suboptimalization
end
maxt = round(SCAN.qrsduration * 0.9);     
T= SCAN.T(:,1:maxt);

%% scan for foci -----------------------------------------------
for inode=1:length(GEOM.VER)
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if ~isempty(initdep) && (shift(inode)==0 || shift(inode) > 0.4 * SCAN.qrsduration)
        continue;
    end   
	dep = shift(inode) + SCAN.DIST2W(:,inode);
    if ~isempty(initdep) 
        depTmp = min(dep,initdep);				
        % choose scale such that the total activation time remains constant
        dep = dep * ( SCAN.qrsduration / max(depTmp) );  
        dep=min(dep,initdep);	
%         fixedNodes = find( ( GEOM.DIST(:,inode)   > 35 ) | ...
%                              initdep > SCAN.qrsduration * 0.5 |...
%                              initdep < initdep(inode) );
        dep(GEOM.DIST(:,inode)   < 35 & initdep > initdep(inode) * 1.35) = ...
        initdep(GEOM.DIST(:,inode)   < 35 & initdep > initdep(inode) * 1.35);
%         dep = iterpolatedDep(SCAN,dep,fixedNodes);
    else
        dep = dep * ( SCAN.qrsduration / max(dep) );
        
        VS(inode) = SCAN.qrsduration / max(dep);
    end    

    deps(:,inode) = dep;
    PSIA =lowpassma(SCAN.AMA*getSmode(T,dep,SCAN.rep,GEOM.pS,[],SCAN.scanmode),SCAN.lpass);

    if SCAN.usecor
        COR=corrcoef(PSIA,SCAN.PSIREF(:,1:maxt));
        cors(inode) = COR(2,1);
    else
        rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA(:,1:maxt),'fro')/SCAN.normphi;
    end
   
end
%%
A=[(1:length(cors))' cors rds VS];
if isempty(initdep) && SCAN.qrsduration < 110
    A(GEOM.Lpurkinjever==0,:)=[];
elseif isempty(initdep) && SCAN.qrsduration < 140
    A(GEOM.endoVER==0,2) = A(GEOM.endoVER==0,2) - 0.1;
    A(GEOM.endoVER==0,3) = A(GEOM.endoVER==0,3) + 0.1;
end
%     A(A(:,3)>=min(A(:,3))+0.05,:)=[];  
% 	if abs(A(1,2) - mean(A(1:3,2))) < 1e-10
% 		return;
% 	end
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
%     if isempty(initdep)
%         vs = A(end);
%         from=getfrom(GEOM.ADJ,bestdep);
%         for i=2:length(from)
%             SCAN.DIST2W(from(i-1),from(i)) = SCAN.DIST2W(from(i-1),from(i)) *vs;
%             SCAN.DIST2W(from(i),from(i-1)) = SCAN.DIST2W(from(i),from(i-1)) *vs;
%         end
%     end
end

   
%% =======================================================================
function [sinkfound,bestdep,sinks,SCAN]=sinkscan(GEOM,SCAN,initdep)

sinkfound   = 0;
cors        = -1 * ones(size(GEOM.VER,1),1);
rds         = 10 * ones(size(GEOM.VER,1),1);
deps        = zeros(size(GEOM.ADJ));		% all depolarization sequences
% initdeprev  = SCAN.qrsduration - initdep;
% 
% DIST2W= SCAN.DIST;%/ (SCAN.MAX_MYO_VS * 1.2 );
% sinks=[];
% for i=1:length(GEOM.VER)
% 	if all(initdep(GEOM.DIST(:,i)~=0 & GEOM.DIST(:,i)<25) < initdep(i) )
% 		sinks=[sinks i];
% 	end
% end


% sinks = find(initdep == max(initdep));
doScan = ones(size(initdep));
% regionRendo = 40;
regionR = 40;
mint = 10;
% adaptFact =0.35;
minTFact = 0.1; %round(SCAN.qrsduration * 0.3)
T = SCAN.T(:,mint:SCAN.qrsduration);
 % surf
PSIA = lowpassma( SCAN.AMA*getSmode(T, initdep, SCAN.rep, GEOM.pS,[],SCAN.scanmode),SCAN.lpass);

if SCAN.usecor
    COR0 = corrcoef(PSIA, SCAN.PSIREF(:,mint:SCAN.qrsduration)); 
    COR0=COR0(2,1);
else
    rd0 = norm(SCAN.PSIREF(:,mint:SCAN.qrsduration) - PSIA,'fro')/SCAN.normphi; 
end
        

%% scan for sinks -----------------------------------------------
for inode=1:length(GEOM.VER)
    fixed=[];
    % only delay parts that are already later
    if doScan(inode) && initdep(inode) > SCAN.qrsduration * minTFact
%         d = initdep;
% 
%         delta = ((SCAN.qrsduration - initdep(inode)) * adaptFact);
%         d(inode) = min(SCAN.qrsduration, initdep(inode) + delta );
       
        fixedNodesSurf = find(GEOM.DISTsurf(:,inode)> regionR  | initdep < SCAN.qrsduration * minTFact);
        tmpdepS = adaptdep(GEOM,SCAN,initdep,fixedNodesSurf);
%         fixedNodesSurf = unique([fixedNodesSurf; inode]);        
        
        minDist = min(GEOM.DIST(inode,GEOM.ADJsurf(inode,:)==0 & GEOM.ADJ(inode,:)~= 0 ));
        if ~isempty(minDist)
%             opsNode  = find(GEOM.DIST(inode,:) == minDist );
%             fixedNodesTrans = find( ( GEOM.DISTsurf(:,opsNode) > regionR  & ...
%                                       GEOM.DISTsurf(:,inode)   > regionR  ) | ...
%                                       initdep < SCAN.qrsduration * minTFact);
            fixedNodesTrans = find( ( GEOM.DIST(:,inode)   > regionR ) | ...
                                      initdep < SCAN.qrsduration * minTFact);                                  
%             d(opsNode) = min(SCAN.qrsduration, initdep(opsNode) + 0.9*delta );
            tmpdepT = adaptdep(GEOM,SCAN,initdep,fixedNodesTrans);
%             fixedNodesTrans = unique([fixedNodesTrans; inode; opsNode]);
        else
            fixedNodesTrans = fixedNodesSurf;
        end
        depSurf = tmpdepS;
        depTrans = tmpdepT;
%         depSurf = iterpolatedDep(SCAN,d,fixedNodesSurf);
%         depTrans = iterpolatedDep(SCAN,d,fixedNodesTrans);
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
                fixed = fixedNodesSurf;
            else
                fixed = fixedNodesTrans;
            end
        else
            rdTr = norm(SCAN.PSIREF(:,mint:SCAN.qrsduration) - PSIA,'fro')/SCAN.normphi; 
            rds(inode) = min(rdSu,rdTr);
            if rdSu < rdTr
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
        
        if cors(inode) < max(cors) - 0.2
%             doScan(GEOM.ADJsurf(inode,:)> 0) = 0;
        end
    end
end
%% select the best sink
A=[(1:length(cors))' cors rds]; 
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
    sinks = A(1);
    SCAN = adaptScanDist(GEOM,SCAN,A(1),regionR,fixed);   
else
    bestdep = initdep;
end
   
%%
function dep = adaptdep(GEOM,SCAN,initdep,fixed)

% from=getfrom(GEOM.ADJ,initdep);
ADJSLOW=GEOM.ADJ2W/(SCAN.reduceVeloFact*SCAN.MAX_MYO_VS);
done=zeros(size(initdep));
done(fixed)=1;
nextDones=zeros(size(initdep));
dep=initdep;
while sum(done) < length(done)
    for i=1:length(initdep)
        if ~done(i) 
            neigh = find(ADJSLOW(:,i) > 0 & done);
            if any(done(neigh)==1)
                dep(i) = min(dep(neigh) + ADJSLOW(neigh,i));
                nextDones(i)=1;
            end
        end
    end
    done = done + nextDones;
end
%%
function SCAN = adaptScanDist(GEOM,SCAN,sink,regionR,fixed)

% ADJSLOW= GEOM.ADJ2W/(SCAN.MAX_MYO_VS * (1 - SCAN.reduceVeloFact));
SCAN.DIST2W(~fixed,:) = SCAN.DIST2W(~fixed,:) /(SCAN.reduceVeloFact*SCAN.MAX_MYO_VS);
SCAN.DIST2W(:,~fixed) = SCAN.DIST2W(:,~fixed) /(SCAN.reduceVeloFact*SCAN.MAX_MYO_VS);

% done=zeros(size(initdep));
% done(fixed)=1;
% nextDones=zeros(size(initdep));
% dep=initdep;
% while sum(done) < length(done)
%     for i=1:length(initdep)
%         if ~done(i) 
%             neigh = find(ADJSLOW(:,i) > 0 & done);
%             if any(done(neigh)==1)
%                 dep(i) = min(dep(neigh) + ADJSLOW(neigh,i));
%                 nextDones(i)=1;
%             end
%         end
%     end
%     done = done + nextDones;
% end


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

% form tranfer matrix T
% T=zeros(nn,n2);
% T(nodes2,:)=eye(n2);
% T(nodes1,:)=INT;
interDep = dep;
interDep(nodes1) = INT * dep(nodes2);


%%
%         initT = initdeprev(inode) - 20;
%         fociD = (DIST2W(foci,inode));
%         dt = 1.5 * initdep(inode) - initdep(foci);
%         vs = fociD / dt;
%         sinkD = min(DIST2W(sinks,inode));
% %         initT = initdeprev(inode) * 0.5;
%         initT = initdeprev(inode) - 20;
%         
% %         vs = 1.5 * (fociD - initT)/SCAN.qrsduration;
%         dep = DIST2W(:,inode) * vs + initT;
%         dep = min(DIST2W(:,inode)*vs + initT,initdeprev);
%         dep= dep.* initdep./(SCAN.qrsduration-dep);
%         dep = min(dep,initdeprev);
% %         dep = min(DIST2W(:,inode) + initdeprev(inode) * 0.5,initdeprev);
% %         % compute vs such that the total activation time remains constant
% %         vs  = max(dep) / SCAN.qrsduration;
% %         dep= (DIST2W(:,inode) + initdeprev(inode) * 0.5) /vs;
% %         dep = min(dep,initdeprev);
%         dep = initdep;
%         inodeFoci = find( GEOM.DIST(inode,foci) == min(GEOM.DIST(inode,foci)) );
%         inodeSink = find( GEOM.DIST(inode,sinks) == min(GEOM.DIST(inode,sinks)) );
%         fact0 = (SCAN.qrsduration - initdep(inode))/ (SCAN.qrsduration + initdep(inode));
%         deltaDepMax = SCAN.qrsduration - initdep(inode);
% %         fact0 = (initdep(sinks(inodeSink)) - initdep(inode))/ (sinks(inodeSink) + initdep(inode));
%         for j = 1:length(initdep)
%             ijFoci = find( GEOM.DIST(j,foci) == min(GEOM.DIST(j,foci)) );
%             ijSink = find( GEOM.DIST(j,sinks) == min(GEOM.DIST(j,sinks)) );
%             if ijFoci == inodeFoci && ...
%                ijSink == inodeSink && ...
%                initdep(j) > initdep(inode) &&...
%                GEOM.DIST(j,inode) < 40 
% %                 fact0 = (initdep(sinks(ijSink)) - initdep(inode))/ (initdep(sinks(ijSink)) + initdep(inode));
% %                 fact =  fact0-(initdep(sinks(ijSink)) - initdep(j))/ (initdep(sinks(ijSink)) + initdep(j));
%                 if sinks(ijSink) ~= inode
%                     fact = GEOM.DIST(j,sinks(ijSink)) / GEOM.DIST(inode,sinks(ijSink));
%                 else
%                     fact = 1;
%                 end
%                 deltaDep = min(deltaDepMax, initdep(sinks(ijSink)));
% %                 disp(num2str(fact))
% %                 dep(j)= dep(j) +fact*dep(j);
% %                 if initdep(j) > initdep(inode) && DIST2W(j,inode) < 40
%                     dep(j)= dep(inode) + fact*deltaDep;%dep(j);
% %                 end
%             end
%         end
%             if abs(dep(j) - dep(inode)) < 0.5 * SCAN.qrsduration && DIST2W(inode,j) < 80 && ~any(foci==j)
%                 dep(j) = (DIST2W(inode,foci) ./ DIST2W(foci,j)) .* initdep(j) + ...
%                          (DIST2W(inode,foci) ./ DIST2W(inode,j)) .* DIST2W(j,inode) * vs;
%             end
%         end


        % reverse the order to obtain the dep with extra sink
%         dep = SCAN.qrsduration - dep;
