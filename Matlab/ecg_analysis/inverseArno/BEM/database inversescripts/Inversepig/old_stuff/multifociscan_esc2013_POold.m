function [bestfoci,bestdep,outp]=multifociscan_esc2013(varargin)

% use all foci found each round
% date:05102012
% identification of one or more focal points of depolarization,
% Peter van Dam
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones
% 20130906 oostep1: added usevers, only scan seection of vertices, 
% TBD not checked for inwall or sinus
% 20130910 oostep1, single velocity scale on QRS if initialvelocity==maxVelocity


GEOM = varargin{1};
SCAN.usecor = 1;
SCAN.sinkscan = 0;
SCAN.initialVelocity = 0.5;
SCAN.maxVelocity = 0.9;
SCAN.isSinus = 0;
SCAN.clusters =5;
SCAN.usetime = 60;
SCAN.regionR = 30;
SCAN.lpass   = 5;
SCAN.scanmode= 1;
SCAN.useleads=[];
SCAN.usever=1:size(GEOM.VER,1);

pp=2;
if length(varargin) < 1
    error('This routine needs at least one parameter');
else
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'usecor'
                    SCAN.usecor      = varargin{pp+1};pp=pp+2;
                case 'sinkscan'
                    SCAN.sinkscan    = varargin{pp+1};pp=pp+2;
                case 'initialvelocity'
                    SCAN.initialVelocity    = varargin{pp+1};pp=pp+2;
                case 'maxvelocity'
                    SCAN.maxVelocity = varargin{pp+1};pp=pp+2;
                case 'issinus'
                    SCAN.isSinus     = varargin{pp+1};pp=pp+2;
                case 'clusters'
                    SCAN.clusters    = max(varargin{pp+1},1);pp=pp+2;
                case 'usetime'
                    SCAN.usetime     = varargin{pp+1};pp=pp+2;
                case 'useleads'
                    SCAN.useleads     = varargin{pp+1};pp=pp+2;
                case 'usever'
                    SCAN.usever     = varargin{pp+1};pp=pp+2; % only consider these vertices
                otherwise
                    error('unknown parameter');
            end
        end
    end
end

SCAN.VER=GEOM.VER;
SCAN.ITRI=GEOM.ITRI;
SCAN.TVER=GEOM.TVER;
SCAN.TITRI=GEOM.TITRI;

SCAN.qrsduration = GEOM.SPECS.qrsduration;

SCAN.DELTA_SHIFT=0.3; %0.35
SCAN.MAX_MYO_SHIFT = 0.6; % 0.8
if strfind(GEOM.type,'atria')
    SCAN.MAX_MYO_SHIFT=0.8;
end

% these are used in the sinkscan
if SCAN.qrsduration > 120
    SCAN.prolongDistFact = 1.75; % might have ischemica areas (dead tissue expected)
else
    SCAN.prolongDistFact = 1.25; % only mild adapations, no dead tissue expected
end



%% prepare
% ECG signals are only used between start QRS and end Twave

SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),SCAN.lpass));
% SCAN.PSIREFINIT = SCAN.PSIREF;

[SCAN.initialActTime,SCAN.termtime] = calcSlowActTimes( SCAN.PSIREF(:,1:GEOM.SPECS.qrsduration) );
if SCAN.initialVelocity==SCAN.maxVelocity
    SCAN.initialActTime=0
end% oostep: force single velocity scaling

SCAN.AMA        = GEOM.AMA;

SCAN.normphi    = norm(SCAN.PSIREF,'fro');
SCAN.rep        = 300 * ones(size(GEOM.VER,1),1) + GEOM.SPECS.time_apexT;%only needed if scanmode ~=1

SCAN.L          = surflapl(GEOM.VER/1000,GEOM.ITRI,0);

SCAN.ADJ      = GEOM.ADJ;
SCAN.ADJsurf  = GEOM.ADJsurf;


SCAN.ADJ2W      = GEOM.ADJ2W;
SCAN.DIST2W     = GEOM.DIST2W;
SCAN.neighADJ   = graphdist(GEOM.ITRI);
[SCAN.normals,~]= trinormals(GEOM.VER,GEOM.ITRI);



%% Scan for foci
dep       = [];
bestdep   = [];
bestfoci  = [];
bestsinks = [];
outp      = [];
bestcor  = -1;
bestrd = 10;
startTime=clock;

nrClust = 0;
sinks=[];
aa=1;
while  nrClust < SCAN.clusters
    [SCAN,focusfound,dep,foci,bestshift]=fociscan(GEOM,SCAN,dep);
    if ~focusfound
        break;
    end
    PSIA =lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * ( 1 : SCAN.qrsduration ),dep,SCAN.rep, GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass);
    COR=corrcoef(PSIA,SCAN.PSIREF(:,1:size(PSIA,2)));
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
%     wrd = sum(rms(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA) ./ (0.0010 + rms(SCAN.PSIREF(:,1:size(PSIA,2))))); % weighted rd
    
    
    % the 0.02 has been increased because in these patients scar tissue is
    % most likely present. The optimization should stop rather earler than
    % later(0.02 = 2%)
    if  ~isempty(bestdep) && ...
            max(bestdep) < SCAN.qrsduration + 5 && ...
            (( SCAN.usecor && cor <= bestcor + 0.02) ||...
            (~SCAN.usecor && rd  >= bestrd  - 0.03 ))
        if SCAN.usecor 
            SCAN.usecor = 0;
            continue;
        else
            break;
        end
    end
    
    bestrd = min(rd,bestrd);
    bestcor = max(cor,bestcor);
    bestdep=dep;
    
    outp=[outp;[bestcor bestrd bestshift max(dep)]];
    bestfoci=[bestfoci foci];
    bestsinks=[bestsinks; sinks];
    
    nrClust=nrClusters(GEOM,unique(bestfoci));
    disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
        '   QRS duration /sim: ' num2str(SCAN.qrsduration,3) '  ' num2str(diff(range(dep)),3)...
        '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2) ...
        '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    if size(PSIA,1) == size(GEOM.LAY,1)-1
        maxAmpl = round(max(max(abs(PSIA))))/2;
        figure('Name','Initial QRS');clf; sigplot(SCAN.PSIREF(:,1:size(PSIA,2)),'',GEOM.LAY,1/maxAmpl,'b',1,0);
        hold on;
        sigplot(PSIA,'',GEOM.LAY,1/maxAmpl,'r',1,0);
%         dd= [pwd '\forThom\'];
%         figure(111);leadv16(GEOM.BSM(:,GEOM.SPECS.onsetqrs:size(PSIA,2)+GEOM.SPECS.onsetqrs),PSIA,'leadsys','nim','paperspeed',200,'max',[-2.5 2.5]); saveas(gcf,[dd 'leadv12_' num2str(aa) '.png']);
%         savemat([dd 'simECG_' num2str(aa) '.mat'],PSIA);
%         if aa==1
%         savemat([dd 'measECG.mat'],GEOM.BSM(:,GEOM.SPECS.onsetqrs:size(PSIA,2)+GEOM.SPECS.onsetqrs));
%         end
    end
    figure('Name','Intial activation');showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.qrsduration]);drawnow;
%     saveasci([dd 'ppd2_dep_' num2str(aa) '.txt'],bestdep);
%     aa=aa+1;
    
end

%% =======================================================================
function [SCAN,focusfound,bestdep,foci,bestshift] = fociscan(GEOM,SCAN,initdepIn)

% init
focusfound  = 0;
foci=[];
bestdep     = initdepIn;
bestshift   = -1;
cors        =-10 * ones(size(GEOM.DIST,1),1);
rds         = 10 * ones(size(GEOM.DIST,1),1);
wrds        = 1000 * ones(size(GEOM.DIST,1),1);
deps        = zeros(size(GEOM.DIST));		% all depolarization sequences

% shifts
shift= 0*ones(size(GEOM.DIST(:,1)));
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
    % in an alternative way of inverse timing estimation, without the
    % posibility to delay activation!
    shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) * SCAN.DELTA_SHIFT;
    shift(GEOM.purkinjever==0) = initdep(GEOM.purkinjever==0) * SCAN.MAX_MYO_SHIFT;
    if SCAN.isSinus
        shift(GEOM.purkinjever==1) = min(0.25 * GEOM.SPECS.qrsduration, shift(GEOM.purkinjever==1) );
    end
    shift(shift < 1) = 0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
corMode = 1;
% for the initial focus only the first half of the QRS is used
if isempty(initdep) && ~SCAN.isSinus
    maxt = floor(SCAN.qrsduration /2);
else
    maxt = SCAN.qrsduration;
end
rep = SCAN.rep;

% for inode=1:length(GEOM.VER)
rds=ones(length(GEOM.VER),1)*inf;
cors=ones(length(GEOM.VER),1)*-inf;
for inode=SCAN.usever
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if isempty(initdep)
        dep = calcDepolarization(SCAN, inode);
    else
        dep = min([shift(inode) + SCAN.DIST2W(:,inode) / SCAN.maxVelocity, initdep],[],2);
        if max(dep) < SCAN.qrsduration
            dep = dep * (SCAN.qrsduration / max(dep));
        end
    end
%     PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
%     cor = compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
%     if  cor > 0.35    
%         meas=inverse_pvd(GEOM,dep ,rep,'maxiter',6,'mudep',1.5e-8,'mode',1);
%         dep =meas.depfinal;
%     end
    
    deps(inode,:) = dep;
    
    if SCAN.scanmode == 6
        rep = initRep(GEOM,dep);
    end
    %     maxt = round(max(dep));
    PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
    %     corsinit(inode) = compCor(PSIA(,1:min(maxt,SCAN.initialActTime)),SCAN.PSIREFINIT(:,1:min(maxt,SCAN.initialActTime)),corMode);
    %     corsterm(inode) = compCor(PSIA(:,SCAN.termtime:maxt),SCAN.PSIREF(:,SCAN.termtime:maxt),corMode);
    %     rdsinit(inode) = norm(SCAN.PSIREF(:,1:min(maxt,SCAN.initialActTime)) - PSIA(:,1:min(maxt,SCAN.initialActTime)),'fro')/norm(SCAN.PSIREF(:,1:min(maxt,SCAN.initialActTime)),'fro');
    rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
    cors(inode) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
%     wrds(inode) = sum(rms(SCAN.PSIREF(:,1:maxt) - PSIA) ./ (0.0010 + rms(SCAN.PSIREF(:,1:maxt)))); % weighted rd
end


if isempty(initdep) && SCAN.isSinus
    cors(GEOM.RendoVER==1) = cors(GEOM.RendoVER==1) - 0.2;
    cors(GEOM.endoVER==0) = cors(GEOM.endoVER==0) - 0.1;

elseif 0&&isempty(initdep)
    %     nodes =[ find(cors==max(cors)) find(GEOM.ADJsurf(cors==max(cors),:) > 0)];
    nodes = find(cors > max(cors)-0.2);
    nodes(cors(nodes) < max(cors) - 0.4) =[];
    beginCors = max(cors);
    deps2=deps;
    for i=1:length(nodes)
        inode = nodes(i);
        nodes2=find(GEOM.ADJ(inode,:) < 30 & GEOM.ADJ(inode,:) > 0 & GEOM.ADJsurf(inode,:) == 0 );
        dep1=deps2(inode,:);
        prevCors = max(cors);
        for j=1:length(nodes2)
            node2=nodes2(j);
            deltadep = dep1(node2);
            for k=0.2: 0.2 : 1.8
                if k < 1
                    dep2 = deps2(node2,:) + (1-k) * deltadep;
                    dep = min([dep2; dep1] )';
                else
                    dep = min([deps2(node2,:); dep1 + (k - 1) * deltadep] )';
                end
                dep = dep - min(dep);
                if max(dep) < SCAN.qrsduration
                    dep = dep * (SCAN.qrsduration / max(dep));
                end
                
                if SCAN.scanmode == 6
                    rep = initRep(GEOM,dep);
                end
                PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass);
                %                 corsinit(inode) = max(corsinit(inode),compCor(PSIA(:,1:min(maxt,SCAN.initialActTime)),SCAN.PSIREFINIT(:,1:min(maxt,SCAN.initialActTime)),corMode));
                %                 corsterm(inode) = max(corsterm(inode),compCor(PSIA(:,SCAN.termtime:maxt),SCAN.PSIREF(:,SCAN.termtime:maxt),corMode));
                %                 rdsinit(inode) = min(rdsinit(inode),  norm(SCAN.PSIREF(:,1:min(maxt,SCAN.initialActTime)) - PSIA(:,1:min(maxt,SCAN.initialActTime)),'fro')/norm(SCAN.PSIREF(:,1:min(maxt,SCAN.initialActTime)),'fro'));
                cors(inode)     = max(cors(inode),compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode));
                rds(inode)     = min(rds(inode),norm(SCAN.PSIREF(:,1:maxt) - PSIA(:,1:maxt),'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro'));
                if max(cors) > beginCors + 0.05 && max(cors) > prevCors
                    disp([ num2str([k/2 beginCors max(cors) GEOM.DIST(inode,node2)  min(rds) ],2) num2str([inode node2],3)]);
                    deps(inode,:)= dep;
                    prevCors = max(cors);
                end
            end
        end
    end
end
if SCAN.isSinus %&& isempty(initdep)
    cors(GEOM.purkinjever == 0) =-10;
end

firstRds = 0;
if SCAN.usecor
    A=[(1:length(cors))' cors rds];
    A(cors==-10,:)=[];
    % outlier removal
    [n,x]=hist(A(:,2),50);
    minVal= max(x(cumsum(n)/length(A(:,2)) < 0.1));
    minVal=-2 % hack till we fix hist
    %lower 10 % is removed
    toRemove = find(A(:,2) < minVal);
    A(toRemove,:)=[];
    disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,2))) minVal] ,2)])
    
    if SCAN.usecor && diff(range(A(:,2))) < 0.15 
        SCAN.usecor =0;
        firstRds = 1;
    end
end
if ~SCAN.usecor
    A=[(1:length(cors))' cors rds];
    A(cors==-10,:)=[];

     % outlier removal
    [n,x]=hist(A(:,3),50);
    maxVal= max(x(cumsum(n)/length(A(:,3)) < 1.5 ));
    toRemove = find(A(:,3) > maxVal);
    A(toRemove,:)=[];
    disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,3))) maxVal] ,2)])
end
if ~isempty(initdepIn) && ~firstRds
    maxDeps = max(deps,[],2);
    toRemove = unique([toRemove; find(maxDeps >= max(GEOM.SPECS.qrsduration + 5, max(initdepIn) - 2  ))]);
%     A(:,maxDeps > max(GEOM.SPECS.qrsduration, max(initdepIn) - 5)) = [];
end
toRemove = unique([toRemove; find(cors==-10)]);
A=[(1:length(cors))' cors rds];
A(toRemove,:)=[];

% the 0.1 (0.05) has been increased because in these patients scar tissue is
% most likely present. The optimization should stop rather earler than
% later

% select focus
if isempty(A)
    stop1=1;
elseif SCAN.usecor && diff(range(A(:,2))) < 0.1
    A=[];
elseif ~SCAN.usecor && min(A(:,3)) < 0.6 && diff(range(A(:,3))) < 0.3
    A=[];
elseif ~SCAN.usecor
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
    bestdep     = deps(A(1),:)';
    foci        = A(1);
    bestshift   = shift(A(1));
end


%% ========================================================================
function dep = calcDepolarization(SCAN,inode)


% the propagtion velocity is assumed to be homogenoeus initially. Probably
% the numerical accuracy of the method used is not high enough to tell the
% difference
dep = SCAN.DIST2W(:,inode) / SCAN.initialVelocity;
depAbove = find(dep >  SCAN.initialActTime );

% dep(depAbove) = SCAN.initialActTime + (dep(depAbove) - SCAN.initialActTime) * SCAN.initialVelocity / min(SCAN.maxVelocity,(max(dep)- SCAN.initialActTime)/ (SCAN.qrsduration - SCAN.initialActTime) );
dep(depAbove)=SCAN.initialActTime + (dep(depAbove)-SCAN.initialActTime) * SCAN.initialVelocity / min(SCAN.maxVelocity,(max(dep) - SCAN.initialActTime) * SCAN.initialVelocity/ (SCAN.qrsduration - SCAN.initialActTime) );

% dep(depAbove) = SCAN.initialActTime + (dep(depAbove) - SCAN.initialActTime) * SCAN.initialVelocity;
% dep(depAbove) = dep(depAbove) / min(SCAN.maxVelocity,(max(dep)- SCAN.initialActTime)/ (SCAN.qrsduration - SCAN.initialActTime) );


%%========================================================================
% function interDep = iterpolatedDep(SCAN,dep,nodes2)
%
% nn=length(dep);
% n2=length(nodes2);
% n1=nn-n2;
% nodes1=zeros(1,n1);
% index=zeros(1,nn);
% index(nodes2)=1;
%
% % create submatrices
% L1=zeros(nn,nn-n2);
% L2=zeros(nn,n2);
% % form L1 and L2
% k=0;
% l=0;
%
% for j=1:nn,
%     if index(j) == 0
%         k=k+1;
%         L1(:,k) = SCAN.L(:,j);
%         nodes1(k)=j;
%     else
%         l=l+1;
%         L2(:,l) = SCAN.L(:,nodes2(l));
%     end
% end
%
% %compute the interpolating matrix
% INT=-inv(L1'*L1)*L1'*L2;
%
% interDep = dep;
% interDep(nodes1) = INT * dep(nodes2);
%
%


%%========================================================================
% input the qrs part of the ECG.

function [initialActTime,termActtime] = calcSlowActTimes( PSIREF )

rrms =rms(PSIREF);
drrms=diffrows(rrms);
initialActTime = 1;
for i = initialActTime : 60
    if abs(drrms(i)) < 0.008
        initialActTime = i;
    else
        break;
    end
end
% initialActTime = max(40,initialActTime );
SCAN.termtime = length(rrms);
for i=length(rrms)-1:-1:length(rrms)-30
    if abs(drrms(i))  < 0.01
        termActtime = i;
    end
end
termActtime = min(termActtime ,length(rrms) - 10);
disp(['useTime ' num2str(initialActTime) ' ms  terminal time ' num2str(length(rrms) - termActtime) ' ms' ])



