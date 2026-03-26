function [bestfoci,bestdep,outp]=multifociscan(varargin)

% use all foci found each round
% date:05102012
% identification of one or more focal points of depolarization,
% Peter van Dam
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones

SCAN.DOPLOT=0;

GEOM = varargin{1};
SCAN.usecor = 1;
SCAN.initialVelocity = 0.5;
SCAN.maxVelocity = 0.9;
SCAN.isSinus = 0;
SCAN.clusters =5;
SCAN.usetime = GEOM.SPECS.qrsduration;
SCAN.regionR = 30;
SCAN.lpass   = 1;
SCAN.scanmode= 1;
SCAN.useleads=[];

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

SCAN.DELTA_SHIFT=0.35;
SCAN.MAX_MYO_SHIFT = 0.8; % 0.8
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

% SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),SCAN.lpass));

SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs + GEOM.SPECS.qrsduration),SCAN.lpass));

% SCAN.PSIREFINIT = SCAN.PSIREF;

[SCAN.initialActTime,SCAN.termtime] = calcSlowActTimes( SCAN.PSIREF(:,1:GEOM.SPECS.qrsduration) );

SCAN.AMA        = GEOM.AMA;

SCAN.normphi    = norm(SCAN.PSIREF,'fro');
SCAN.rep        = 300 * ones(size(GEOM.VER,1),1) + GEOM.SPECS.time_apexT;%only needed if scanmode ~=1

SCAN.ADJ      = GEOM.ADJ;
SCAN.ADJsurf  = GEOM.ADJsurf;


SCAN.ADJ2W      = GEOM.ADJ2W;
SCAN.DIST2W     = GEOM.DIST2W;
SCAN.neighADJ   = graphdist(GEOM.ITRI);


%% Scan for foci
dep       = [];
bestdep   = 1000;
bestfoci  = [];
bestsinks = [];
outp      = [];
bestcor  = -1;
bestrd = 10;
startTime=clock;

nrClust = 0;
sinks=[];
useQRS = 1;
delta = 5;
% the automaic detection algorithm tends to detect the Jpoint somewhat
% earlier, therefore the scan can stop earlier
if SCAN.isSinus
    delta = -5;
elseif SCAN.qrsduration > 130
    delta = 20;
end

while  (nrClust < 2 ||  max(bestdep) > SCAN.qrsduration + delta ) && nrClust < SCAN.clusters
    [SCAN,focusfound,dep,foci,bestshift]=fociscan(GEOM,SCAN,dep,useQRS);
    if ~focusfound
        break;
    end
    PSIA =lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * ( 1 : SCAN.qrsduration ),dep,SCAN.rep, GEOM.SPECS,SCAN.scanmode,GEOM.neigh),SCAN.lpass);
    COR=corrcoef(PSIA,SCAN.PSIREF(:,1:size(PSIA,2)));
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
    
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
    
    SCAN.usecor = SCAN.usecor & length(unique(bestfoci)) < 10;
    
    bestrd = min(rd,bestrd);
    bestcor = max(cor,bestcor);
    bestdep=dep;
    
    outp=[outp;[bestcor bestrd bestshift max(dep)]];
    bestfoci=[bestfoci foci];
    bestsinks=[bestsinks; sinks];
    
    nrClust=nrClusters(GEOM,unique(bestfoci));
    disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
        '   QRS duration /sim: ' num2str(SCAN.qrsduration,3) '  ' num2str(max(dep),3) '  start QRS ' num2str(min(dep),3)...
        '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ...
        '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    if SCAN.DOPLOT
        if size(PSIA,1) == size(GEOM.LAY,1)-1
            maxAmpl = round(max(max(abs(SCAN.PSIREF))))/2;
            figure(100);clf; sigplot(SCAN.PSIREF(:,1:size(PSIA,2)),'',GEOM.LAY,1/maxAmpl,'b',1,0);
            hold on;
            sigplot(PSIA,'',GEOM.LAY,1/maxAmpl,'r',1,0);
            %         dd= [pwd '\forThom\'];
            %         figure(111);leadv16(GEOM.BSM(:,GEOM.SPECS.onsetqrs:size(PSIA,2)+GEOM.SPECS.onsetqrs),PSIA,'leadsys','nim','paperspeed',200,'max',[-2.5 2.5]); saveas(gcf,[dd 'leadv12_' num2str(aa) '.png']);
            %         savemat([dd 'simECG_' num2str(aa) '.mat'],PSIA);
            %         if aa==1
            %         savemat([dd 'measECG.mat'],GEOM.BSM(:,GEOM.SPECS.onsetqrs:size(PSIA,2)+GEOM.SPECS.onsetqrs));
            %         end
        end
%     figure(101);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.qrsduration]);%drawnow;
    %     saveasci([dd 'ppd2_dep_' num2str(aa) '.txt'],bestdep);
    %     aa=aa+1;
    end
end

%% =======================================================================
function [SCAN,focusfound,bestdep,foci,bestshift] = fociscan(GEOM,SCAN,initdepIn,useQRS)

% init
focusfound  = 0;
foci=[];
bestdep     = initdepIn;
bestshift   = -1;
cors        =-10 * ones(size(GEOM.DIST,1),1);
corsAll     = cors;
rds         = 10 * ones(size(GEOM.DIST,1),1);
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
    mindep = min(initdepIn);
    shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) * SCAN.DELTA_SHIFT + mindep;
    shift(GEOM.purkinjever==0) = initdep(GEOM.purkinjever==0) * SCAN.MAX_MYO_SHIFT + mindep;
    if SCAN.isSinus
        shift(GEOM.purkinjever==1) = min(0.25 * GEOM.SPECS.qrsduration, shift(GEOM.purkinjever==1) );
    end
    shift(shift < 1) = 0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
corMode = 1;
% for the initial focus only the first half of the QRS is used
if isempty(initdep) && ~SCAN.isSinus && ~useQRS
    rrms= rms(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs+GEOM.SPECS.qrsduration));
    maxt = find(rrms==max(rrms));%  (SCAN floor(SCAN.qrsduration *0.5);
else
    maxt = SCAN.qrsduration;
end
rep = SCAN.rep;
bestInitCor= 0.1;
for inode=1:length(GEOM.VER)
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if isempty(initdep)
        dep = calcDepolarization(SCAN, GEOM, maxt,inode,bestInitCor);
    else
        dep = min([shift(inode) + SCAN.DIST2W(:,inode) / SCAN.maxVelocity, initdep],[],2);
        if max(dep) < SCAN.qrsduration
            dep = dep * (SCAN.qrsduration / max(dep));
        end
    end
    deps(inode,:) = dep;
    
    if SCAN.scanmode == 6
        rep = initRep(GEOM,dep);
    end
    PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM.neigh),SCAN.lpass);
    if isempty( initdep )
        PSIA2 =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:SCAN.qrsduration),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM.neigh),SCAN.lpass);
        corsAll(inode) = compCor(SCAN.PSIREF(:,1:size(PSIA2,2)),PSIA2,corMode);
    end
    
    rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
    cors(inode) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
    bestInitCor = max(bestInitCor, cors(inode));
end


if isempty(initdep)
    if corsAll(cors==max(cors)) < 0 %%max(cors) - 0.5
        cors = corsAll;
    end
    if SCAN.DOPLOT
        figure(3000);clf;showpatch(GEOM.VER,GEOM.ITRI, cors,'nodes',find(cors==max(cors)),'range',[-1 max(cors)]);
    end
end

if isempty(initdep) && SCAN.isSinus
    cors(GEOM.RendoVER==1) = cors(GEOM.RendoVER==1) - 0.4;
    cors(GEOM.endoVER==0) = cors(GEOM.endoVER==0) - 0.2;    
end
if SCAN.DOPLOT
    figure(3001);clf;showpatch(GEOM.VER,GEOM.ITRI, cors,'nodes',find(cors==max(cors)),'range',[max(cors)-0.3 max(cors)]);
end
if SCAN.isSinus %&& isempty(initdep)
    cors(GEOM.purkinjever == 0) =-10;
    if isfield(GEOM,'typ')
        dist = min(GEOM.DIST(GEOM.typ > 3,:));
        rds( dist < 20 ) = rds( dist < 20 ) + 0.1;
    end
    toRemove =[];
else
    firstRds = 0;
    if SCAN.usecor
        A=[(1:length(cors))' cors rds];
        A(cors==-10,:)=[];
        % outlier removal
        [n,x]=hist(A(:,2),50);
        minVal= max(x(cumsum(n)/length(A(:,2)) < 0.1));
        %lower 10 % is removed
        toRemove = find(A(:,2) < minVal);
        A(toRemove,:)=[];
        disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,2))) minVal] ,2)])
        maxDeps = max(deps,[],2);
        if ~isempty(initdepIn)
            toRemove = unique([toRemove; find(maxDeps >= max([GEOM.SPECS.qrsduration + 5 max(initdepIn)-2]) )]);
        end
        if SCAN.usecor && diff(range(A(:,2))) < 0.15 || length(toRemove) == length(cors)
            SCAN.usecor =0;
            firstRds = 1;
        end
    end
    if ~SCAN.usecor
        A=[(1:length(cors))' cors rds];
        A(cors==-10,:)=[];
        toRemove =[];
        %      % outlier removal
        %     [n,x]=hist(A(:,3),50);
        %     maxVal= max(x(cumsum(n)/length(A(:,3)) < 1.5 ));
        %     toRemove = find(A(:,3) > maxVal);
        %     A(toRemove,:)=[];
        %     disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,3))) maxVal] ,2)])
    end
    if ~isempty(initdepIn) %&& ~firstRds
        maxDeps = max(deps,[],2);
        toRemove = unique([toRemove; find(maxDeps >= max(GEOM.SPECS.qrsduration + 5, max(initdepIn) - 2  ))]);
        %     A(:,maxDeps > max(GEOM.SPECS.qrsduration, max(initdepIn) - 5)) = [];
    end
    toRemove = unique([toRemove; find(cors==-10)]);
    if ~SCAN.isSinus %&& isempty(initdep)
        if length(toRemove) == length(cors)
            toRemove = find(maxDeps > min(maxDeps));
        end
    end
end
A=[(1:length(cors))' cors rds];

A(toRemove,:)=[];

% the 0.1 (0.05) has been increased because in these patients scar tissue is
% most likely present. The optimization should stop rather earler than
% later

% select focus
if isempty(A)
    disp('no focus found');
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
function bestdep = calcDepolarization(SCAN,GEOM,maxt,inode,bestInitCor)

% the propagtion velocity is assumed to be homogenoeus initially. Probably
% the numerical accuracy of the method used is not high enough to tell the
% difference

dep = SCAN.DIST2W(:,inode) / SCAN.initialVelocity;
depAbove = find(dep >  SCAN.initialActTime );

dep(depAbove)=SCAN.initialActTime + (dep(depAbove)-SCAN.initialActTime) * SCAN.initialVelocity / min(SCAN.maxVelocity,(max(dep) - SCAN.initialActTime) * SCAN.initialVelocity/ (SCAN.qrsduration - SCAN.initialActTime) );
rep= dep + 400;
PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(SCAN.VER),1) * (0:maxt-1),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM.neigh ),SCAN.lpass);
rd = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
cor =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);
bestdep = dep;
if cor > bestInitCor
    
    for dt=2:2:SCAN.initialActTime
        %         initialActTime = SCAN.initialActTime -dt;
        dep = SCAN.DIST2W(:,inode) / SCAN.initialVelocity;
        dep = dep + dt;
        depAbove = find(dep >  SCAN.initialActTime );
        dep(depAbove) = SCAN.initialActTime + (dep(depAbove)-SCAN.initialActTime) * SCAN.initialVelocity / min(SCAN.maxVelocity,(max(dep) - SCAN.initialActTime) * SCAN.initialVelocity/ (SCAN.qrsduration - SCAN.initialActTime) );
        
        PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM.neigh),SCAN.lpass);
        rdnew = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
        cornew =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);
        
        if min(dep) ~= dt
            stop=1;
        end
        if cornew > cor
            cor = cornew;
%             disp(num2str([inode cor dt rdnew]))
            bestdep = dep;
        end
    end
    disp(['better node ' num2str([inode cor rd])])
end



%%========================================================================
% input the qrs part of the ECG.

function [initialActTime,termActtime] = calcSlowActTimes( PSIREF )

rrms =rms(PSIREF);
drrms=diffrows(rrms);
initialActTime = 1;
for i = initialActTime : min(150,length(rrms)-30)
    if abs(drrms(i)) < 0.008
        initialActTime = i;
    else
        break;
    end
end

termActtime = 0;
for i=length(rrms)-1:-1:length(rrms)-30
    if abs(drrms(i))  < 0.01
        termActtime = i;
    end
end
termActtime = min(termActtime ,length(rrms) - 10);
disp(['initial Time ' num2str(initialActTime) ' ms  terminal time ' num2str(length(rrms) - termActtime) ' ms' ])

