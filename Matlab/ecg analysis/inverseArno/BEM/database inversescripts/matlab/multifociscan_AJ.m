function [bestfoci,bestdep,outp,corsinit,rdsinit] = multifociscan_AJ(varargin)

% use all foci found each round
% date:05102013
% identification of one or more focal points of depolarization,
% Peter van Dam
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones
%
% copy multifociscan_publicationPigsPO.m [20150108]
% version 1: 8-jan-2015 AJ:
% (GEOM,'clusters',S.clusters,'usecor',1,'issinus',S.issinus,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'showplots',showplots,'blmode',S.blmode);

%% initialization: check/set input
GEOM                    = varargin{1};
SCAN.usecor             = 1;
SCAN.sinkscan           = 0;
SCAN.initialVelocity    = 0;
SCAN.maxVelocity        = 2.5;
SCAN.isSinus            = 0;
SCAN.clusters           = 1;
SCAN.usetime            = 60;
SCAN.regionR            = 30;
SCAN.lpass              = 5;            % Ook zonder hakkel filteren
SCAN.scanmode           = 1;
SCAN.useleads           = [];
SCAN.blmode             = 1;
useFocus                = -1;
showplots               = true;
useAntiHakkel           = true;

pp = 2;
if length(varargin) < 1
    error('This routine needs at least one parameter');
else
    while pp <= nargin
        if ischar(varargin{pp})
            key = lower(varargin{pp});
            switch key
                case 'usecor'
                    SCAN.usecor             = varargin{pp+1};        pp = pp+2;
                case 'sinkscan'
                    SCAN.sinkscan           = varargin{pp+1};        pp = pp+2;
                case 'initialvelocity'
                    SCAN.initialVelocity    = varargin{pp+1};        pp = pp+2;
                case 'maxvelocity'
                    SCAN.maxVelocity        = varargin{pp+1};        pp = pp+2;
                case 'issinus'
                    SCAN.isSinus            = varargin{pp+1};        pp = pp+2;
                case 'clusters'
                    SCAN.clusters           = max(varargin{pp+1},1); pp = pp+2;
                case 'usetime'
                    SCAN.usetime            = varargin{pp+1};        pp = pp+2;
                case 'useleads'
                    SCAN.useleads           = varargin{pp+1};        pp = pp+2;
                case 'focus'
                    useFocus                = varargin{pp+1};        pp = pp+2;
                case 'showplots'
                    showplots               = varargin{pp+1};        pp = pp+2;
                case 'lpass'
                    SCAN.lpass              = varargin{pp+1};        pp = pp+2;
                case 'blmode'
                    SCAN.blmode             = varargin{pp+1};        pp = pp+2;
                case 'scanmode'
                    SCAN.scanmode           = varargin{pp+1};        pp = pp+2;
                otherwise
                    error('unknown parameter');
            end
        end
    end
end

SCAN.VER                = GEOM.VER;                             % get vertices ventricles
SCAN.ITRI               = GEOM.ITRI;                            % get triangles ventricles
SCAN.TVER               = GEOM.TVER;                            % get vertices thorax
SCAN.TITRI              = GEOM.TITRI;                           % get triangles thorax
SCAN.qrsduration        = GEOM.SPECS.qrsduration;               % total duration QRS-complex
SCAN.DELTA_SHIFT        = 0.35;                                 % 0.35
SCAN.MAX_MYO_SHIFT      = 0.8;                                  % 0.8
SCAN.REGOP              = surflapl(GEOM.VER/1000,GEOM.ITRI,1);  % keep track of the rescaling factor 1000 --> check if meters or milimeters??
SCAN.prolongDistFact    = 1.25; 
SCAN.AMA                = GEOM.AMA;                                                 % copy A-matrix from GEOM
SCAN.rep                = 300*ones(size(GEOM.VER,1),1) + GEOM.SPECS.time_apexT;     % only needed if scanmode ~=1: peak-value of T-wave + 300 ms
SCAN.ADJ                = GEOM.ADJ;                                                 % copy neighbour matrix
SCAN.ADJsurf            = GEOM.ADJsurf;                                             % copy neighbour matrix over surface
SCAN.anisotropyRatio    = GEOM.anisotropyRatio;                                     % copy anisotropy ratio
SCAN.ADJ2W              = GEOM.ADJ2W;                                               % copy neighbour matrix aniso
SCAN.DIST2W             = GEOM.DIST2W;                                              % copy distance matrix aniso
SCAN.DIST               = GEOM.DIST;                                                % copy distance matrix
SCAN.neighADJ           = graphdist(GEOM.ITRI);
[SCAN.normals, ydum]    = trinormals(GEOM.VER,GEOM.ITRI);

% If atria are studied, introduce shift:
if strfind(GEOM.type,'atria'), SCAN.MAX_MYO_SHIFT = 0.8; end

% these are used in the sinkscan: above 120 ms, 
% might have ischemica areas (dead tissue expected), 
% else only mild adapations, no dead tissue expected
if SCAN.qrsduration > 120, SCAN.prolongDistFact = 1.75; end

%% prepare
% ECG signals are only used between start QRS and end Twave: filtering and
% baseline corection. SCAN.blmode determines order.
if SCAN.blmode == 1,
    SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),SCAN.lpass));
else
    bsmt        = baselinecor(GEOM.BSM,GEOM.SPECS.time_Vstim-10,GEOM.SPECS.endtwave);
    SCAN.PSIREF = lowpassma(bsmt(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),5);
end

% norm of bl-corrected and filtered BSM signal [Frobenius norm, sqrt(sum(diag(X'*X)))]
SCAN.normphi            = norm(SCAN.PSIREF,'fro');  

% Determine part of RMS that belongs to 'early' QRS: This can be BEFORE
% defined qrs-onset
[SCAN.initialActTime,SCAN.termtime] = calcSlowActTimes(SCAN.PSIREF(:,1:GEOM.SPECS.qrsduration));

% correct SCAN.initialActTime, if SCAN.initialVelocity = 0:
if SCAN.initialVelocity == 0, SCAN.initialActTime = 0; SCAN.initialVelocity = 1; end % !!!! --> should be 1.

disp(['useTime ' num2str(SCAN.initialActTime) ' ms  terminal time ' num2str(length(SCAN.PSIREF(:,1:GEOM.SPECS.qrsduration)) - SCAN.termtime) ' ms' ])                             
 
%% Scan for foci
bestdep     = [];   % best depolarization pattern on ventricles
bestfoci    = [];   % best foci found in inital estimation
bestsinks   = [];   % ...
outp        = [];   % output file with [best correlation, best RD, max depolarization time, maximum velocity, initial activation time]
bestcor     = -1;   % initial best correlation value [no foci found yet!]
bestrd      = 10;   % initial best RD value [no foci found yet!]
startTime   = clock;
nrClust     = 0;    % initial number of foci found --> begin is zero [no foci found yet!]
sinks       = [];   % ...
dep         = [];   % depolarization pattern
MaxVelo     = 0;    % initial maxVelo [no foci found yet!]
corsinit    = [];
rdsinit     = [];

% search number of foci untill # SCAN.clusters is reached
while nrClust < SCAN.clusters, 
    
    [SCAN,focusfound,dep,foci,bestshift,maxvelo]    = fociscan(GEOM,SCAN,dep,useFocus);
    SCAN.maxVelocity                                = min(maxvelo, SCAN.maxVelocity);
    MaxVelo                                         = max(MaxVelo,maxvelo);
    corsinit                                        = [corsinit SCAN.corsinit];
    rdsinit                                         = [rdsinit SCAN.rdsinit];
    
    % When no focus is found [focusfound ~= 1], search for focus is stopped!
    if ~focusfound, break; end
    
    % Calculate modelled RMS [PSIA] for dep values determined with
    % fociscan: initial focus!
    % At this moment, AntiHakkel is NOT USED: see getSmode.m!!
    if useAntiHakkel,
        PSIA = lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1)*(1:SCAN.qrsduration),dep,SCAN.rep,GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass);
    else
        PSIA = lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1)*(1:SCAN.qrsduration),dep,SCAN.rep,GEOM.SPECS,SCAN.scanmode),SCAN.lpass);
    end

    % calculate correlation and RD values for initial focus found with fociscan
    COR = corrcoef(PSIA,SCAN.PSIREF(:,1:size(PSIA,2)));
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
    
    % the 0.02 has been increased because in these patients scar tissue is most likely present. 
    % The optimization should stop rather earler than later(0.02 = 2%)
    if  ~isempty(bestdep) && max(bestdep) < SCAN.qrsduration + 5 && (( SCAN.usecor && cor <= bestcor + 0.02) || (~SCAN.usecor && rd  >= bestrd  - 0.03 )),
        if SCAN.usecor, SCAN.usecor = 0; warning('switched to rd'); continue; else break; end
    end
    
    bestrd      = min(rd,bestrd);
    bestcor     = max(cor,bestcor);
    bestdep     = dep;
        
    outp        = [outp;[bestcor bestrd bestshift max(dep) MaxVelo SCAN.initialActTime]];
    bestfoci    = [bestfoci foci];
    bestsinks   = [bestsinks; sinks];
    
    nrClust     = nrClusters(GEOM,unique(bestfoci));
    disp(['nr foci: ' num2str(length(unique(bestfoci)) ,3) '. clusters: ' num2str(nrClust,3) ...
        '. QRS duration: ' num2str(SCAN.qrsduration,3) '. sim: ' num2str([min(dep) max(dep)],3)...
        '. corr/rd/init: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2) ...
        '. delay:  ' datestr(datenum(clock)-datenum(startTime),'HH.MM.SS')])
    
    if size(PSIA,1) == size(GEOM.LAY,1)-1
        maxAmpl = round(max(max(abs(SCAN.PSIREF(:,1:size(PSIA,2))))))/2; % determine maximum amplitude in all BSM signals [all channels, abs(signal)/2]
        
        % visualize measured BSM in blue and simulated BSM in red
        if showplots
            figure(100); 
            clf; 
            sigplot(SCAN.PSIREF(:,1:size(PSIA,2)),'',GEOM.LAY,1/maxAmpl,'b',1,0);
            hold on;
            sigplot(PSIA,'',GEOM.LAY,1/maxAmpl,'r',1,0);
        end
    end
    
    % visualize depolarization pattern on ventricles
    if showplots
        figure(101);
        showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.qrsduration]);
        drawnow;
    end
end

%% =======================================================================
function [SCAN,focusfound,bestdep,foci,bestshift,maxvelo] = fociscan(GEOM,SCAN,initdepIn,useFocus)

focusfound      = 0;                                % no focus found yet at start
foci            = [];                               % empty foci array
bestdep         = initdepIn;                        % initially this is []. See start main loop multifociscan.m
bestshift       = -1;                               % ?
cors            = -10*ones(size(GEOM.DIST,1),1);    % correlation measure PSIA-PSIREF for each vertex
rds             = 10*ones(size(GEOM.DIST,1),1);     % norm (PSIA-PSIREF) / norm(PSIREF) for each vertex
deps            = zeros(size(GEOM.DIST));           % all depolarization sequences
depsrd          = zeros(size(GEOM.DIST));           % all depolarization sequences
maxvelo         = SCAN.maxVelocity;                 % maximum velocity activation times 
SCAN.corsinit   = [];
SCAN.rdsinit    = [];

% shifts: ??
shift       = 0*ones(size(GEOM.DIST(:,1)));         % why times 0?
shiftrd     = shift;
initdep     = initdepIn;

% Is only applicable when more than one foci are searched for: [NOT WITH FIRST FOCUS!]
if  ~isempty(initdep),
    % initial shifts are DELTA_SHIFT of the previous depolarization time for purkinje nodes and MAX_MYO_SHIFT for myocardial nodes. In case
    % the previous depolarization duration is longer than the QRS complex also the purkinje nodeshift is allowed on myocardial nodes. Most
    % probably there is more than one focus in this situation. The reason for making a diffrence in myocardial nodes and purkinje nodes is that
    % under almost all circumcantces the activation wave ends on myocardial nodes. Optimizing the depolarization sequence without changing this
    % latest site of activation causes only local small minima, resulting in an alternative way of inverse timing estimation, without the
    % posibility to delay activation!
    
    mint            = min(initdep);
    shift           = (initdep - mint)*SCAN.MAX_MYO_SHIFT+mint;
%     if SCAN.isSinus, shift(GEOM.purkinjever==1) = min(0.25*GEOM.SPECS.qrsduration,(initdep(GEOM.purkinjever==1)-mint)*SCAN.DELTA_SHIFT)+mint; end
    shift(shift<1)  = 0; % prevent suboptimalization
    
end

% scan for foci -----------------------------------------------
corMode = 1;

% % for the initial focus only the first half of the QRS is used:
% if isempty(initdep) && ~SCAN.isSinus, 
%     rrms = rms(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs+GEOM.SPECS.qrsduration));
%     maxt = find(rrms==max(rrms));
% else
    maxt = SCAN.qrsduration;
% end

% repolarization times:
rep = SCAN.rep; 

% In case of stimulation: initial focus/foci is known and fixed!
if useFocus > 0, 
    scanrange   = useFocus;
    cors        = -10*ones(length(GEOM.VER),1);
    rds         = inf(length(GEOM.VER),1);
else
    scanrange   = 1:length(GEOM.VER);           % number of possible foci locations on ventricles mesh
end

%% For all nodes that are possible foci:
for inode = scanrange,
    
    % For one / first focus, first part of if-loop, for second and later
    % focus second part of if-loop.
    if isempty(initdep),
        [dep, maxvelo]  = calcDepolarization(SCAN,GEOM,maxt,inode);                                                                         % Calculate depolarization times for focus on inode and maximum dep-time [maxt]
        deps(inode,:)   = dep;                                                                                                              % File depolarization times for focus on inode.
        if SCAN.scanmode == 6, rep = initRep(GEOM,dep); end                                                                                 % Only for scanmode == 6.
        PSIA        = lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1)*(1:maxt),dep,rep,GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);    % Calculate modelled RMS  for focus on inode.
        rds(inode)  = norm(SCAN.PSIREF(:,1:maxt)-PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');                                             % Determine RD-value for measured RMS and modelled RMS with focus on inode
        cors(inode) = compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);                                                                  % Determine correlation-value for measured RMS and modelled RMS with focus on inode  
    else
        nref        = norm(SCAN.PSIREF(:,1:maxt),'fro');
 
        % make earlier than first focus:
        shiftvec    = min(initdep):10:initdep(inode)-10;
        if isempty(shiftvec)
            shiftvec = shift(inode);
        end
        rdshift     = NaN(1,length(shiftvec));
        corshift    = NaN(1,length(shiftvec));
        depshift    = NaN(length(shiftvec),size(GEOM.VER,1));
        ishift      = 1;
        for shiftl = shiftvec
            dep     = min([shiftl + SCAN.DIST2W(:,inode) / SCAN.maxVelocity, initdep],[],2); % 
%             SCAN.maxvelocity set to maxv of first focus

            if max(dep) < SCAN.qrsduration,
                dt0             = max([min(dep),0, SCAN.initialActTime]); % do not scale below 0, initialActTime (or first dep)
                depAbove        = find(dep>dt0);
                dep(depAbove)   = dt0 + (dep(depAbove) - dt0) * (SCAN.qrsduration - dt0)/(max(dep)-dt0);
                %                 dep(depAbove) = t0+(dep-t0) * ((SCAN.qrsduration-t0) / max(dep-t0)); % leave first dep alone.
                % ?? scale from shiftl???, scale from 0 or min(initdep)??
            end
            
            depshift(ishift,:) = dep;
            if SCAN.scanmode == 6
                rep = initRep(GEOM,dep);
            end
            PSIA                = lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
            rdshift(ishift)     = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/nref;
            corshift(ishift)    = compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
            ishift              = ishift+1;
        end
        [xdum,imaxcor]  = max(corshift);
        [xdum,iminrd]   = min(rdshift);
        
        if SCAN.usecor,
            cors(inode)     = corshift(imaxcor);
            rds(inode)      = rdshift(imaxcor);
            dep             = depshift(imaxcor,:);
            deps(inode,:)   = dep;
            shift(inode)    = shiftvec(imaxcor);
            depsrd(inode,:) = depshift(iminrd,:);
            shiftrd(inode)  = shiftvec(iminrd);
        else
            [xdum,imin]     = min(rdshift);
            cors(inode)     = corshift(imin);
            rds(inode)      = rdshift(imin);
            dep             = depshift(imin,:);
            deps(inode,:)   = dep;
            shift(inode)    = shiftvec(imin);
        end
    end
end

% only for sinus-beat:
%if SCAN.isSinus, cors(GEOM.purkinjever == 0) = -10; end

% file correlation and RD values
SCAN.corsinit   = cors;
SCAN.rdsinit    = rds;

%% USE COR-VALUE [SCAN.usecor == 1] or RD-VALUE [SCAN.usecor ~= 1]:
firstRds = 0;

if SCAN.usecor,
    A                   = [(1:length(cors))' cors rds];         % file all values in one matrix
    A(cors == -10,:)    = [];                                   % exclude 'original' -10 cor values --> no cor value calculated for this vertex
    
    % outlier removal
    [n,x]           = hist(A(:,2),50);                          % construct 50 bins with cor-values
    minVal          = max(x(cumsum(n)/length(A(:,2)) < 0.1));   %
    toRemove        = find(A(:,2) < minVal);                    % lower 10 percent is removed
    A(toRemove,:)   = [];
    disp(['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,2))) minVal] ,2)])
    
    if ~isempty(diff(range(A(:,2)))), % oostep1, not for single focus evaluation: adapted by AJ 9-jan-2015
        if SCAN.usecor && diff(range(A(:,2))) < 0.15 && useFocus <= 0
            warning('Switched from cor to rd');
            SCAN.usecor = 0;
            firstRds    = 1;
            shift       = shiftrd;
            deps        = depsrd;
        end
    end  
elseif ~SCAN.usecor,
    A               = [(1:length(cors))' cors rds];
    A(cors==-10,:)  = [];
    
    % outlier removal
    [n,x]           = hist(A(:,3),50);
    maxVal          = max(x(cumsum(n)/length(A(:,3)) < 1.5 ));
    toRemove        = find(A(:,3) > maxVal);
    A(toRemove,:)   = [];
    disp(['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,3))) maxVal] ,2)])    
end


if ~isempty(initdepIn) && ~firstRds,
    maxDeps     = max(deps,[],2);
    toRemove    = unique([toRemove; find(maxDeps >= max(GEOM.SPECS.qrsduration + 5, max(initdepIn) - 2  ))]);
end

toRemove        = unique([toRemove; find(cors == -10)]);
A               = [(1:length(cors))' cors rds];
A(toRemove,:)   = [];
 
%% select focus:
if isempty(A),
    stop1   = 1;
elseif ~isempty(diff(range(A(:,2)))), % oostep1, not for single focus evaluation: adapted by AJ 9-jan-2015
    if SCAN.usecor && diff(range(A(:,2))) < 0.1 && useFocus<=0,
        A       = [];
    end
elseif ~isempty(diff(range(A(:,3)))), % oostep1, not for single focus evaluation: adapted by AJ 9-jan-2015
    if ~SCAN.usecor && min(A(:,3)) < 0.6 && diff(range(A(:,3))) < 0.3 && useFocus<=0,
        A       = [];
    end
elseif ~SCAN.usecor
    A(A(:,3)>=10,:) = [];
    A               = sortrows(A,3);
else
    A(A(:,2)< 0,:)  = [];               % cor values below zero, removed!
    A               = sortrows(A,2);    % sort values
    A               = A(end:-1:1,:);
end

% If A is not empty, select vertex with highest correlation as focus:
if ~isempty(A),
    
    A_RD = A(find(A(:,2)>max(A(:,2))*0.98),:);
    A_RD = sortrows(A_RD,3);
    
    focusfound  = 1;
    
    % A           = A(1,:);           % vertex with highest correlation
    % bestdep     = deps(A(1),:)';
    % foci        = A(1);
    % bestshift   = shift(A(1));
    
    A           = A_RD(1,:);           % vertex with highest correlation
    bestdep     = deps(A_RD(1),:)';
    foci        = A_RD(1);
    bestshift   = shift(A_RD(1));
    
end

%% ========================================================================
function [bestdep, maxvelo]= calcDepolarization(SCAN,GEOM,maxt,inode)

% maxt  = maximum depolarization time --> end SCAN.qrsduration.
% inode = node in ventricles mesh 

dt0max      = 60;                                                           % limit initial activation time
dt0         = min(dt0max,SCAN.initialActTime);                              % additional check over function calcSlowActTimes --> already max 60!! 0 ms, if SCAN.initialVelocity = 0;

dep         = GEOM.DIST25(:,inode)/SCAN.initialVelocity;                    % distance [anisotropic 2.5] devided by initial velocity!   [m]/[m s^-1] = [s]. If SCAN.initialVelocity is defined as 0, it was made 1!
depw        = SCAN.DIST2W(:,inode)/SCAN.initialVelocity;                    % distance [anisotropic] devided by initial velocity!       [m]/[m s^-1] = [s]
depAbove    = find(dep>dt0);                                                % find vertices with dep times later than SCAN.initialActTime.
qrsdurcor   = (SCAN.qrsduration-dt0);                                       % total qrs-duration, minus SCAN.initialActTime.
maxv        = (max(depw)-dt0)*SCAN.initialVelocity/qrsdurcor;               % [m s^-1]: distance to location most far away devided by time QRS-complex
maxvelo     = maxv;

% only for sinus-beat:
if SCAN.isSinus, maxv = min(SCAN.maxVelocity,maxv); end

dep(depAbove)   = dt0+(depw(depAbove)-dt0)*SCAN.initialVelocity/maxv;       % for vertices with dep times later than SCAN.initialActTime: determine dep-times based on constant velocity maxv after initial part.

bestdep         = dep;
rep             = dep + 400;                                                % adjust repolarization time, based on activation time. standard = 400 ms.

PSIA            = lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1)*(1:maxt),dep,rep,GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass); % lowpassma(A-matrix * output getSmode)
cor             = compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);                                                                    % compare PSIA and SCAN.PSIREF

if cor > 0,
    for dt = 0:2:SCAN.initialActTime,
        % the propagation velocity is assumed to be homogenoeus initially. Probably
        % the numerical accuracy of the method used is not high enough to tell the
        % difference
        
        dep         = GEOM.DIST25(:,inode)/SCAN.initialVelocity+dt;                 % distance [anisotropic 2.5] devided by initial velocity + dt [which can be negative]
        depw        = SCAN.DIST2W(:,inode)/SCAN.initialVelocity+dt;                 % distance devided by initial velocity + dt [which can be negative]
        
        if dt0 == 0,depAbove = 1:size(dep,1); else depAbove = find(dep>dt0);  end       % find vertices with dep times later than SCAN.initialActTime.
%         maxv        = (max(depw)-dt0)*SCAN.initialVelocity/(SCAN.qrsduration-dt0);  % different maxv than above! max(depw) is lower with negative dt. --> higher maxv
        
        % only for sinus-beat:
        if SCAN.isSinus, maxv = min(SCAN.maxVelocity,maxv); end

        % qrssteps = ceil(SCAN.qrsduration*.1/5);
        
        for qrsd = SCAN.qrsduration, % (SCAN.qrsduration-5*qrssteps):5:(SCAN.qrsduration+5*qrssteps),
            
            maxv            = (max(depw-dt)-dt0)*SCAN.initialVelocity/(qrsd - dt0 -dt);                                                                    % adjust maximum velocity stepwise
           
            dep(depAbove)   = (depw(depAbove)-dt0-dt)*SCAN.initialVelocity/maxv + dt0 + dt;                                                                   %
            
            PSIA            = lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1)*(1:maxt),dep,rep, GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass);    %
            rdnew           = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');                                           %
            cornew          = compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);                                                                        %
            
            if cornew > cor,
                cor     = cornew;
                bestdep = dep;
                maxvelo = maxv;
            end
            
        end
    end
end

%%========================================================================

% input the qrs part of the ECG.
function [initialActTime,termActtime] = calcSlowActTimes(PSIREF)

rrms            = rms(PSIREF);      % average QRS of BSM
drrms           = diffrows(rrms);   % simple diff of rrms
initialActTime  = 1;
fms             = 60;               % assuming QRS onset has to occur within the first 60 ms [limit activation time]
iat             = 0.008;            % value determined with trial and error --> [mV/ms]
tms             = 1;                % because at least [tms = 10] ms after/of QRS is assumed to have a slow velocity               
qst             = 0.01;
maxms           = 80;               % maximum part of QRS with slow velocity at end is [maxms = 30] ms [a bit arbitrary].

% determine slope of QRS above iat
for i = initialActTime:fms,      
    if abs(drrms(i)) < iat,    
        initialActTime = i;
    else
        break;
    end
end

% determine slope of QRS is below qst:
termActtime     = length(rrms) - tms;
SCAN.termtime   = length(rrms);
for i = length(rrms)-1:-1:max(1,length(rrms)-maxms), 
    if abs(drrms(i))  < qst, termActtime = i; end
end

termActtime = min(termActtime ,length(rrms) - tms);