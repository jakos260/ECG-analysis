function [bestfoci,bestdep,outp] = multifociscan_publicationPigsPO(varargin)

% use all foci found each round
% date:05102013
% identification of one or more focal points of depolarization,
% Peter van Dam
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones


GEOM = varargin{1};
SCAN.usecor = 1;
SCAN.sinkscan = 0;
SCAN.initialVelocity = 0.5;
SCAN.maxVelocity = 0.9;
SCAN.isSinus = 0;
SCAN.clusters =5;
SCAN.usetime = 60;
SCAN.regionR = 30;
SCAN.lpass   = 5; % Ook zonder hakkel filteren
SCAN.scanmode= 1;
SCAN.useleads=[];
SCAN.blmode=1;
useFocus = -1;
showplots=true;
useAntiHakkel=true;
% warning('Hakkelversie');
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
                    SCAN.useleads    = varargin{pp+1};pp=pp+2;
                case 'focus'
                    useFocus       = varargin{pp+1};pp=pp+2;
                case 'showplots'
                    showplots       = varargin{pp+1};pp=pp+2;
                case 'lpass'
                    SCAN.lpass       = varargin{pp+1};pp=pp+2;
                case 'blmode'
                    SCAN.blmode       = varargin{pp+1};pp=pp+2;
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
SCAN.MAX_MYO_SHIFT = 0.8; % 0.8
if strfind(GEOM.type,'atria')
    SCAN.MAX_MYO_SHIFT=0.8;
end


SCAN.REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
% these are used in the sinkscan
if SCAN.qrsduration > 120
    SCAN.prolongDistFact = 1.75; % might have ischemica areas (dead tissue expected)
else
    SCAN.prolongDistFact = 1.25; % only mild adapations, no dead tissue expected
end



%% prepare
% ECG signals are only used between start QRS and end Twave

if SCAN.blmode==1
    SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),SCAN.lpass));
else
    bsmt=baselinecor(GEOM.BSM,GEOM.SPECS.time_Vstim-10,GEOM.SPECS.endtwave);
    SCAN.PSIREF=lowpassma(bsmt(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),5);
end


% SCAN.PSIREF = GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave);
% SCAN.PSIREF = baselinecor(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave));
% SCAN.PSIREF = lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),5); % fixed filter only
% bsmt=baselinecor(GEOM.BSM(:,1:GEOM.SPECS.endtwave));
% SCAN.PSIREF =lowpassma(bsmt(:,GEOM.SPECS.onsetqrs:end),5); % baseline op P-T, filt.

% SCAN.PSIREFINIT = SCAN.PSIREF;

[SCAN.initialActTime,SCAN.termtime] = calcSlowActTimes( SCAN.PSIREF(:,1:GEOM.SPECS.qrsduration) );

if SCAN.initialVelocity==0
    SCAN.initialActTime=0; % For testing single velocity
    SCAN.initialVelocity=1;
end

SCAN.AMA        = GEOM.AMA;

SCAN.normphi    = norm(SCAN.PSIREF,'fro');
SCAN.rep        = 300 * ones(size(GEOM.VER,1),1) + GEOM.SPECS.time_apexT;%only needed if scanmode ~=1

% SCAN.L          = surflapl(GEOM.VER/1000,GEOM.ITRI,1);

SCAN.ADJ      = GEOM.ADJ;
SCAN.ADJsurf  = GEOM.ADJsurf;



SCAN.anisotropyRatio= GEOM.anisotropyRatio;
SCAN.ADJ2W      = GEOM.ADJ2W;
SCAN.DIST2W     = GEOM.DIST2W;
SCAN.DIST     = GEOM.DIST;
SCAN.neighADJ   = graphdist(GEOM.ITRI);
[SCAN.normals,~]= trinormals(GEOM.VER,GEOM.ITRI);



%% Scan for foci

bestdep   = [];
bestfoci  = [];
bestsinks = [];
outp      = [];
bestcor  = -1;
bestrd = 10;
startTime=clock;

nrClust = 0;
sinks=[];

dep       = [];

MaxVelo=0;
% while  nrClust < 1 || max(bestdep) > SCAN.qrsduration + 10 %nrClust < SCAN.clusters ||
%     while  nrClust < 1 || (max(bestdep) > SCAN.qrsduration + 10  && nrClust < SCAN.clusters) % oostep1
while  nrClust < SCAN.clusters % oostep1 20131211
    
    [SCAN,focusfound,dep,foci,bestshift,maxvelo] = fociscan(GEOM,SCAN,dep,useFocus);
    SCAN.maxVelocity=min(maxvelo, SCAN.maxVelocity);
    MaxVelo=max(MaxVelo,maxvelo);
    if ~focusfound
        break;
    end
    
    if useAntiHakkel
        PSIA =lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * ( 1 : SCAN.qrsduration ),dep,SCAN.rep, GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass);
    else
        PSIA =lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * ( 1 : SCAN.qrsduration ),dep,SCAN.rep, GEOM.SPECS,SCAN.scanmode),SCAN.lpass);
    end

 %   PSIANormal=lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * ( 1 : SCAN.qrsduration ),dep,SCAN.rep, GEOM.SPECS,SCAN.scanmode),5);

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
            warning(' switched to rd');
            continue;
        else
            break;
        end
    end
    
    bestrd = min(rd,bestrd);
    bestcor = max(cor,bestcor);
    bestdep=dep;
    
    outp=[outp;[bestcor bestrd bestshift max(dep) MaxVelo SCAN.initialActTime]];
    bestfoci=[bestfoci foci];
    bestsinks=[bestsinks; sinks];
    
    nrClust=nrClusters(GEOM,unique(bestfoci));
    disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
        '   QRS duration /sim: ' num2str(SCAN.qrsduration,3) '  ' num2str([min(dep) max(dep)],3)...
        '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2) ...
        '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    if size(PSIA,1) == size(GEOM.LAY,1)-1
        maxAmpl = round(max(max(abs(SCAN.PSIREF(:,1:size(PSIA,2))))))/2;
        if showplots
            figure(100);clf; sigplot(SCAN.PSIREF(:,1:size(PSIA,2)),'',GEOM.LAY,1/maxAmpl,'b',1,0);
            hold on;
            sigplot(PSIA,'',GEOM.LAY,1/maxAmpl,'r',1,0);
            
%                 save('Pig09RecPig.mat','PSIA');
            if 0 % for publication figure
                tmp=load('Pig09RecPig.mat','PSIA');
                sigplot(tmp.PSIA,'',GEOM.LAY,1/maxAmpl,'k',1,0);
            end
%             pause
        end
    end
    if showplots
        figure(101);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.qrsduration]);drawnow;
    end
end

%% =======================================================================
function [SCAN,focusfound,bestdep,foci,bestshift,maxvelo] = fociscan(GEOM,SCAN,initdepIn, useFocus)

% init
focusfound  = 0;
foci=[];
bestdep     = initdepIn;
bestshift   = -1;
cors        = -10 * ones(size(GEOM.DIST,1),1);
rds         = 10 * ones(size(GEOM.DIST,1),1);
deps        = zeros(size(GEOM.DIST));		% all depolarization sequences
depsrd      = zeros(size(GEOM.DIST));		% all depolarization sequences
maxvelo=SCAN.maxVelocity;

% shifts
shift= 0*ones(size(GEOM.DIST(:,1)));
shiftrd=shift;
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
    mint=min(initdep);
    %     shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) * SCAN.DELTA_SHIFT + mint;
    shift = (initdep - mint) * SCAN.MAX_MYO_SHIFT + mint;
    if SCAN.isSinus
        shift(GEOM.purkinjever==1) = min(0.25 * GEOM.SPECS.qrsduration, (initdep(GEOM.purkinjever==1) -mint) * SCAN.DELTA_SHIFT ) + mint;
    end
    shift(shift < 1) = 0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
corMode = 1;
% for the initial focus only the first half of the QRS is used
% if isempty(initdep) && ~SCAN.isSinus
%     rrms= rms(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs+GEOM.SPECS.qrsduration));
%     maxt = find(rrms==max(rrms));
% else
    maxt = SCAN.qrsduration;
% end
rep = SCAN.rep;

% oostep1
if useFocus>0
    scanrange=useFocus;
    cors=-10*ones(length(GEOM.VER),1);
    rds=inf(length(GEOM.VER),1);
else
    scanrange=1:length(GEOM.VER);
end


for inode=scanrange %1:length(GEOM.VER)
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if isempty(initdep)
        [dep, maxvelo] = calcDepolarization(SCAN, GEOM,maxt,inode);
        deps(inode,:) = dep;
        if SCAN.scanmode == 6
            rep = initRep(GEOM,dep);
        end
        PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
        rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
        cors(inode) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
        
    else
        nref=norm(SCAN.PSIREF(:,1:maxt),'fro');
 
        % make earlier than first focus??
        shiftvec=min(initdep):10:initdep(inode)-10;
        if isempty(shiftvec)
            shiftvec=shift(inode);
        end
        rdshift=NaN(1,length(shiftvec));
        corshift=NaN(1,length(shiftvec));
        depshift=NaN(length(shiftvec),size(GEOM.VER,1));
        ishift=1;
        for shiftl=shiftvec
            dep = min([shiftl + SCAN.DIST2W(:,inode) / SCAN.maxVelocity, initdep],[],2); % 
%             SCAN.maxvelocity set to maxv of first focus

            if max(dep) < SCAN.qrsduration
                dt0=max([min(dep),0, SCAN.initialActTime]); % do not scale below 0, initialActTime (or first dep)
                depAbove=find(dep>dt0);
                dep(depAbove) = dt0 + (dep(depAbove) - dt0) * (SCAN.qrsduration - dt0)/(max(dep)-dt0);
                %                 dep(depAbove) = t0+(dep-t0) * ((SCAN.qrsduration-t0) / max(dep-t0)); % leave first dep alone.
                % ?? scale from shiftl???, scale from 0 or min(initdep)??
            end
            
            depshift(ishift,:) = dep;
            if SCAN.scanmode == 6
                rep = initRep(GEOM,dep);
            end
            PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
            rdshift(ishift) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/nref;
            corshift(ishift) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
            ishift=ishift+1;
        end
        [~,imaxcor]=max(corshift);
        [~,iminrd]=min(rdshift);
        
        if SCAN.usecor
            cors(inode)=corshift(imaxcor);
            rds(inode)=rdshift(imaxcor);
            dep=depshift(imaxcor,:);
            deps(inode,:)=dep;
            shift(inode)=shiftvec(imaxcor);
            depsrd(inode,:)=depshift(iminrd,:);
            shiftrd(inode)=shiftvec(iminrd);
        else
            [~,imin]=min(rdshift);
            cors(inode)=corshift(imin);
            rds(inode)=rdshift(imin);
            dep=depshift(imin,:);
            deps(inode,:)=dep;
            shift(inode)=shiftvec(imin);
        end
%         if SCAN.scanmode == 6
%             rep = initRep(GEOM,dep);
%         end
%         
    end
end

if SCAN.isSinus
    cors(GEOM.purkinjever == 0) =-10;
end

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
    
    if SCAN.usecor && diff(range(A(:,2))) < 0.15 && useFocus<=0 % oostep1, not for single focus evaluation
        warning('Switched from cor to rd');
        SCAN.usecor =0;
        firstRds = 1;
        shift=shiftrd;
        deps=depsrd;
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
elseif SCAN.usecor && diff(range(A(:,2))) < 0.1 && useFocus<=0 % oostep1, not for single focus evaluation
    A=[];
elseif ~SCAN.usecor && min(A(:,3)) < 0.6 && diff(range(A(:,3))) < 0.3 && useFocus<=0 % oostep1, not for single focus evaluation
    A=[];
elseif ~SCAN.usecor
    A(A(:,3)>=10,:)=[];
    A=sortrows(A,3);
else
    if useFocus<=0 % oostep1, not for single focus evaluation
        A(A(:,2)< 0,:)=[];
    end
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
function [bestdep, maxvelo]= calcDepolarization(SCAN,GEOM,maxt,inode)

    % initial ansiotropy 2.5
    % the propagtion velocity is assumed to be homogenoeus initially. Probably
    % the numerical accuracy of the method used is not high enough to tell the
    % difference
    maxvelo=nan;
    dt0max=60; % limit initial activation time
    dt0=min(dt0max, SCAN.initialActTime);
    
    
%     Tstim=GEOM.SPECS.time_Vstim -  GEOM.SPECS.onsetqrs; % Ruben exp.
%     dep = GEOM.DIST25(:,inode) / SCAN.initialVelocity-Tstim;
%     depw = SCAN.DIST2W(:,inode) / SCAN.initialVelocity-Tstim;

    dep = GEOM.DIST25(:,inode) / SCAN.initialVelocity;
    depw = SCAN.DIST2W(:,inode) / SCAN.initialVelocity;

    depAbove = find(dep >  dt0 );
    maxv = (max(depw) - dt0) * SCAN.initialVelocity/ (SCAN.qrsduration - dt0);
    if SCAN.isSinus
        maxv=min(SCAN.maxVelocity, maxv);            
    end
    
%     dep(depAbove) = dt0 + (depw(depAbove) - dt0) * SCAN.initialVelocity / maxv;
    dep(depAbove) = dt0 + (depw(depAbove) - dt0) * SCAN.initialVelocity / maxv;

    bestdep = dep;
    maxvelo = maxv;

    rep= dep + 400;

    PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
    cor =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);

    if cor > 0
        for dt=-16:4:SCAN.initialActTime %-16 %2:2:SCAN.initialActTime
            % the propagtion velocity is assumed to be homogenoeus initially. Probably
            % the numerical accuracy of the method used is not high enough to tell the
            % difference
            
            dep = GEOM.DIST25(:,inode) / SCAN.initialVelocity + dt;
            depw = SCAN.DIST2W(:,inode) / SCAN.initialVelocity + dt;
            depAbove = find(dep >  dt0 );
            
            maxv = (max(depw) - dt0) * SCAN.initialVelocity/ (SCAN.qrsduration - dt0);
            if SCAN.isSinus
                maxv=min(SCAN.maxVelocity, maxv);            
            end
            %             for velo = min(SCAN.maxVelocity, maxv);% No veloscan, Until
            %             peak or whole QRS??
            qrssteps=ceil(SCAN.qrsduration*.1/5);
            %             for velo = maxv -0.3: 0.1 : maxv +0.2
            for qrsd=(SCAN.qrsduration-5*qrssteps):5:(SCAN.qrsduration+5*qrssteps)
                
                %                 dep(depAbove) = dt0 + (depw(depAbove) - dt0) * SCAN.initialVelocity / velo;
                maxv = (max(depw) - dt0) * SCAN.initialVelocity/ (qrsd - dt0);
                dep(depAbove) = dt0 + (depw(depAbove) - dt0) * SCAN.initialVelocity / maxv;
                PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
                rdnew = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
                cornew =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);
                
                if cornew > cor
                    cor = cornew;
                    bestdep = dep;
                    maxvelo=maxv;
                end
            end
        end
    end

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
% initialActTime = max(20,initialActTime ) % oostep1 TEST
termActtime=length(rrms) - 10;
SCAN.termtime = length(rrms);
for i=length(rrms)-1:-1:max(1,length(rrms)-30)
    if abs(drrms(i))  < 0.01
        termActtime = i;
    end
end
termActtime = min(termActtime ,length(rrms) - 10);
disp(['useTime ' num2str(initialActTime) ' ms  terminal time ' num2str(length(rrms) - termActtime) ' ms' ])

