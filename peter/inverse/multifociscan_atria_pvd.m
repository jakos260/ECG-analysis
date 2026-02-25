function [bestfoci,bestdep,outp]=multifociscan_atria_pvd(varargin)

% use all foci found each round
% date:05102012
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
SCAN.lpass   = 10;
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

SCAN.PwaveDuration = min(125,GEOM.SPECS.endP - GEOM.SPECS.onsetP);

SCAN.DELTA_SHIFT=0.3; %0.35
SCAN.MAX_MYO_SHIFT = 0.6; % 0.8
if strfind(GEOM.type,'atria')
    SCAN.MAX_MYO_SHIFT=0.8;
end

% these are used in the sinkscan
if SCAN.PwaveDuration > 120
    SCAN.prolongDistFact = 1.75; % might have ischemica areas (dead tissue expected)
else
    SCAN.prolongDistFact = 1.25; % only mild adapations, no dead tissue expected
end



%% prepare
% ECG signals are only used between start QRS and end Twave
if ~isempty(strfind(GEOM.type,'ventricles'))
    SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),SCAN.lpass));
    [SCAN.initialActTime,SCAN.termtime] = calcSlowActTimes( SCAN.PSIREF(:,1:GEOM.SPECS.PwaveDuration) );
else
    SCAN.PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endP),SCAN.lpass));
%     SCAN.PSIREF = SCAN.PSIREF(:,1:GEOM.SPECS.endP- GEOM.SPECS.onsetP);
    SCAN.initialActTime = 0;
    SCAN.termtime = 0;
end
% SCAN.PSIREFINIT = SCAN.PSIREF;



SCAN.AMA        = GEOM.AMA;

SCAN.normphi    = norm(SCAN.PSIREF,'fro');
SCAN.rep        = 300 * ones(size(GEOM.VER,1),1) + GEOM.SPECS.time_apexT;%only needed if scanmode ~=1

SCAN.L          = surflapl(GEOM.VER/1000,GEOM.ITRI,0);

SCAN.ADJ      = GEOM.ADJ;
SCAN.ADJsurf  = GEOM.ADJsurf;


SCAN.ADJ2W      = GEOM.ADJ2W;
SCAN.DIST2W     = GEOM.DIST2W;
SCAN.neighADJ   = graphdist(GEOM.ITRI);
% SCAN.BMAT       = forwardDipole(GEOM);
% [SCAN.normals,~]= trinormals(GEOM.VER,GEOM.ITRI);



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
useQRS = 0;
delta = 5;
% the automaic detection algorithm tends to detect the Jpoint somewhat
% earlier, therefore the scan can stop earlier
if SCAN.PwaveDuration > 100 
    delta = 20;
end

while  nrClust < 1 
    [SCAN,focusfound,dep,foci,bestshift]=fociscan(GEOM,SCAN,dep,useQRS);
    if ~focusfound
        break;
    end
    PSIA =lowpassma(SCAN.AMA * getSmode(ones(length(GEOM.VER),1) * ( 1 : SCAN.PwaveDuration ),dep,SCAN.rep, GEOM.SPECS,SCAN.scanmode,GEOM),SCAN.lpass);
    COR=corrcoef(PSIA,SCAN.PSIREF(:,1:size(PSIA,2)));
    cor = COR(2,1);
    rd  = norm(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:size(PSIA,2)),'fro');
    %     if ~useQRS && cor < 0.3
    %         useQRS =1;
    %         continue;
    %     end
    %     wrd = sum(rms(SCAN.PSIREF(:,1:size(PSIA,2)) - PSIA) ./ (0.0010 + rms(SCAN.PSIREF(:,1:size(PSIA,2))))); % weighted rd
    
    
    % the 0.02 has been increased because in these patients scar tissue is
    % most likely present. The optimization should stop rather earler than
    % later(0.02 = 2%)
    if  ~isempty(bestdep) && ...
            max(bestdep) < SCAN.PwaveDuration + 5 && ...
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
        '   QRS duration /sim: ' num2str(SCAN.PwaveDuration,3) '  ' num2str(max(dep),3) '  start QRS ' num2str(min(dep),3)...
        '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ...
        '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
    if size(PSIA,1) == size(GEOM.LAY,1)-1
        maxAmpl = (max(max(abs(PSIA))))/2;
%         figure(100);clf; sigplot(SCAN.PSIREF(:,1:size(PSIA,2)),'',GEOM.LAY,1/maxAmpl,'b',1,0);
%         hold on;
%         sigplot(PSIA,'',GEOM.LAY,1/maxAmpl,'r',1,0);
        %         dd= [pwd '\forThom\'];
        %         figure(111);leadv16(GEOM.BSM(:,GEOM.SPECS.endP:size(PSIA,2)+GEOM.SPECS.endP),PSIA,'leadsys','nim','paperspeed',200,'max',[-2.5 2.5]); saveas(gcf,[dd 'leadv12_' num2str(aa) '.png']);
        %         savemat([dd 'simECG_' num2str(aa) '.mat'],PSIA);
        %         if aa==1
        %         savemat([dd 'measECG.mat'],GEOM.BSM(:,GEOM.SPECS.endP:size(PSIA,2)+GEOM.SPECS.endP));
        %         end
    end
    figure(101);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'onodes',bestsinks,'range',[0 SCAN.PwaveDuration]);drawnow;
    
    
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
    mindep = min(initdepIn);
    shift(GEOM.purkinjever==1) = initdep(GEOM.purkinjever==1) * SCAN.DELTA_SHIFT + mindep;
    shift(GEOM.purkinjever==0) = initdep(GEOM.purkinjever==0) * SCAN.MAX_MYO_SHIFT+ mindep;
    if SCAN.isSinus
        shift(GEOM.purkinjever==1) = min(0.25 * SCAN.PwaveDuration, shift(GEOM.purkinjever==1) );
    end
    shift(shift < 1) = 0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
corMode = 1;
% for the initial focus only the first half of the QRS is used
% if isempty(initdep) && ~SCAN.isSinus && ~useQRS
%     rrms= rms(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endP));
%     maxt = find(rrms==max(rrms));%  (SCAN floor(SCAN.qrsduration *0.5);
% else
    maxt = SCAN.PwaveDuration;
% end
rep = SCAN.rep;
bestInitCor= 0.1;
for inode=1:length(GEOM.VER)
    % foci should be found in the first part of the QRS complex, if shift
    % is 0 it cannot be adpated anymore
    if isempty(initdep)
        dep = calcDepolarization(SCAN, GEOM, maxt,inode,bestInitCor);
    else
        dep = min([shift(inode) + SCAN.DIST2W(:,inode) / SCAN.maxVelocity, initdep],[],2);
        if max(dep) < SCAN.PwaveDuration
            dep = dep * (SCAN.PwaveDuration / max(dep));
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
    PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
    if isempty( initdep )
        PSIA2 =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:SCAN.PwaveDuration),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
        corsAll(inode) = compCor(SCAN.PSIREF(:,1:size(PSIA2,2)),PSIA2,corMode);
    end
    
    %     corsinit(inode) = compCor(PSIA(,1:min(maxt,SCAN.initialActTime)),SCAN.PSIREFINIT(:,1:min(maxt,SCAN.initialActTime)),corMode);
    %     corsterm(inode) = compCor(PSIA(:,SCAN.termtime:maxt),SCAN.PSIREF(:,SCAN.termtime:maxt),corMode);
    %     rdsinit(inode) = norm(SCAN.PSIREF(:,1:min(maxt,SCAN.initialActTime)) - PSIA(:,1:min(maxt,SCAN.initialActTime)),'fro')/norm(SCAN.PSIREF(:,1:min(maxt,SCAN.initialActTime)),'fro');
    rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
    cors(inode) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
    bestInitCor = max(bestInitCor, cors(inode));
    %     wrds(inode) = sum(rms(SCAN.PSIREF(:,1:maxt) - PSIA) ./ (0.0010 + rms(SCAN.PSIREF(:,1:maxt)))); % weighted rd
end


if isempty(initdep)
    if corsAll(cors==max(cors)) < 0 %%max(cors) - 0.5
        cors = corsAll;
    end
    figure(3000);clf;showpatch(GEOM.VER,GEOM.ITRI, cors,'nodes',find(cors==max(cors)),'range',[-1 max(cors)]);
end
if 0 && isempty(initdep)
    
    %'KORS et al.'
    %[VR,VL VF V1-V6]
    KI_9=[  -0.5267    0.1634   -0.2867   -0.1300    0.0500   -0.0100    0.1400    0.0600    0.5400;   % Vx
        -0.8633   -0.0733    0.9266    0.0600   -0.0200   -0.0500    0.0600   -0.1700    0.1300;   % Vy
        0.3300    0.3200   -0.0200   -0.4300   -0.0600   -0.1400   -0.2000   -0.1100    0.3100;];  % Vz
    %     KI_9=  [-KI_9(3,:);  KI_9(1,:) ;  -KI_9(2,:);]; %  variant: frontal,  left sagital,  headward
    ii = find(cors> max(cors)-0.5);
    for ij=1:length(ii)
        %         dipLoc = GEOM.VER(cors==max(cors),:);
        i = ii(ij);
        dipLoc = GEOM.VER(ii(ij),:);
        dip = estimateDipole(SCAN,GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endP),SCAN.useleads,dipLoc);
%         initialvecI = find(norm3d(dip')/max(norm3d(dip')) > 0.1);
%         initialvecI = initialvecI(1);
        initialvecI = SCAN.initialActTime;
        %         disp(['initialvect time' num2str(initialvecI)])
        %         d = dipLoc + mean(dip(:,19:20)'); % this represents the initial vector
        d = dipLoc + dip(:,initialvecI)'; % this represents the initial vector
        d = d/ norm3d(d);
        D=[dipLoc; dipLoc + d *50];
        initialVector = diff(D);
        initialVector = initialVector / norm(initialVector);
        
        angleN = dot(SCAN.normals(i,:), initialVector);
        if angleN > -0.5
            cors(i) = 0;
            rds(i) = 2;
        else
            dV = GEOM.VER - ones(length(GEOM.VER),1)*dipLoc;
            dist= norm3d(dV);
            dV = dV./ (norm3d(dV)* ones(1,3));
            angles = dot(dV',(ones(length(GEOM.VER),1)*initialVector)');
            angles(dist > 40) = [];
            angles(angles< 0.6) =[];
            angles(isnan(angles))=[];
            if length(angles) / length(dV) < 0.04
                cors(i) = 0;
                rds(i) = 2;
            end
        end
    end
    dipLoc = GEOM.VER(cors==max(cors),:);
    dip = estimateDipole(SCAN,GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endP),SCAN.useleads,dipLoc);
    initialvecI = find(norm3d(dip')/max(norm3d(dip')) > 0.1);
    initialvecI = initialvecI(1);
    d = dipLoc + dip(:,initialvecI)';
    d = d/ norm3d(d);
    D=[dipLoc; dipLoc + d*50];
    
    figure(3001);clf;showpatch(GEOM.VER,GEOM.ITRI, cors,'nodes',find(cors==max(cors)),'range',[0 max(cors)]);hold on;plot3(D(:,1),D(:,2),D(:,3),'k');plot3(D(1,1),D(1,2),D(1,3),'ok','markersize',10)
    
end


if isempty(initdep) && SCAN.isSinus
    cors(GEOM.RendoVER==1) = cors(GEOM.RendoVER==1) - 0.2;
    cors(GEOM.endoVER==0) = cors(GEOM.endoVER==0) - 0.1;
    
elseif 0 && isempty(initdep)
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
                if max(dep) < SCAN.PwaveDuration
                    dep = dep * (SCAN.PwaveDuration / max(dep));
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
    %lower 10 % is removed
    toRemove = find(A(:,2) < minVal);
    A(toRemove,:)=[];
    disp( ['cor/rd ' num2str([max(cors) min(rds)]) '  std ' num2str([std(cors) std(rds)  diff(range(A(:,2))) minVal] ,2)])
    maxDeps = max(deps,[],2);
    if ~isempty(initdepIn)
        toRemove = unique([toRemove; find(maxDeps >= max([SCAN.PwaveDuration + 5 max(initdepIn)-2]) )]);
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
    toRemove = unique([toRemove; find(maxDeps >= max(SCAN.PwaveDuration + 5, max(initdepIn) - 2  ))]);
end
toRemove = unique([toRemove; find(cors==-10)]);
if length(toRemove) == length(cors)
    toRemove = find(maxDeps > min(maxDeps));
end
A=[(1:length(cors))' cors rds];
A(toRemove,:)=[];

% the 0.1 (0.05) has been increased because in these patients scar tissue is
% most likely present. The optimization should stop rather earler than
% later

% select focus
if isempty(A)
    stop1=1;
    % elseif SCAN.usecor && diff(range(A(:,2))) < 0.1
    %     A=[];
    % elseif ~SCAN.usecor && min(A(:,3)) < 0.6 && diff(range(A(:,3))) < 0.3
    %     A=[];
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

dep(depAbove)=SCAN.initialActTime + (dep(depAbove)-SCAN.initialActTime) * SCAN.initialVelocity / min(SCAN.maxVelocity,(max(dep) - SCAN.initialActTime) * SCAN.initialVelocity/ (SCAN.PwaveDuration - SCAN.initialActTime) );
rep= dep + 400;
PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(SCAN.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
rd = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
cor =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,1);
bestdep = dep;
if cor > bestInitCor
    
    for dt=2:2:SCAN.initialActTime
        %         initialActTime = SCAN.initialActTime -dt;
        dep = SCAN.DIST2W(:,inode) / SCAN.initialVelocity;
        dep = dep + dt;
        depAbove = find(dep >  SCAN.initialActTime );
        dep(depAbove) = SCAN.initialActTime + (dep(depAbove)-SCAN.initialActTime) * SCAN.initialVelocity / min(SCAN.maxVelocity,(max(dep) - SCAN.initialActTime) * SCAN.initialVelocity/ (SCAN.PwaveDuration - SCAN.initialActTime) );
        
        PSIA =lowpassma(SCAN.AMA*getSmode(ones(length(GEOM.VER),1) * (1:maxt),dep,rep, GEOM.SPECS, SCAN.scanmode,GEOM),SCAN.lpass);
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
for i = initialActTime : 150
    if abs(drrms(i)) < 0.008
        initialActTime = i;
    else
        break;
    end
end
% initialActTime = max(40,initialActTime );
termActtime = 0;
for i=length(rrms)-1:-1:length(rrms)-30
    if abs(drrms(i))  < 0.01
        termActtime = i;
    end
end
termActtime = min(termActtime ,length(rrms) - 10);
disp(['useTime ' num2str(initialActTime) ' ms  terminal time ' num2str(length(rrms) - termActtime) ' ms' ])


%%
function BMAT = forwardDipole(GEOM)


% dipLoc = mean(GEOM.VER);
% DIPS = [dipLoc 1 0 0;
%         dipLoc 0 1 0;
%         dipLoc 0 0 1;];
% nrsc = size(DIPS,1);

nvers=0;
% ntri=0;
VER = GEOM.TVER;
ITRI= GEOM.TITRI;
ns = 1;
nver = length(GEOM.TVER);
PNTSPE(1,:)=[nvers nver-nvers+1 nver];

SIGMAS = 1;

% compute B matrix
B=zeros(nver,nver);
VERS= GEOM.TVER;
ITRIS = GEOM.TITRI;
for is=1:ns,
    for i=PNTSPE(is,2):PNTSPE(is,3),
        [B(i,PNTSPE(js,2):PNTSPE(js,3)),jsing] = rowforw(VERS,ITRIS,VER(i,:));
        k=jsing + PNTSPE(js,2)-1;
        if jsing~=0 && js~=is,
            B(i,k) = B(i,k) - 1;
        end
    end
end

% compute surface containments: SURCON(i,j)= 1 for i==j, else,
%                                          = 2 if Sj contains Si; else
%                                          = 0
SURCON=eye(ns);
for is=1:ns,
    i1=PNTSPE(is,2);
    for js=1:ns,
        if js~=is,
            i2=PNTSPE(js,2);
            i3=PNTSPE(js,3);
            test=sum(B(i1,i2:i3));
            if test > 1.5
                SURCON(is,js)=2;
            end
        end
    end
end

% compute deflation factors
[SIGMAS,KAPPA,DEFL,GAMMA]=deflat(SURCON,SIGMAS);

% deflate the B matrix
C=zeros(nver,nver);
for is=1:ns,
    for i=PNTSPE(is,2):PNTSPE(is,3),
        for js=1:ns,
            for j=PNTSPE(js,2):PNTSPE(js,3),
                C(i,j)=B(i,j)*KAPPA(is,js)-DEFL(is,js)/PNTSPE(js,1);
            end
        end
    end
end

BMAT=inv(eye(nver)-C);



%%
function G=cgmadipole(OBS,DIPS)

nobs = size(OBS,1);
nsrc = size(DIPS,1);
G=zeros(nobs,nsrc);
for k=1:nsrc,
    R=OBS-ones(nobs,1)*DIPS(k,1:3);
    r=norm3d(R);
    ss=find(r>eps);
    if isempty(ss)==0,
        G(ss,k)=dots(R(ss,:),DIPS(k,4:6))./r(ss).^3;
    end
end

G = G/(4*pi);

%%

function dip = estimateDipole(SCAN,ECG,useleads,dipLoc)

% compute the infinite medium potentials in a medium of unit conductivity
% from dipoles;  position and strength:
% GG=zeros(nver,nsrc);

DIPS = [dipLoc 1 0 0;
    dipLoc 0 1 0;
    dipLoc 0 0 1;];
G = cgmadipole(SCAN.TVER,DIPS);

% scale by conductivities
%homogeneous
% for is=1:ns,
%     for i = PNTSPE(is,2):PNTSPE(is,3),
%         G(i,:)=2*G(i,:)/(SIGMAS(is,1)+SIGMAS(is,2));
%     end
% end


% PSI = inv(eye(nver)-C) * G;
PSI=SCAN.BMAT * G;


next=0;
if next==1, % reflate the PSI matrix
    for k=1:size(DIPS,1)
        sums=zeros(ns,1);
        for is=1:ns,
            if SIGMAS(is,2)~=0,
                sums(is)=sums(is)+sum(PSI(PNTSPE(is,2):PNTSPE(is,3),k));
            end
        end
        for is=1:ns,
            theta=0;
            for js=1:ns,
                if SIGMAS(js,2)~=0,
                    theta=theta+(SIGMAS(js,1)-...
                        SIGMAS(js,2))*GAMMA(is,js)*sums(js)/(2*SIGMAS(js,2)*PNTSPE(js,1));
                end
            end
            if SIGMAS(is,2)~=0,
                for i=PNTSPE(is,2):PNTSPE(is,3),
                    PSI(i,k)=PSI(i,k)+theta;
                end
            end
        end
    end
end

% PSI = forwardDipole(GEOM,dipLoc);
PHI_dips = PSI(1:length(SCAN.TVER),:);
PHI_dips = PHI_dips - ones(length(PHI_dips ),1)*mean(PHI_dips(1:3,:));

L = PHI_dips(useleads,:);

% L=zeromean(PHI_dips(4:12,:)'); % lead vectors
% L=L';


% direct pseudo inverse
dip = L \ ECG;


