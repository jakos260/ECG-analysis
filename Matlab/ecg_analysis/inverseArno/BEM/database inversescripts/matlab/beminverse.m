%% %% MAIN SCRIPT FOR INVERSE ANALYSES %%
%
% Adaptation of script CreatePigs0_2.m , version 6-jan-2015, PvO.
% Version 1.0 AJ

%% Initialize general parameters inverse analyses %%
clearvars
clear all
close all

measwrd         = [];
measrd          = [];
[bp]            = bempath;                  % main directory BEM analyses
S.DATE          = datestr(now,'yyyymmdd');
S.geomdir       = [bp 'Data/Geometries/'];  % directory geometries
casedir         = '';
S.dirout        = fullfile('./results/');
group           = [];
sub             = {9};                      % {'b01', 'b03'} or { 9 10 7 8} or ... , can be an array
type            = 'ventricles';             % on which surface should the inverse calculations be performed {ventricles or atria}
saveresult      = 1;                        % [1] save results of initial estimate and optimization
saveCase        = 1;                        % [2] rewrite ecgspecs, [1] read ecgspecs when available
S.skipexisting  = 0;                        % [1] skip if result files already present in target folder (recover from crash)
S.doinit        = true;                     % if not true read initial estimate from file
S.doinv         = true;                     % perform inverse after intial estimate
if ~S.doinit, S.skipexisting = false; end

S.initmode      = 0;                        % methods for intial estimate: [0] use normal fociscan, [1] test with fixed initial estimate focus at stimulation site, [n] (n = 'number') test nearest n points at same surface, 'negative number': test with focus on opposite surface
S.issinus       = 0;                        % will change beat folder also

leadset         = 'all';                    % which leads should be used in the initial esitimate
S.refLeads      = 1;                        % [1] zeromean BSM and A-matrix, [2] apply WCT to A-matrix and BSM, [-2] 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured.
S.lpass         = 5;                        % # samples in lowpassma used to filter fw results and BSM!!!
S.clusterDist   = 30;                       % used in multifociscan
S.clusters      = 1;                        % Number of initial activation locations [foci] determined for inverse calculation
S.interpol      = 0;                        % [1] removed leads are interpolated, no interpolation for bucket
wrd             = 0;                        % [0,1]; % use weighed rd in inverse
S.useScaling    = 1;
S.blmode        = 1;                        % [1] blcorrection qrsonset-T, [2] (Vstim-10)-T. Applied in multifociscan_AJ.m
S.dostampl      = false;                    % fit AP amplitude on ST segment.

showplots       = true;                     % visualize results
plotegm         = true;                     % plot egm on catherer positions
S.sinkScan      = 0;

S.INITIALVELOCITY   = [0.01,0.4,1.0];       % velocity of first part QRS complex, vector with undefined length
S.maxvelocity       = 2.5;                  % maximum propagation velocity
velo                = 2.0;                  % ??
S.ANIS              = 1;                    % anisotropy in propagation of activation [the higher the value, the slower the propagation!], vector with undefined length. 1 = isotrope
S.INITANIS          = [0.01,1,2.5];         % anisotropy of first part QRS complex, vector with undefined length

%% INVERSE ANALYSES FOR ALL SUBJECTS%%
for sub = sub % for all defined subjects
    
    %% create/define directories
    [subject, bsmdir,layfile]   = subjectdir(sub, bp, S.issinus);       % set subject filename, bsm directory and layfile
    S.dirout                    = [bp 'Data/results/' subject '/'];     % set output results directory
    S.diroutOrg                 = S.dirout;
    if strcmp(subject,''), return; else;end
    bsmfiles                    = dir(fullfile(bsmdir, '*.selecg'));    % select all bsm files
    RESULTS                     = [];                                   % clear results structure
    
    %% COMPLETE INVERSE ANALYSES FOR 1 SUBJECT: ALL BEATS%%
    for ibeat = 1 : length(bsmfiles),
        %% INVERSE ANALYSES FOR 1 SUBJECT, 1 BEAT: posible multiple de/repolarization times%%
        for wrd = wrd,
            if wrd, MUDEP = 5e-5; murep = 5e-5; else; MUDEP = 5e-6; murep = 5e-5; end
            
            clear meas;         % why here?
            clear measinit;     % ?
            
            bsmfile     = fullfile(bsmdir, bsmfiles(ibeat).name);   % set pathname of #ibeat from bsmfiles
            beat        = bsmfiles(ibeat).name(1:end-7);            % select #ibeat from bsmfiles and remove (.selecg)
            orgname     = beat(1:strfind(beat,'beat')-2);           % remove (_beat...) from beat
            GEOM.beat   = beat;
            
            if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file')
                T       = load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat']);   % how is this file created??
                remove  = T.DATA.remove(1:64);                                          % why only 64???
                clear T;
            else
                error('Export file not available');
            end

            %% Load all data from selected subject
            disp('===============================================================')
            disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
            GEOM = invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',S.ANIS(1),'group',group,'basedir',S.geomdir);
            
            S.dirout    = S.diroutOrg;
            heartpath   = fullfile(S.geomdir, group,  subject,  subject );
            
            % create directory for subject in results folder:
            if ~exist(S.dirout,'dir'), mkdir(S.dirout); end
            
            %% save settings:
            if saveresult, save(fullfile(S.dirout, [subject '_' S.DATE '_Settings.mad']),'S'); end
            
            %% detect bad leads and remove or interpolate these leads:
            useLeads    = find(remove==0)';    % BSM leads to be used
            remLeads    = find(remove==1);     % BSM leads to be removed
            
            % interpolate OR remove missing leads!
            if S.interpol,
                disp('Warning: S.interpolating missing leads!');
                T                       = intripol(GEOM.TVER,GEOM.TITRI,useLeads);  % input geometry thorax and 'good' leads BSM: T is transfer to all nodes geometry
                Trem                    = T(remLeads,:);                            % select T for removed leads only!
                GEOM.BSM(remLeads,:)    = Trem*GEOM.BSM(useLeads,:);                % interpolate 'good' leads to removed leads
                remove(remLeads)        = 0;                                        % indicate removed lead as good leads now.
                useLeads                = find(remove==0)';                         % update useLeads
                remLeads                = find(remove==1);                          % update remLeads
            end
            
            % Updating layout after lead are removed:
            L                                   = GEOM.LAY(2:end,:);
            L(ismember(L(:,3),find(remove)),:)  = [];
            for i = 1:size(L,1), L(i,3) = L(i,3) - sum(remove (1:L(i,3))); end      % BSM numbering after removing channels
            GEOM.LAY                            = [GEOM.LAY(1,:); L];
            
            %%
            GEOM        = selectLeads(GEOM,useLeads,S.refLeads);                                                                % adjust GEOM.BSM and use zeromean or wct    
            GEOM.SPECS  = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim',1); % determine parameters AP: e.g. activation, repolarization {slopes!!! not times --> GetSmode.m}
            
            % fit AP amplitude on ST segment.
            if S.dostampl,
                rmsbsm          = rms(GEOM.BSM);                                                    % take rms of BSM
                [ap_min,ist]    = min(rmsbsm((25:100)+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint)); % take sample point minimum value rms of ST-segment
                ist             = ist + GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint;                 % take time point minimum value rms of ST-segment
                Ampl            = sqrAmpl(GEOM,ist);                                                % inverse determination of timing of Equivalent Double Layer source EDL
                
                figure; plot(rms(GEOM.BSM))
                hold on
                plot(ist,rmsbsm(ist),'r*');
                
                % show results:
                if showplots,
                    qtriplot('reset');
                    qtriplot(GEOM.VER,GEOM.ITRI,'Ventricles');
                    qtriplot(Ampl);
                    qtriplot([-16.8	27.3	888.3],[],'LVLatENdoII');
                end
            end
            
            close all;
            
            %% INITIAL ESTIMATION & INVERSE CALCAULATION: for multiple initial velocities, initial anisotropies, overall anisotropies and depolarization parameters (mu)
            for initialvelocity = S.INITIALVELOCITY,
                for initanis    = S.INITANIS,
                    for anis    = S.ANIS,
                        for mudep = MUDEP,
                            
                            clear meas;
                            clear measinit;

                            %% CONSTRUCT INITIAL ESTIMATION, INVERSE AND OUTPUT DIRECTORIES
                            if 1, S.dirout = fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), sprintf('%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2fmudep%0.2e',S.DATE,S.initmode,wrd,initialvelocity,initanis,anis,mudep)); end
                            
                            % make output directory (when it does not exist yet)
                            if ~exist(S.dirout,'dir'), mkdir(S.dirout); end
                            
                            fnameout = sprintf('%s_%s_%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f',subject,beat,type,S.initmode,wrd,initialvelocity,initanis,anis);
                            
                            % check if their is already an output file present. If not, do analyses [big if-loop]:
                            if ~S.skipexisting || length(dir(fullfile(S.dirout, [fnameout '*.mat']))) < 2,
                                
                                % check if the initial anisotropy is equal to isotropy, else load anisotropy distance matrix
                                if initanis == 1, GEOM.DIST25 = GEOM.DIST; else GEOM.DIST25 = loadmat(fullfile(S.geomdir,subject,[subject num2str(initanis) '_' type '.dist'])); end
                                if isempty(GEOM.DIST25), error('requested distance file not found'); end
                                
                                % check if the general anisotropy is equal to isotropy, else load anisotropy distance matrix
                                if anis == 1, GEOM.DIST2W = GEOM.DIST; else GEOM.DIST2W = loadmat(fullfile(S.geomdir,subject,[subject num2str(anis) '_' type '.dist'])); end
                                GEOM.anisotropyRatio = anis;
                                
                                % initial scan
                                disp(['anisotropyRatio: ' num2str(GEOM.anisotropyRatio)]);
                                surfnames = {'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % Simplified
                                
                                S.diroutnew = S.dirout;
                                ddir        = dir(fullfile(S.dirout, [fnameout '*init.mat']));
                                
                                if isempty(ddir),
                                    trydir = dir(fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), sprintf('%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f*',S.DATE,S.initmode,wrd,initialvelocity,initanis,anis)));
                                    trydir = trydir(cell2mat({trydir.isdir}));
                                    for td = trydir',
                                        ddir = dir(fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters),td.name, [fnameout '*init.mat']));
                                        if ~isempty(ddir), S.diroutnew = fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters),td.name); break; end;
                                    end
                                end
                                
                                if isempty(ddir),
                                    if wrd, wrdin1 = 'wrd1'; wrdin2 = 'wrd0'; else wrdin1 = 'wrd0'; wrdin2 = 'wrd1'; end
                                    S.diroutnew = regexprep(S.dirout,wrdin1,wrdin2);
                                    ddir        = dir(fullfile(S.diroutnew, [regexprep(fnameout,wrdin1,wrdin2) '*init.mat']));
                                end
                                
                                %% INITITIAL ESTIMATION!
                                if ~S.doinit && ~isempty(ddir),
                                    clear measinit;             % why this one also?
                                    load(fullfile(S.diroutnew,ddir(1).name),'measinit');
                                else
                                    % initial estimation based on unknown initial activation site [S.initmode = 0] or stimulation focus [S.initmode = 1]:
                                    if S.initmode == 0,
                                        mfclusters  = S.clusters;
                                        mfusecor    = 1;
                                        mfissinus   = S.issinus;
                                        mffocus     = -1;
                                    else % --> .m-files not present!!
                                        [stimpos,cathsurf]                                  = GetCatheterPosition2(sub{1},bsmfiles(ibeat).name);
                                        if S.initmode > 0 % scan initmode vertices nearest to stimulation site
                                            [pnearest, distp, trinearest, focus, distver]   = findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,S.initmode);
                                        else % scan initmode vertices nearest to vertex opposite to stimulation site
                                            oppsurf                                         = unique(surfnames(~strcmp(cathsurf,surfnames)));
                                            [pnearest, distp, trinearest, focus, distver]   = findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,-S.initmode);
                                        end
                                        mfclusters  = 1;
                                        mfusecor    = 1;
                                        mfissinus   = 0;
                                        mffocus     = focus;
                                    end
                                    
                                    % actual initial estimation
                                    [measinit.foci,measinit.dep,measinit.outp]  = multifociscan_AJ(GEOM,'clusters',mfclusters,'usecor',mfusecor,'issinus',mfissinus, ...
                                                                                    'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'focus',mffocus, ...
                                                                                    'showplots',showplots,'blmode',S.blmode);                                 
                                    measinit.anisotropyRatio                    = GEOM.anisotropyRatio;         % save anisotropy value used for initial estimation
                                    measinit.cor                                = measinit.outp(end,1);         % save correlation value for initial estimation
                                    measinit.rd                                 = measinit.outp(end,2);         % save RD value for initial estimation
                                    measinit.rep                                = initRep(GEOM,measinit.dep);   % determine repolarization times
                                    disp(['min deptime ' num2str(min(measinit.dep))])
                                end
                                
                                %% INVERSE CALCULATION!
                                if ~S.doinv % if no inverse is required, copy initial guess to inverse results
                                    
                                    meas.depfinal   = measinit.dep;
                                    meas.repfinal   = measinit.rep;
                                    meas.corfinal   = measinit.cor;
                                    meas.rdfinal    = measinit.rd;
                                    meas.log        = '';
                                    meas.iterfinal  = NaN;
                                    
                                else % inverse calculation!        
                                    estamp  = 0;
                                    maxit   = 400;
                                    mrd     = 0.15;
                                    imode   = 1;
                                    
                                    % inverse:
                                    meas = inverse_aj(GEOM,measinit.dep,measinit.rep,'estimateampl',estamp,'casedir',S.dirout,...
                                        'repopt','apd','maxiter',maxit,'mudep',mudep,'murep',murep,'minrd',mrd,'mode',imode,...
                                        'weighedrd',wrd,'blmode',S.blmode,'lpass',S.lpass);
                                    
                                    % visualize the results of the inverse calculations:
                                    if showplots,
                                        figure('Name','Final activation'); 
                                        showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',[find(meas.depfinal==min(meas.depfinal(GEOM.typ==1))), find(meas.depfinal==min(meas.depfinal(GEOM.typ==2)))]); 
                                        view(0,0);
                                        
                                        figure('Name','Final repolarization');
                                        showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal)
                                        
                                        figure('Name','Final ARI: rep-dep');
                                        showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal)
                                    end
                                    
                                    RESULTS = [RESULTS; GEOM.VER(meas.depfinal==min(meas.depfinal),:)]; % writeresults
                                    t       = 0:GEOM.SPECS.qrstduration-1;                              % construct time-vector QRST
                                    T       = ones(length(GEOM.VER),1)*t;                               % matrix vertices mesh x time vector
                                    
                                    % visualize the results of the inverse calculations:
                                    if showplots,
                                        
                                        % find init opposite:
                                        initdep                                 = measinit.dep;
                                        [xdum,verfocus]                         = min(measinit.dep);
                                        [xdum,xdum2,xdum2,altfocus]             = findnearest(GEOM.VER(verfocus,:),GEOM.VER,GEOM.ITRI,GEOM.typ,surfnames{GEOM.typ(verfocus)},-1);
                                        [xdum,altdep]                           = multifociscan_AJ(GEOM,'clusters',S.clusters,'usecor',1,'issinus',0,'focus',altfocus,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'showplots',0,'blmode',S.blmode);
                                        altdep(altdep<min(altdep(GEOM.typ==1))) = min(altdep(GEOM.typ==1)); % nothing before epi
                                        altrep                                  = initRep(GEOM,altdep,0);
                                        AltPSIAinitnohakkel                     = lowpassma(GEOM.AMA*getSmode(T,altdep,altrep,GEOM.SPECS,1, GEOM),S.lpass);                                       
                                        PSIA                                    = lowpassma(GEOM.AMA*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,1,GEOM),S.lpass);
                                        initrep                                 = initRep(GEOM,initdep,0);
                                        PSIAinitnohakkel                        = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,1, GEOM),S.lpass);
                                        PSIAinit                                = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,1),S.lpass);

                                        % Visualize final solution inverse QRST:
                                        figure('Name','Final QRST');
                                        clf;
                                        sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),'',GEOM.LAY,1.5,'b',1,0);
                                        hold on
                                        sigplot(PSIA,'',GEOM.LAY,1.5,'r',1,0);
                                        
                                        % Visualize final solution inverse QRS:
                                        figure('Name','Final QRS');
                                        clf;
                                        sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:(GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint-1)),'',GEOM.LAY,1.5,'b',1,0);
                                        hold on
                                        sigplot(PSIA(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'r',1,0);
                                        sigplot(PSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'m',1,0);
                                        hold on
                                        sigplot(AltPSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'k',1,0);
                                        
                                        PSIprev = PSIA;
                                        hold on
                                                                         
                                        [mindep,imindep]    = min(meas.depfinal);
                                        idepearly           = find(meas.depfinal<=(mindep+2) & meas.depfinal>(mindep));
                                        
                                        % visualize initial estimation and final solution inverse with qtriplot [DEPOLARIZATION]:
                                        qtriplot('reset');
                                        qtriplot('horizontal 2'); 
                                        qtriplot(GEOM.VER,GEOM.ITRI,'init');
                                        qtriplot(measinit.dep);
                                        qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit');
                                        qtriplot('color white');
                                        qtriplot(GEOM.VER,GEOM.ITRI,'final');
                                        qtriplot(meas.depfinal);
                                        qtriplot(GEOM.VER(imindep,:),[],'focusfinal');
                                        qtriplot('color green');
                                        qtriplot(GEOM.VER(idepearly,:),[],'focifinal');
                                        qtriplot('color white');
                                        qtriplot('panel final=2 1');
                                        qtriplot('panel focusfinal=2 1');
                                        qtriplot('panel focifinal=2 1');
                                        qtriplot('funcolor tim');
                                        qtriplot('funscale autocol');
                                        
                                        pause  
                                        
                                        qtriplot('exit');
                                    end    
                                end
                                
                                close all;
                                
                                meas.velocity       = velo;
                                meas.anis           = anis;
                                meas.mudep          = mudep;
                                meas.murep          = murep;
                                meas.specs          = GEOM.SPECS;
                                measinit.specs      = GEOM.SPECS;
                                measinit.velocity   = velo;
                                measinit.anis       = anis;
                                
                                if saveresult,
                                    save(fullfile(S.dirout, [fnameout '_' S.DATE 'init.mat']),'measinit')
                                    save(fullfile(S.dirout ,[fnameout '_' S.DATE '.mat']),'meas')
                                    if S.dostampl, save(fullfile(S.dirout, [GEOM.subject '_' beat '_' type '_IST' num2str(ist) '_' num2str(anis) '_' S.DATE 'stampl.mat']),'Ampl'); end
                                end
                                
                                if wrd, measwrd = meas; else measrd = meas; end
                                if ~isempty(measwrd) && ~isempty(measrd), Comparefun(measrd.depfinal,measwrd.depfinal,GEOM.VER,GEOM.ITRI); end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end