%% %% MAIN SCRIPT FOR INVERSE ANALYSES OF PVC AMC STUDY%%
%
% Adaptation of script CreatePigs0_2.m, version 6-jan-2015, PvO, & beminverse.m 24-jan-2015 AJ.
% Version 1.0 AJ: Adapted for PVC AMC STUDY: NICE. 20150225
%
% The initial  values for this study are:
% mudep                 = 1e-5
% murep                 = 1e-5
% anisotropy            = 2.
% initial anisotropy    = 2.
% initial velocity      = 0 or 0.4 [m/s].
% maximum velocity      = 2.5 [m/s]. This should be high enough
% S.clusters            = 1; --> ectopic beat should only have one!
% rd-value              = 0.25.
% S.scanmode            = 1. Do not include T-wave

%% Initialize general parameters inverse analyses %%
clearvars
clear all
close all

% Define subjects for analyses:
subjects        = [1 2];                    % can be a vector
sel_nr_beats    = [];                       % define vector with beats to analyze. If empty, all beats are included.

measwrd         = [];
measrd          = [];
[bp]            = '/Users/arnojanssen/Documents/STW/PVCs/'; % main directory NICE analyses
S.DATE          = datestr(now,'yyyymmdd');
casedir         = '';
group           = [];
type            = 'ventricles';             % on which surface should the inverse calculations be performed {ventricles or atria}
saveresult      = 1;                        % [1] save results of initial estimate and optimization
saveCase        = 1;                        % [2] rewrite ecgspecs, [1] read ecgspecs when available
S.skipexisting  = 0;                        % [1] skip if result files already present in target folder (recover from crash)
S.doinit        = true;                     % if not true read initial estimate from file
S.doinv         = true;                     % perform inverse after intial estimate
S.scanmode      = 1;                        % use depolarization and repolarization to fit signal!

if ~S.doinit, S.skipexisting = false; end

S.initmode      = 0;                        % methods for intial estimate: [0] use normal fociscan, [1] test with fixed initial estimate focus at stimulation site, [n] (n = 'number') test nearest n points at same surface, 'negative number': test with focus on opposite surface
S.issinus       = 0;                        % will change beat folder also

leadset         = 'all';                    % which leads should be used in the initial esitimate
S.refLeads      = 1;                        % [1] zeromean BSM and A-matrix, [2] apply WCT to A-matrix and BSM, [-2] 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured.
S.lpass         = 5;                        % # samples in lowpassma used to filter fw results and BSM!!!
S.clusterDist   = 30;                       % used in multifociscan
S.clusters      = 1;                        % Number of initial activation locations [foci] determined for inverse calculation
S.interpol      = 0;                        % [1] removed leads are interpolated, no interpolation for bucket
wrd             = 0;                        % [0,1]; % use weighted rd in inverse
S.useScaling    = 1;
S.blmode        = 1;                        % [1] blcorrection qrsonset-T, [2] (Vstim-10)-T. Applied in multifociscan_AJ.m
S.dostampl      = false;                    % fit AP amplitude on ST segment.

showplots       = true;                     % visualize results
plotegm         = true;                     % plot egm on catherer positions
S.sinkScan      = 0;

% Inverse parameters:
S.INITIALVELOCITY   = [0.4];              % velocity of first part QRS complex, vector with undefined length --> 0 = no initial velocity!!
S.maxvelocity       = 2.5;              % maximum propagation velocity
S.ANIS              = 2;                % [2] for NICE STUDY! Anisotropy in propagation of activation [the higher the value, the slower the propagation!], vector with undefined length. 1 = isotrope
S.INITANIS          = 2;                % [2] for NICE STUDY! Anisotropy of first part QRS complex, vector with undefined length
maxit               = 400;              % maximum number of iterations
mrd                 = 0.25;             % Target RD-value
MUDEP               = 1e-5; 
murep               = 5e-5;

%% INVERSE ANALYSES FOR ALL SUBJECTS%%
for sub = subjects, % for all defined subjects
    
    %% create/define directories
    [subject, bsmdir,layfile]   = subjectdir(sub, bp, S.issinus);       % set subject filename, bsm directory and layfile
    S.geomdir                   = [bp subject '/'];    % directory geometries
    S.dirout                    = [bp subject '/results/'];             % set output results directory
    S.diroutOrg                 = S.dirout;
    if strcmp(subject,''), return; else end
    bsmfiles                    = dir(fullfile(bsmdir, '*.selecg'));    % select all bsm files
    RESULTS                     = [];                                   % clear results structure
    if isempty(sel_nr_beats), nrbeats = 1:length(bsmfiles); else nrbeats = sel_nr_beats; end;
    
    %% COMPLETE INVERSE ANALYSES FOR 1 SUBJECT: ALL BEATS%%
    for ibeat = nrbeats,
        %% INVERSE ANALYSES FOR 1 SUBJECT, 1 BEAT: posible multiple de/repolarization times%%
        for wrd = wrd,
            %if wrd, MUDEP = 5e-5; murep = 5e-5; else MUDEP = 5e-6; murep = 5e-5; end
            
            clear meas;         % why here?
            clear measinit;     % ?
            
            bsmfile     = fullfile(bsmdir, bsmfiles(ibeat).name);   % set pathname of #ibeat from bsmfiles
            beat        = bsmfiles(ibeat).name(1:end-7);            % select #ibeat from bsmfiles and remove (.selecg)
            orgname     = beat(1:strfind(beat,'beat')-2);           % remove (_beat...) from beat
            GEOM.beat   = beat;
            
            if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file'),
                T       = load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat']);   % .MAT file from sel_beat_aj.m
                remove  = T.DATA.remove;
                clear T;
            else
                error('Export file not available');
            end

            %% Load all data from selected subject
            disp('===============================================================')
            disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
            GEOM        = invInit_aj('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',S.ANIS(1),'group',group,'basedir',S.geomdir);
            GEOM.beat   = beat;
            S.dirout    = S.diroutOrg;
            
            % select geometry directory for specific subject:
            submodel    = [subject, '_model'];
            heartpath   = fullfile(S.geomdir, group, submodel, submodel);
            
            % create directory for subject in results folder:
            if ~exist(S.dirout,'dir'), mkdir(S.dirout); end
            
            %% save settings:
            if saveresult, save(fullfile(S.dirout, [subject '_' S.DATE '_Settings.mad']),'S'); end
            
            %% detect bad leads and remove or interpolate these leads:
            useLeads    = find(remove == 0)';    % BSM leads to be used
            remLeads    = find(remove == 1);     % BSM leads to be removed
            
            % interpolate OR remove missing leads! --> intripol is method of AvO
            if S.interpol,
                disp('Warning: S.interpolating missing leads!');
                T                       = intripol(GEOM.TVER,GEOM.TITRI,useLeads);  % input geometry thorax and 'good' leads BSM: T is transfer to all nodes geometry
                Trem                    = T(remLeads,:);                            % select T for removed leads only!
                GEOM.BSM(remLeads,:)    = Trem*GEOM.BSM(useLeads,:);                % interpolate 'good' leads to removed leads
                remove(remLeads)        = 0;                                        % indicate removed lead as good leads now.
                useLeads                = find(remove == 0)';                       % update useLeads
                remLeads                = find(remove == 1);                        % update remLeads
            end
            
            % Updating layout after lead are removed:
            L                                   = GEOM.LAY(2:end,:);
            L(ismember(L(:,3),find(remove)),:)  = [];
            for i = 1:size(L,1),    L(i,3)      = L(i,3) - sum(remove (1:L(i,3))); end % BSM numbering after removing channels
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
                    qtriplot([-16.8	27.3 888.3],[],'LVLatENdoII');
                end
            end
            
            % save figures:
            % create directory to save figures:
            if exist([S.dirout '/beat_figures']), else mkdir([S.dirout '/beat_figures']); end
            saveas(figure(301), [S.dirout '/beat_figures/' GEOM.beat '_RMS.jpg']);
            saveas(figure(302), [S.dirout '/beat_figures/' GEOM.beat '_DOM-T.jpg']);
            
            close all;
            
            %% INITIAL ESTIMATION & INVERSE CALCAULATION: for multiple initial velocities, initial anisotropies, overall anisotropies and depolarization parameters (mu)
            for initialvelocity = S.INITIALVELOCITY,
                for initanis    = S.INITANIS,
                    for anis    = S.ANIS,
                        for mudep = MUDEP,
                            
                            clear meas;
                            clear measinit;

                            %% CONSTRUCT INITIAL ESTIMATION, INVERSE AND OUTPUT DIRECTORIES
                            if 1, S.dirout = fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), sprintf('%s_mode%i_rd%0.2f_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2fmudep%0.2e',S.DATE,S.scanmode,mrd,S.initmode,wrd,initialvelocity,initanis,anis,mudep)); end
                            
                            % make output directory (when it does not exist yet)
                            if ~exist(S.dirout,'dir'), mkdir(S.dirout); end
                            
                            fnameout = sprintf('%s_%s_%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f',subject,beat,type,S.initmode,wrd,initialvelocity,initanis,anis);
                            
                            % create directory to save figures:
                            if exist([S.dirout '/figures']), else mkdir([S.dirout '/figures']); end
                            
                            % figure name:
                            fnameout_fig = sprintf('%s_%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f',beat,type,S.initmode,wrd,initialvelocity,initanis,anis);
                                                        
                            % check if their is already an output file present. If not, do analyses [big if-loop]:
                            if ~S.skipexisting || length(dir(fullfile(S.dirout, [fnameout '*.mat']))) < 2,
                                
                                % check if the initial anisotropy is equal to isotropy, else load anisotropy distance matrix
                                if initanis == 1, GEOM.DIST25 = GEOM.DIST; else GEOM.DIST25 = loadmat(fullfile(S.geomdir,submodel,[submodel num2str(initanis) '_' type '.dist'])); end
                                if isempty(GEOM.DIST25), error('requested distance file not found'); end
                                
                                % check if the general anisotropy is equal to isotropy, else load anisotropy distance matrix
                                if anis == 1, GEOM.DIST2W = GEOM.DIST; else GEOM.DIST2W = loadmat(fullfile(S.geomdir,submodel,[submodel num2str(anis) '_' type '.dist'])); end
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
                                    [measinit.foci,measinit.dep,measinit.outp,corsinit,rdsinit]  = multifociscan_AJ(GEOM,'clusters',mfclusters,'usecor',mfusecor,'issinus',mfissinus, ...
                                                                                    'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'focus',mffocus, ...
                                                                                    'showplots',showplots,'blmode',S.blmode);                                 
                                    measinit.anisotropyRatio                    = GEOM.anisotropyRatio;         % save anisotropy value used for initial estimation
                                    measinit.cor                                = measinit.outp(end,1);         % save correlation value for initial estimation
                                    measinit.rd                                 = measinit.outp(end,2);         % save RD value for initial estimation
                                    measinit.rep                                = initRep(GEOM,measinit.dep);   % determine repolarization times
                                    disp(['min deptime ' num2str(min(measinit.dep))])
                                end
                                                                
                                close all;
                                
                                if size(measinit.foci,2) == 1,
                                    
                                    % visualize correlation plot:
                                    qtriplot('reset');
                                    qtriplot('horizontal 2');
                                    qtriplot('bgdcolor white'); pause(0.1);
                                    qtriplot('size 1000 800'); pause(0.1);
                                    
                                    % panel 1,1
                                    qtriplot(GEOM.VER,GEOM.ITRI,'surf');
                                    qtriplot(corsinit);
                                    qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit');
                                    qtriplot('color green');
                                    qtriplot(GEOM.VER(find(corsinit>max(corsinit)-0.02),:),[])
                                    qtriplot('color white');
                                    
                                    % panel 1,2
                                    qtriplot(GEOM.VER,GEOM.ITRI,'surf_2');
                                    qtriplot(rdsinit);
                                    qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_2');
                                    qtriplot('color green');
                                    qtriplot(GEOM.VER(find(rdsinit<min(rdsinit)+0.02),:),[],'fociinit_3')
                                    qtriplot('color white');
                                    qtriplot('panel surf_2=2 1');
                                    qtriplot('panel fociinit_2=2 1');
                                    qtriplot('panel fociinit_3=2 1');
                                    
                                    qtriplot('funcolor uniheat');
                                    qtriplot(['funscale 0 ' num2str(max(corsinit))]);
                                    qtriplot('step 0.02');
                                    qtriplot('scale 60');
                                    qtriplot('text 0.1 0.96 white dots < 2% of max correlation ~0.8'); pause(0.1);
                                    
                                    pause
                                    
                                    % save figure:
                                    qtriplot(['png ' S.dirout '/figures/init_' fnameout_fig  '.png 1000 1000']); pause(0.2);
                                    
                                    qtriplot('exit');
                                    
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
                                    
                                    % inverse:
                                    meas = inverse_aj(GEOM,measinit.dep,measinit.rep,'estimateampl',estamp,'casedir',S.dirout,...
                                        'repopt','apd','maxiter',maxit,'mudep',mudep,'murep',murep,'minrd',mrd,'mode',S.scanmode,...
                                        'weighedrd',wrd,'blmode',S.blmode,'lpass',S.lpass);
                                                                        
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
                                        
                                        AltPSIAinitnohakkel                     = lowpassma(GEOM.AMA*getSmode(T,altdep,altrep,GEOM.SPECS,S.scanmode, GEOM),S.lpass);                                       
                                        
                                        PSIA                                    = lowpassma(GEOM.AMA*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,S.scanmode,GEOM),S.lpass);
                                        initrep                                 = initRep(GEOM,initdep,0);
                                        
                                        PSIAinitnohakkel                        = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,S.scanmode, GEOM),S.lpass);
                                        
                                        PSIAinit                                = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,S.scanmode),S.lpass);

                                        % Visualize final solution inverse QRST:
                                        figure('Name','Final QRST');
                                        clf;
                                        sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), '' ,GEOM.LAY,1.5,'b',1,0);
                                        hold on
                                        sigplot(PSIA,'blue = meas, red = sim',GEOM.LAY,1.5,'r',1,0);

                                        
                                        % save figure:
                                        saveas(gcf, [S.dirout '/figures/QRST_' fnameout_fig '.jpg'])
                                        
                                        % Visualize final solution inverse QRS:
                                        figure('Name','Final QRS');
                                        clf;
                                        sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:(GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint-1)),'',GEOM.LAY,1.5,'b',1,0);
                                        hold on
                                        sigplot(PSIA(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'r',1,0);
                                        hold on
                                        sigplot(PSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'blue = meas, red = sim, black = initnoH',GEOM.LAY,1.5,'k',1,0);

                                        PSIprev = PSIA;
                                        hold on
                                        
                                        % save figure:
                                        saveas(gcf, [S.dirout '/figures/QRS_' fnameout_fig '.jpg'])
                                        
                                        [mindep,imindep]    = min(meas.depfinal);                                       % first site of depolarization
                                        idepearly           = find(meas.depfinal<=(mindep+2) & meas.depfinal>(mindep)); % sites within 2 ms of first site
                                        
                                        % determine difference between earliest endo and epicardiac depolarization times
                                        measendo    = find(meas.depfinal == min(meas.depfinal(GEOM.endoVER==1)));
                                        measepi     = find(meas.depfinal == min(meas.depfinal(GEOM.endoVER~=1)));
                                        diffendoepi = meas.depfinal(measendo) - meas.depfinal(measepi);
                                        
                                        measendo_init       = find(measinit.dep == min(measinit.dep(GEOM.endoVER==1)));
                                        measepi_init        = find(measinit.dep == min(measinit.dep(GEOM.endoVER~=1)));
                                        diffendoepi_init    = measinit.dep(measendo_init) - measinit.dep(measepi_init);
                                        
                                        %% visualize initial estimation and final solution inverse with qtriplot [DEPOLARIZATION]:
                                        aptg    = 'angle 160 20 10';
                                        sz      = 'step 5';
                                        
                                        qtriplot('reset');
                                        qtriplot('horizontal 2');
                                        qtriplot('vertical 2');
                                        qtriplot('bgdcolor white'); pause(0.1);
                                        qtriplot('size 1000 1000'); pause(0.1);
                                        
                                        % panel 1,1
                                        qtriplot(GEOM.VER,GEOM.ITRI,'init_1');
                                        qtriplot(measinit.dep);
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(measendo_init,:),[],'fociendoinit_1');
                                        qtriplot('color cyan');
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(measepi_init,:),[],'fociepiinit_1');
                                        qtriplot('color cyan');
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_1');
                                        qtriplot('color white');
                                        qtriplot(aptg)
                                        
                                        % panel 1,2
                                        qtriplot(GEOM.VER,GEOM.ITRI,'init_2');
                                        qtriplot(measinit.dep);
                                        qtriplot(GEOM.VER(measendo_init,:),[],'fociendoinit_2');
                                        qtriplot('color cyan');
                                        qtriplot(GEOM.VER(measepi_init,:),[],'fociepiinit_2');
                                        qtriplot('color cyan');
                                        qtriplot(GEOM.VER(measinit.foci,:),[],'fociinit_2');
                                        qtriplot('color white');
                                        qtriplot('panel init_2=1 2');
                                        qtriplot('panel fociendoinit_2=1 2');
                                        qtriplot('panel fociepiinit_2=1 2');
                                        qtriplot('panel fociinit_2=1 2');
                                        
                                        % panel 2,1
                                        qtriplot(GEOM.VER,GEOM.ITRI,'final_1');
                                        qtriplot(meas.depfinal);
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(measendo,:),[],'fociendo_1');
                                        qtriplot('color cyan');
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(measepi,:),[],'fociepi_1');
                                        qtriplot('color cyan');
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(imindep,:),[],'focusfinal_1');
                                        qtriplot('color green');
                                        qtriplot(aptg)
                                        qtriplot(GEOM.VER(idepearly,:),[],'focifinal_1');
                                        qtriplot('color white');
                                        qtriplot(aptg)
                                        qtriplot('panel final_1=2 1');
                                        qtriplot('panel focifinal_1=2 1');
                                        qtriplot('panel fociendo_1=2 1');
                                        qtriplot('panel fociepi_1=2 1');
                                        qtriplot('panel focusfinal_1=2 1');
                                        
                                        
                                        % panel 2,2
                                        qtriplot(GEOM.VER,GEOM.ITRI,'final');
                                        qtriplot(meas.depfinal);
                                        qtriplot(GEOM.VER(measendo,:),[],'fociendo');
                                        qtriplot('color cyan');
                                        qtriplot(GEOM.VER(measepi,:),[],'fociepi');
                                        qtriplot('color cyan');
                                        qtriplot(GEOM.VER(imindep,:),[],'focusfinal');
                                        qtriplot('color green');
                                        qtriplot(GEOM.VER(idepearly,:),[],'focifinal');
                                        qtriplot('color white');
                                        qtriplot('panel final=2 2');
                                        qtriplot('panel focifinal=2 2');
                                        qtriplot('panel fociendo=2 2');
                                        qtriplot('panel fociepi=2 2');
                                        qtriplot('panel focusfinal=2 2');
                                        
                                        % color scale
                                        qtriplot('funcolor tim');
                                        qtriplot('funscale autocol');
                                        qtriplot(sz);
                                        
                                        % show values:
                                        qtriplot(['text 0.10 0.98 RD    = ' num2str(measinit.rd) '~0.5']); pause(0.1);
                                        qtriplot(['text 0.10 0.96 COR   = ' num2str(measinit.cor) '~0.5']); pause(0.1);
                                        qtriplot(['text 0.10 0.94 FOCUS = ' num2str(measinit.foci) '~0.5']); pause(0.1);
                                        qtriplot(['text 0.25 0.98 ENDO-EPI = ' num2str(diffendoepi_init) ' ms ~0.5']); pause(0.1);
                                        
                                        qtriplot(['text 0.50 0.98 RD   = ' num2str(meas.rdfinal) '~0.5']); pause(0.1);
                                        qtriplot(['text 0.50 0.96 COR  = ' num2str(meas.corfinal) '~0.5']); pause(0.1);
                                        qtriplot(['text 0.50 0.94 ITER = ' num2str(meas.iterfinal) '~0.5']); pause(0.1);
                                        qtriplot(['text 0.65 0.98 FOCUS = ' num2str(imindep) '~0.5']); pause(0.1);
                                        
                                        if size(idepearly,1) > 5, idepearly = idepearly(1:5); end
                                        qtriplot(['text 0.65 0.96 < 2ms = ' num2str(idepearly') '~0.5']); pause(0.1);
                                        qtriplot(['text 0.65 0.94 ENDO-EPI = ' num2str(diffendoepi) ' ms ~0.5']); pause(0.1);
                                        
                                        qtriplot(['text 0.80 0.06 target RD = ' num2str(mrd) '~0.3']); pause(0.1);
                                        qtriplot(['text 0.80 0.05 mudep = ' num2str(mudep) '~0.3']); pause(0.1);
                                        qtriplot(['text 0.80 0.04 murep = ' num2str(murep) '~0.3']); pause(0.1);
                                        qtriplot(['text 0.80 0.03 INIT ANIS = ' num2str(initanis) '~0.3']); pause(0.1);
                                        qtriplot(['text 0.90 0.06 ANIS = ' num2str(anis) '~0.3']); pause(0.1);
                                        qtriplot(['text 0.90 0.05 MODE = ' num2str(S.scanmode) '~0.3']); pause(0.1);
                                        qtriplot(['text 0.90 0.04 INIT VEL = ' num2str(initialvelocity) '~0.3']); pause(0.1);
                                        qtriplot('scale 70');
                                                                                
                                        pause
                                        
                                        % save figure:
                                        qtriplot(['png ' S.dirout '/figures/focus_' fnameout_fig  '.png 1000 1000']); pause(0.2);
                                        
                                        qtriplot('exit');
                                    end
                                end
                                
                                close all;
                                
                                meas.velocity       = '';
                                meas.anis           = anis;
                                meas.mudep          = mudep;
                                meas.murep          = murep;
                                meas.specs          = GEOM.SPECS;
                                measinit.specs      = GEOM.SPECS;
                                measinit.velocity   = '';
                                measinit.anis       = anis;
                                
                                if saveresult,
                                    save(fullfile(S.dirout, [fnameout '_' S.DATE 'init.mat']),'measinit', 'GEOM')
                                    save(fullfile(S.dirout ,[fnameout '_' S.DATE '.mat']),'meas', 'GEOM')
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