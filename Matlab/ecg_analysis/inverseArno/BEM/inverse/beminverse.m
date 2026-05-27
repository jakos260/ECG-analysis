%% %% MAIN SCRIPT FOR INVERSE ANALYSES %%
%
% This is a template for the main script to run your inverse calculations.
%
% The first part contains all the variables that can be changed.
%
% The general parameter values for human studies are:
% anisotropy            = 2.
% initial anisotropy    = 2.
% initial velocity      = 0.4 [m/s].
% maximum velocity      = 2.5 [m/s]. This should be high enough
% S.clusters            = 2; --> ectopic beats are determined with 2!
% rd-value              = 0.25.
% mudep                 = [5e-6 1e-5 5e-5]
% murep                 = 1e-5
% S.scanmode            = 1 [Do not include T-wave] or 4 [include repolarization]
% S.lpass               = 1;
%
% For a new study, the datapaths in the matlab script [inversepath] have to
% adjusted.
%
% Version 1: 01-apr-2015

%% General parameters to run inverse analyses %%
clearvars
clear all
close all

% Define subjects for analyses:
subject_numbers = 1;                                                                        % can be a vector
sel_nr_beats    = [];                                                                           % define vector with beats to analyze. If empty, all beats are included.

measwrd         = [];
measrd          = [];
[bp]            = inversepath(1);                                                               % main directory analyses
S.DATE          = datestr(now,'yyyymmdd');                                                      % date-stamp
casedir         = '';
group           = [];
type            = 'ventricles';                                                                 % on which surface should the inverse calculations be performed {ventricles or atria}
saveresult      = 1;                                                                            % [1] save results of initial estimate and optimization
saveCase        = 1;                                                                            % [2] rewrite ecgspecs, [1] read ecgspecs when available
S.skipexisting  = 0;                                                                            % [1] skip if result files already present in target folder (recover from crash)
S.doinit        = true;                                                                         % if not true read initial estimate from file
S.doinv         = true;                                                                         % perform inverse after intial estimate
S.scanmode      = 4;                                                                            % use depolarization and repolarization to fit signal!
S.savefigures   = 1;                                                                            % [1] create subfolder for figures and save results as .jpg & .png
S.init_cor_rd   = 0;                                                                            % [1] show and save correlation and RD plot with qtriplot of initial estimation;

if ~S.doinit, S.skipexisting = false; end

S.initmode      = 0;                                                                            % methods for intial estimate: [0] use normal fociscan, [1] test with fixed initial estimate focus at stimulation site, [n] (n = 'number') test nearest n points at same surface, 'negative number': test with focus on opposite surface
S.issinus       = 0;                                                                            % will change beat folder also

leadset         = 'all';                                                                        % which leads should be used in the initial esitimate
S.refLeads      = 1;                                                                            % [1] zeromean BSM and A-matrix, [2] apply WCT to A-matrix and BSM, [-2] 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured.
S.lpass         = 1;                                                                            % # samples in lowpassma used to filter fw results and BSM!!!
S.clusters      = 1;                                                                            % Number of initial activation locations [foci] determined for inverse calculation
S.interpol      = 0;                                                                            % [1] removed leads are interpolated, no interpolation for bucket
wrd             = 0;                                                                            % [0,1]; % use weighted rd in inverse
S.useScaling    = 1;
S.blmode        = 1;                                                                            % baseline correction [1] qrsonset-T, [2] (Vstim-10)-T. Order is also turned around with lowpassma.m! Applied in multifociscan.m
S.dostampl      = false;                                                                        % fit AP amplitude on ST segment.
S.bsm_range     = 25:100;                                                                       % samples after J-point where to find minimum in ST-segment with S.dostampl
showplots       = true;                                                                         % visualize results
plotegm         = true;                                                                         % plot egm on catherer positions
S.sinkScan      = 0;
estamp          = 0;                                                                            % [1] use optimizeAmpl.m instead of optimizeDepRep.m
sigplotscale    = 1.5;                                                                          % Scale use for sigplot.m at the end of this script.

%% Standard inverse parameters:
S.INITIALVELOCITY   = 0.4;                                                                      % velocity of first part QRS complex, vector with undefined length --> 0 = no initial velocity!!
S.maxvelocity       = 2.5;                                                                      % maximum propagation velocity
S.ANIS              = 2;                                                                        % [2] for NICE STUDY! Anisotropy in propagation of activation [the higher the value, the slower the propagation!], vector with undefined length. 1 = isotrope
S.INITANIS          = 2;                                                                        % [2] for NICE STUDY! Anisotropy of first part QRS complex, vector with undefined length
maxit               = 400;                                                                      % maximum number of iterations
mrd                 = 0.25;                                                                     % Target RD-value
MUDEP               = [5e-6 1e-5 5e-5];                                                         % Weight of regularization with Laplacian for depolarization
murep               = 1e-5;                                                                     % Weight of regularization with Laplacian for repolarization

%% INVERSE ANALYSES FOR ALL SUBJECTS%%
for sub = subject_numbers,                                                                      % for all defined subjects
    
    %% create/define directories
    [subject, bsmdir,layfile]   = subjectdir(sub, bp, S.issinus);                               % set subject filename, bsm directory and layfile
    S.geomdir                   = [bp subject '/'];                                             % directory geometries
    S.dirout                    = [bp subject '/results/'];                                     % set output results directory
    S.diroutOrg                 = S.dirout;
    
    if strcmp(subject,''), return; else end
    
    bsmfiles                    = dir(fullfile(bsmdir, '*.selecg'));                            % select all bsm files
    RESULTS                     = [];                                                           % clear results structure
    
    if isempty(sel_nr_beats), nrbeats = 1:length(bsmfiles); else nrbeats = sel_nr_beats; end;   % check which beats to analyze
    
    %% COMPLETE INVERSE ANALYSES FOR 1 SUBJECT: ALL BEATS%%
    for ibeat = nrbeats,
        
        clear meas;
        clear measinit;
        
        bsmfile     = fullfile(bsmdir, bsmfiles(ibeat).name);                                   % set pathname of #ibeat from bsmfiles
        beat        = bsmfiles(ibeat).name(1:end-7);                                            % select #ibeat from bsmfiles and remove (.selecg)
        orgname     = beat(1:strfind(beat,'beat')-2);                                           % remove (_beat...) from beat
        GEOM.beat   = beat;                                                                     % add beat file name to GEOM-structure
        
        if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file'),
            T       = load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat']);               % .MAT file from ViewandExport.m
            remove  = T.DATA.remove;
            clear T;
        else
            error('Export file not available');
        end
        
        %% Load all data from selected subject
        disp('===============================================================')
        disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
        
        GEOM = invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type, ...       % compile all start data and information into one structure
            'anisotropyRatio',S.ANIS(1),'group',group,'basedir',S.geomdir);
        
        GEOM.beat   = beat;                                                                     % add beat file name to GEOM-structure
        S.dirout    = S.diroutOrg;
        
        submodel    = [subject, '_model'];                                                      % select geometry directory for specific subject:
        heartpath   = fullfile(S.geomdir, group, submodel, submodel);
        
        if ~exist(S.dirout,'dir'), mkdir(S.dirout); end                                         % create directory for subject in results folder:
        
        if saveresult, save(fullfile(S.dirout, [subject '_' S.DATE '_Settings.mad']),'S'); end  % save settings
        
        % detect bad leads and remove or interpolate these leads:
        useLeads    = find(remove == 0)';                                                       % BSM leads to be used
        remLeads    = find(remove == 1);                                                        % BSM leads to be removed
        
        % interpolate OR remove missing leads! --> intripol is method of AvO
        if S.interpol,
            disp('Warning: S.interpolating missing leads!');
            T                       = intripol(GEOM.TVER,GEOM.TITRI,useLeads);                  % input geometry thorax and 'good' leads BSM: T is transfer to all nodes geometry
            Trem                    = T(remLeads,:);                                            % select T for removed leads only!
            GEOM.BSM(remLeads,:)    = Trem*GEOM.BSM(useLeads,:);                                % interpolate 'good' leads to removed leads
            remove(remLeads)        = 0;                                                        % indicate removed lead as good leads now.
            useLeads                = find(remove == 0)';                                       % update useLeads
            remLeads                = find(remove == 1);                                        % update remLeads
        end
        
        % Updating layout after lead(s) are removed:
        L                                   = GEOM.LAY(2:end,:);
        L(ismember(L(:,3),find(remove)),:)  = [];
        for i = 1:size(L,1),    L(i,3)      = L(i,3) - sum(remove (1:L(i,3))); end              % BSM numbering after removing channels
        GEOM.LAY                            = [GEOM.LAY(1,:); L];
        
        %% If neccessary, create new .ecgspecs & .selecg:
        GEOM        = selectLeads(GEOM,useLeads,S.refLeads);                                    % adjust GEOM.BSM and use zeromean or wct
        GEOM.SPECS  = prepareECG(GEOM.BSM, GEOM.LAY, GEOM.neigh, 'documsum', 1, 'filename', ... % determine parameters AP: e.g. activation, repolarization {slopes!!! not times --> GetSmode.m}
            bsmfile(1:end-7), 'dosave', saveCase, 'dovstim', 1);
        
        close all;                                                                              % close figures opened with prepareECG.m
        
        %% Fit AP amplitude on ST segment when S.dostampl == true:
        if S.dostampl,
            rmsbsm          = rms(GEOM.BSM);                                                                % take rms of BSM
            [ap_min,ist]    = min(rmsbsm((S.bsm_range) + GEOM.SPECS.onsetqrs + GEOM.SPECS.time_Jpoint));    % take sample point minimum value rms of ST-segment
            ist             = ist + GEOM.SPECS.onsetqrs + GEOM.SPECS.time_Jpoint + min(S.bsm_range);        % take time point minimum value rms of ST-segment
            
            % visualize result detection time point minimum value rms of ST-segment:
            figure; 
            plot(rms(GEOM.BSM))
            hold on
            plot(ist,rmsbsm(ist),'r*');
            
            % inverse determination of timing of Equivalent Double Layer source EDL:
            Ampl            = sqrAmpl(GEOM,ist);
            
            % show results:
            if showplots,
                qtriplot('reset');
                qtriplot(GEOM.VER,GEOM.ITRI,'Ventricles');
                qtriplot(Ampl);
                qtriplot([-16.8	27.3 888.3],[],'LVLatENdoII');
            end
        end
        
        %% INITIAL ESTIMATION & INVERSE CALCAULATION: for multiple initial velocities, initial anisotropies, overall anisotropies and depolarization parameters (mu)
        for initialvelocity = S.INITIALVELOCITY,
            for initanis    = S.INITANIS,
                for anis    = S.ANIS,
                    for mudep = MUDEP,
                        
                        clear meas;
                        clear measinit;
                        
                        %% CONSTRUCT INITIAL ESTIMATION, INVERSE AND OUTPUT DIRECTORIES
                        S.dirout =  fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), ...
                            sprintf('%s_mode%i_rd%0.2f_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2fmudep%0.2e', ...
                            S.DATE,S.scanmode,mrd,S.initmode,wrd,initialvelocity,initanis,anis,mudep));
                        
                        % make output directory (when it does not exist yet)
                        if ~exist(S.dirout,'dir'), mkdir(S.dirout); end
                        
                        fnameout = sprintf('%s_%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f',beat,type,S.initmode,wrd,initialvelocity,initanis,anis);
                        
                        % create directory and filename to save figures:
                        if exist([S.dirout '/figures'], 'dir') || S.savefigures ~= 1,
                        else
                            mkdir([S.dirout '/figures']);
                        end
                        
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
                            surfnames = {'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'};
                            
                            S.diroutnew = S.dirout;
                            ddir        = dir(fullfile(S.dirout, [fnameout '*init.mat']));
                            
                            % what is happening here??
                            if isempty(ddir),
                                trydir = dir(fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), sprintf('%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f*',S.DATE,S.initmode,wrd,initialvelocity,initanis,anis)));
                                trydir = trydir(cell2mat({trydir.isdir}));
                                for td = trydir',
                                    ddir = dir(fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters),td.name, [fnameout '*init.mat']));
                                    if ~isempty(ddir), S.diroutnew = fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters),td.name); break; end;
                                end
                            end
                            
                            if isempty(ddir),
                                if wrd == 1, wrdin1 = 'wrd1'; wrdin2 = 'wrd0'; else wrdin1 = 'wrd0'; wrdin2 = 'wrd1'; end
                                S.diroutnew = regexprep(S.dirout,wrdin1,wrdin2);
                                ddir        = dir(fullfile(S.diroutnew, [regexprep(fnameout,wrdin1,wrdin2) '*init.mat']));
                            end
                            
                            %% INITITIAL ESTIMATION!
                            if ~S.doinit && ~isempty(ddir),
                                clear measinit;
                                load(fullfile(S.diroutnew,ddir(1).name),'measinit');                        % load previous initial estimation
                            else
                                % initial estimation based on unknown initial activation site [S.initmode = 0] or stimulation focus [S.initmode = 1]:
                                if S.initmode == 0,
                                    mfclusters  = S.clusters;
                                    mfusecor    = 1;
                                    mfissinus   = S.issinus;
                                    mffocus     = -1;
                                else
                                    [stimpos,cathsurf]                                  = GetCatheterPosition2(sub{1},bsmfiles(ibeat).name);    % get stimulation position and surface of stimulation:
                                    if S.initmode > 0                                                                                        	% scan initmode vertices nearest to stimulation site
                                        [pnearest, distp, trinearest, focus, distver]   = findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,S.initmode);
                                    else                                                                                                        % scan initmode vertices nearest to vertex opposite to stimulation site
                                        oppsurf                                         = unique(surfnames(~strcmp(cathsurf,surfnames)));
                                        [pnearest, distp, trinearest, focus, distver]   = findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,-S.initmode);
                                    end
                                    mfclusters  = 1;
                                    mfusecor    = 1;
                                    mfissinus   = 0;
                                    mffocus     = focus;
                                end
                                
                                % actual initial estimation: include scanmode?
                                [measinit.foci, measinit.dep, measinit.outp, corsinit, rdsinit]  =  multifociscan(GEOM,'clusters',mfclusters,'usecor',mfusecor,'issinus',mfissinus, ...
                                    'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'focus',mffocus, 'showplots',showplots,'blmode',S.blmode);
                                
                                
                                
                                measinit.anisotropyRatio                    = GEOM.anisotropyRatio;         % store anisotropy value used for initial estimation
                                measinit.cor                                = measinit.outp(end,1);         % store correlation value for initial estimation
                                measinit.rd                                 = measinit.outp(end,2);         % store RD value for initial estimation
                                measinit.rep                                = initRep(GEOM,measinit.dep);   % determine repolarization times
                                disp(['min deptime ' num2str(min(measinit.dep))])
                            end
                            
                            close all;
                            
                            %% Visualize the correlation and RD values for the first focus:
                            if S.init_cor_rd == 1,
                                
                                % Define data:
                                data.corsinit   = corsinit;
                                data.rdsinit    = rdsinit;
                                data.measinit   = measinit;
                                qtrimode        = 1;
                                figname         = [S.dirout '/figures/init_' fnameout_fig];
                                                                
                                % Call visualization function:
                                invqtri(qtrimode, figname, GEOM, data)
                            end
                            
                            %% INVERSE CALCULATION!
                            if ~S.doinv % if no inverse is required, copy initial guess to inverse results
                                
                                meas.depfinal   = measinit.dep;
                                meas.repfinal   = measinit.rep;
                                meas.corfinal   = measinit.cor;
                                meas.rdfinal    = measinit.rd;
                                meas.log        = '';
                                meas.iterfinal  = NaN;
                                
                            else
                                % The actual inverse calculation:
                                meas = inversefunc(GEOM,measinit.dep,measinit.rep,'estimateampl',estamp,'casedir',S.dirout,...
                                       'repopt','apd','maxiter',maxit,'mudep',mudep,'murep',murep,'minrd',mrd,'mode',S.scanmode,...
                                       'weighedrd',wrd,'blmode',S.blmode,'lpass',S.lpass);
                                
                                RESULTS = [RESULTS; GEOM.VER(meas.depfinal==min(meas.depfinal),:)];                                                     % write results
                                t       = 0:GEOM.SPECS.qrstduration-1;                                                                                  % construct time-vector QRST
                                T       = ones(length(GEOM.VER),1)*t;                                                                                   % matrix vertices mesh x time vector
                                
                                % Visualize the results of the inverse calculations:
                                if showplots,
                                    
                                    % Define data:
                                    data            = [];
                                    data.T          = T;
                                    data.meas       = meas;
                                    data.measinit   = measinit;
                                    
                                    % Visualize ECG's and save results with sigplot
                                    plotecg(GEOM, data, [S.dirout '/figures/' fnameout_fig], S.scanmode, S.lpass)
                                    
                                    % Define data:
                                    data            = [];
                                    data.meas       = meas;
                                    data.measinit   = measinit;
                                    qtrimode        = 2;
                                    figname         = [S.dirout '/figures/focus_' fnameout_fig];
                                    
                                    % Visualize activation times and save results with qtriplot
                                    invqtri(qtrimode, figname, GEOM, data)
                                    
                                    if S.scanmode == 4, figname = [S.dirout '/figures/rep_' fnameout_fig]; invqtri(qtrimode, figname, GEOM, data, S.scanmode); end

                                end
                            end
                            
                            close all;
                            
                            % Store parameters used to produce results:
                            meas.velocity       = '';           % File velocity
                            meas.anis           = anis;         % File anisotropy value 
                            meas.mudep          = mudep;        % File mudep value
                            meas.murep          = murep;        % File murep value
                            meas.specs          = GEOM.SPECS;   % File specs
                            measinit.specs      = GEOM.SPECS;   % File specs
                            measinit.velocity   = '';           % 
                            measinit.anis       = '';           % 
                            
                            % Save results as matlab files:
                            if saveresult,
                                save(fullfile(S.dirout, [fnameout '_' S.DATE 'init.mat']),'measinit', 'GEOM')
                                save(fullfile(S.dirout ,[fnameout '_' S.DATE '.mat']),'meas', 'GEOM')
                                if S.dostampl, save(fullfile(S.dirout, [GEOM.subject '_' beat '_' type '_IST' num2str(ist) '_' num2str(anis) '_' S.DATE 'stampl.mat']),'Ampl'); end
                            end                            
                        end
                    end
                end
            end
        end
    end
end