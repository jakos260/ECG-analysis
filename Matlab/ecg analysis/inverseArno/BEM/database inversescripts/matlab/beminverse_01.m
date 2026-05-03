%% %% MAIN SCRIPT FOR INVERSE ANALYSES %%
%
% Adaptation of script CreatePigs0_2.m , version 6-jan-2015, PvO.
% Version 1.0 AJ


%% Initialize general parameters inverse analyses %%
clearvars

measwrd = [];   % global measwrd
measrd  = [];   % global measrd

close all

[bp]            = bempath;                  % main directory BEM analyses
S.geomdir       = [bp 'Data/Geometries/'];  % directory geometries
S.refLeads      = 1;                        % [1] zeromean BSM and A-matrix, [2] apply WCT to A-matrix and BSM, [-2] 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured.
casedir         = '';
S.dirout        = fullfile('./results/');
S.lpass         = 5;                        % # samples in lowpassma used to filter fw results and BSM!!!
S.clusterDist   = 30;                       % used in multifociscan

% which leads should be used in the initial esitimate
leadset         = 'all';
group           = [];
sub             = {9};      % {'b01', 'b03'} or { 9 10 7 8} or ... , can be an array
saveCase        = 1;        % [2] rewrite ecgspecs, [1] read ecgspecs when available
saveresult      = 1;        % [1] save results of initial estimate and optimization
S.skipexisting  = 0;        % [1] skip if result files already present in target folder (recover from crash)
S.doinit        = true;     % if not true read initial estimate from file

if ~S.doinit, S.skipexisting = false; end

S.doinv         = true;     % perform inverse after intial estimate
S.interpol      = 0;        % [1] removed leads are interpolated, no interpolation for bucket

showplots       = true;
plotegm         = true;     % plot egm on catherer positions
S.sinkScan      = 0;
S.clusters      = 1;        % Number of initial activation locations [foci] determined for inverse calculation
wrd             = 0;        %[0,1]; % use weighed rd in inverse
S.useScaling    = 1;
S.blmode        = 1;        % [1] blcorrection qrsonset-T, [2] (Vstim-10)-T. Applied in multifociscan_AJ.m

%initialtype
S.initmode      = 0 ;       % methods for intial estimate: [0] use normal fociscan, [1] test with fixed initial estimate focus at stimulation site, [n] (n = 'number') test nearest n points at same surface, 'negative number': test with focus on opposite surface
S.issinus       = 0;        % will change beat folder also

type            = 'ventricles'; %% What is this????

S.dostampl      = false;    % fit AP amplitude on ST segment.

S.DATE          = datestr(now,'yyyymmdd');

%% INVERSE ANALYSES FOR ALL SUBJECTS%%
for sub = sub % for all defined subjects
    
    %% set initial velocities and anisotropies, and create/define directories
    S.INITIALVELOCITY   = [0,0.4,1.0];  % velocity of first part QRS complex, vector with undefined length
    S.maxvelocity       = 2.5;          % maximum propagation velocity
    velo                = 2.0;          % ??
    S.ANIS              = 1; %[0.01,1,2.5]; % anisotropy in propagation of activation [the higher the value, the slower the propagation!], vector with undefined length. 1 = isotrope
    S.INITANIS          = [0.01,1,2.5]; % anisotropy of first part QRS complex, vector with undefined length
    
    [subject, bsmdir,layfile]   = subjectdir(sub, bp, S.issinus);       % set subject filename, bsm directory and layfile
    S.dirout                    = [bp 'Data/results/' subject '/'];     % set output results directory
    S.diroutOrg                 = S.dirout;
    if strcmp(subject,''), return; else end
    bsmfiles                    = dir(fullfile(bsmdir, '*.selecg'));    % select all bsm files
    
    RESULTS = []; % clear results structure
    
    %% COMPLETE INVERSE ANALYSES FOR 1 SUBJECT: ALL BEATS%%
    tic % time analyses: ends with toc.
    for ibeat = 1 : length(bsmfiles),
        %% INVERSE ANALYSES FOR 1 SUBJECT, 1 BEAT: posible multiple de/repolarization times%%
        for wrd = wrd,
            if wrd,
                MUDEP = 5e-5; % can be a vector!
                murep = 5e-5;
            else
                MUDEP = 5e-6;
                murep = 5e-5;
            end
            
            clear meas;         % why here?
            clear measinit;     % ?
            
            bsmfile = fullfile(bsmdir, bsmfiles(ibeat).name);   % set pathname of #ibeat from bsmfiles
            beat    = bsmfiles(ibeat).name(1:end-7);            % select #ibeat from bsmfiles and remove (.selecg)
            orgname = beat(1:strfind(beat,'beat')-2);           % remove (_beat...) from beat
            
            if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file')
                T       = load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat']);
                remove  = T.DATA.remove(1:64);                  % why only 64???
                clear T;
            else
                error('Export file not available');
            end
            
            % remove additional data for specific datasets!!
            if ~ischar(sub{1}), warning('NOT skipping front legs'); else remove([32 63 64 ]) = 1; end
            
            %% Load all data from selected subject
            disp('===============================================================')
            disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
            GEOM = invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',S.ANIS(1),'group',group,'basedir',S.geomdir);
            
            S.dirout    = S.diroutOrg; % ??
            heartpath   = fullfile(S.geomdir, group,  subject,  subject );
            
            % create directory for subject in results folder:
            if ~exist(S.dirout,'dir'), mkdir(S.dirout); end
            
            %% save settings:
            if saveresult, save(fullfile(S.dirout, [subject '_' S.DATE '_Settings.mad']),'S'); end
            
            %%
            GEOM.beat = beat;
            
            useLeads = find(remove==0)';    % BSM leads to be used
            remLeads = find(remove==1);     % BSM leads to be removed
            
            % interpolate OR remove missing leads!
            if S.interpol
                warning('S.interpolating missing leads!');
                T                       = intripol(GEOM.TVER,GEOM.TITRI,useLeads);
                Trem                    = T(remLeads,:);
                GEOM.BSM(remLeads,:)    = Trem*GEOM.BSM(useLeads,:);
                remove(remLeads)        = 0;
                useLeads                = find(remove==0)';
                remLeads                = find(remove==1);
            end
            
            % removing bad leads
            L                                   = GEOM.LAY(2:end,:);
            L(ismember(L(:,3),find(remove)),:)  = [];
            for i = 1:size(L,1), L(i,3) = L(i,3) - sum(remove (1:L(i,3))); end  % BSM numbering after removing channels
            GEOM.LAY                            = [GEOM.LAY(1,:); L];
            
            % adjust GEOM.BSM and use zeromean or wct
            GEOM                                = selectLeads(GEOM,useLeads,S.refLeads);
            
            % determine parameters AP: e.g. activation, repolarization
            GEOM.SPECS                          = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim',1);
            
            % fit AP amplitude on ST segment.
            if S.dostampl,
                % Amplitude
                rmsbsm          = rms(GEOM.BSM);
                [ap_min,ist]    = min(rmsbsm((25:100)+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint));
                ist             = ist+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint;
                figure; plot(rms(GEOM.BSM))
                hold on
                plot(ist,rmsbsm(ist),'r*');
                %         ist=350; % defaultposition for ST amplitude estimation.
                Ampl = sqrAmpl(GEOM,ist);
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
                                    if wrd,
                                        S.diroutnew = regexprep(S.dirout,'wrd1','wrd0');
                                        ddir        = dir(fullfile(S.diroutnew, [regexprep(fnameout,'wrd1','wrd0') '*init.mat']));
                                    else
                                        S.diroutnew = regexprep(S.dirout,'wrd0','wrd1');
                                        ddir        = dir(fullfile(S.diroutnew, [regexprep(fnameout,'wrd0','wrd1') '*init.mat']));
                                    end
                                end
                                
                                %% INITITIAL ESTIMATION!
                                if ~S.doinit && ~isempty(ddir),
                                    clear measinit;
                                    load(fullfile(S.diroutnew,ddir(1).name),'measinit');
                                else
                                    % initial estimation based on .... 2 options!
                                    if S.initmode == 0,
                                        [measinit.foci,measinit.dep,measinit.outp]  = multifociscan_AJ(GEOM,'clusters',S.clusters,'usecor',1,'issinus',S.issinus,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'showplots',showplots,'blmode',S.blmode); % actual initial estiation
                                    else
                                        [stimpos,cathsurf]                          = GetCatheterPosition2(sub{1},bsmfiles(ibeat).name); % --> .m-file not present!!
                                        if S.initmode > 0 % scan initmode vertices nearest to stimulation site
                                            [pnearest, distp, trinearest, focus, distver]   = findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,S.initmode);
                                        else % scan initmode vertices nearest to vertex opposite to stimulation site
                                            oppsurf                                         = unique(surfnames(~strcmp(cathsurf,surfnames)));
                                            [pnearest, distp, trinearest, focus, distver]   = findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,-S.initmode);
                                        end
                                        [measinit.foci,measinit.dep,measinit.outp] = multifociscan_AJ(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'focus',focus,'showplots',showplots,'blmode',S.blmode);
                                    end
                                    
                                    % 
                                    mfclusters  = S.clusters;
                                    mfusecor    = 1;
                                    mfissinus   = S.issinus;
                                    mffocus     = -1;
                                    
                                    %
                                    mfclusters  = 1;
                                    mfusecor    = 1;
                                    mfissinus   = 0;
                                    mffocus     = focus;
                                    
                                    [measinit.foci,measinit.dep,measinit.outp] = multifociscan_AJ(GEOM,'clusters',mfclusters,'usecor',mfusecor,'issinus',mfissinus,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'focus',mffocus,'showplots',showplots,'blmode',S.blmode);
                                    
                                    disp(['min deptime ' num2str(min(measinit.dep))])
                                    
                                    measinit.anisotropyRatio    = GEOM.anisotropyRatio;
                                    measinit.cor                = measinit.outp(end,1);
                                    measinit.rd                 = measinit.outp(end,2);
                                    measinit.rep                = initRep(GEOM,measinit.dep);   % what is happening here?
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
                                    pp = 30; % ??
                                    
                                    meas = inverse_aj(GEOM,measinit.dep,measinit.rep,'estimateampl',0,'casedir',S.dirout,...
                                        'repopt','apd','maxiter',400,'mudep',mudep,'murep',murep,'minrd',0.15,'mode',1,...
                                        'weighedrd',wrd,'blmode',S.blmode,'lpass',S.lpass);
                                    
                                    % visualize the results of the inverse calculations:
                                    if showplots,
                                        figure('Name','Final activation'); 
                                        showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',[find(meas.depfinal==min(meas.depfinal(GEOM.typ==1))), find(meas.depfinal==min(meas.depfinal(GEOM.typ==2)))]); 
                                        view(0,0);
                                        
                                        figure('Name','Final repolarization');
                                        showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal)
                                        
                                        figure('Name','Final ARI');
                                        showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal)
                                    end
                                    
                                    % writeresults;
                                    RESULTS = [RESULTS; GEOM.VER(meas.depfinal==min(meas.depfinal),:)];

                                    t = 0:GEOM.SPECS.qrstduration-1;
                                    T = ones(length(GEOM.VER),1)*t;
                                    
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
                                        
                                        % hack nothing before epi breakthrough
                                        % initdep(initdep<min(initdep(GEOM.typ==1))) = min(initdep(GEOM.typ==1)); % nothing before epi
                                        
                                        initrep                                 = initRep(GEOM,initdep,0);
                                        PSIAinitnohakkel                        = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,1, GEOM),S.lpass);
                                        PSIAinit                                = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,1),S.lpass);
                                        % PSIAinitendowait                        = lowpassma(GEOM.AMA*getSmode(T,initdep,initrep,GEOM.SPECS,4, GEOM),S.lpass);
                                        % load('/Users/peteroosterhoff/Documents/MATLAB/Emmer/actint.mat','ActInt');
                                        % PSIANeedle                              = lowpassma(GEOM.AMA*getSmode(T,ActInt,ActInt+meas.repfinal-meas.depfinal,GEOM.SPECS,4,GEOM),S.lpass);
                                        
                                        figure('Name','Final QRST');
                                        clf;
                                        sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),'',GEOM.LAY,1.5,'b',1,0);
                                        hold on
                                        sigplot(PSIA,'',GEOM.LAY,1.5,'r',1,0);
                                        % sigplot(PSIANeedle,'',GEOM.LAY,1.5,'k',1,0);
                                        
                                        figure('Name','Final QRS');
                                        clf;
                                        sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:(GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint-1)),'',GEOM.LAY,1.5,'b',1,0);
                                        hold on
                                        sigplot(PSIA(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'r',1,0);
                                        % sigplot(PSIANeedle(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'k',1,0);
                                        sigplot(PSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'m',1,0);
                                        hold on
                                        sigplot(AltPSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'k',1,0);
                                        % sigplot(PSIAinitendowait(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'m:',1,0);
                                        
                                        % if exist('PSIAprev','var'),
                                        %    sigplot(PSIAprev,'',GEOM.LAY,1.5,'g',1,0);
                                        % end
                                        
                                        PSIprev = PSIA;
                                        hold on
                                        
                                        if 0,
                                            figure
                                            plot(rms(GEOM.BSM(:,GEOM.SPECS.onsetqrs:end)));
                                            hold all
                                            plot(rms(PSIAinit));
                                            figure('Name','Inverse/bsm');
                                            bsmb        = GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave-1);
                                            ratio       = mean(PSIA./bsmb,1);
                                            plot(ratio,'-r');
                                            hold on
                                            ratioinit   = mean(PSIAinit./bsmb,1);
                                            plot(ratioinit,'-g');
                                            fprintf('average ratio init:%0.1f, Inv:%0.1f\n',mean(ratioinit),mean(ratio));
                                            ylim([-3 3]);
                                        end
                                        
                                        if 0, % endo epi plot init, altinit, depfinal
                                            qtriplot({{'reset'};{'horizontal 2'};{'vertical 3'};{'funcolor tim'};{'mouse fun'}});
                                            EndoEpiPlot(GEOM,ActInt,'needles',[1,1;2,1],'BOTH');
                                            EndoEpiPlot(GEOM,measinit.dep,'init',[1,1;2,1],'BOTH');
                                            EndoEpiPlot(GEOM,altdep,'altinit',[1,2;2,2],'BOTH');
                                            EndoEpiPlot(GEOM,meas.depfinal,'final',[1,3;2,3],'BOTH');
                                        end
                                        
                                        [mindep,imindep]    = min(meas.depfinal);
                                        idepearly           = find(meas.depfinal<=(mindep+2) & meas.depfinal>(mindep));
                                        qtriplot({{'reset'};{'horizontal 2'};{GEOM.VER,GEOM.ITRI,'init'};{measinit.dep};...
                                            {GEOM.VER(measinit.foci,:),[],'fociinit'};{'color white'};...
                                            {GEOM.VER,GEOM.ITRI,'final'};{meas.depfinal};{GEOM.VER(imindep,:),[],'focusfinal'};{'color green'};...
                                            {GEOM.VER(idepearly,:),[],'focifinal'};{'color white'};...
                                            {'panel final=2 1'};{'panel focusfinal=2 1'};{'panel focifinal=2 1'};...
                                            {'funcolor tim'}});
                                        pause
                                                                                
                                    end
                                    %     PSIA              = PSIA(:,1:GEOM.qrsduration);
                                    %     BSM               = GEOM.BSM(:,GEOM.specs(2):GEOM.specs(2)+size(PSIA,2)-1);
                                    %     meas.rdfinalqrs   =  norm(BSM - PSIA,'fro')/norm(BSM,'fro');
                                    
                                end
                                
                                meas.velocity       = velo;
                                meas.anis           = anis;
                                meas.mudep          = mudep;
                                meas.murep          = murep;
                                meas.specs          = GEOM.SPECS;
                                measinit.specs      = GEOM.SPECS;
                                % meas.vt             = vt;
                                measinit.velocity   = velo;
                                measinit.anis       = anis;
                                % measinit.vt         = vt;
                                
                                if saveresult,
                                    save(fullfile(S.dirout, [fnameout '_' S.DATE 'init.mat']),'measinit')
                                    save(fullfile(S.dirout ,[fnameout '_' S.DATE '.mat']),'meas')
                                    if S.dostampl, save(fullfile(S.dirout, [GEOM.subject '_' beat '_' type '_IST' num2str(ist) '_' num2str(anis) '_' S.DATE 'stampl.mat']),'Ampl'); end
                                end
                                
                                % A             = GEOM.VER(find(meas.depfinal == min(meas.depfinal)),:);
                                % xlswrite(GEOM.subject,A,['velocity' num2str(velo) 'anis' num2str(anis)],['A' num2str(ibeat)]);
                                
                                % geomheart     = 'C:\Inge\invedl\dataVarken08\model\Pig08\Pig08_ventricles.tri';
                                % [VER, ITRI]   = loadtri(geomheart);
                                % [p,q]         = min(meas.depfinal);
                                % xlswrite(num2str(GEOM.beat),[num2cell(GEOM.beat),' ',num2str(measinit.rd),' ',num2str(measinit.cor),' ',num2str(measinit.foci),' ',num2str(GEOM.VER(measinit.foci,:)),' ',num2str(meas.rdfinal),' ',num2str(meas.corfinal),' ',q,' ',num2str(GEOM.VER(q,:))]);
                                
                                if wrd, measwrd = meas; else measrd = meas; end
                                if ~isempty(measwrd) && ~isempty(measrd), Comparefun(measrd.depfinal,measwrd.depfinal,GEOM.VER,GEOM.ITRI); end
                                
                            end
                        end
                    end
                end
            end
        end
    end
    toc
    
    %% VISUALIZE RESULTS ... WITH QTRIPLOT
    if 0, % replace with main parameter at top of this script!
        qtriplot({{'reset'} {'horizontal 2'} {'vertical 1'}});
        reg     = zeros(length(GEOM.VER),1);
        regop   = surflapl(GEOM.VER/1000,GEOM.ITRI,1);
        for i = 1:length(GEOM.VER), reg(i) = rms(regop*GEOM.DIST25(:,i)/0.8); end
        qtriplot(GEOM.VER,GEOM.ITRI,'v0.8');
        qtriplot(reg,'v0.8');
        qtriplot('panel v0.8 = 1 1');
        for i = 1:length(GEOM.VER), reg(i) = rms(regop*GEOM.DIST25(:,i)/1.5); end
        qtriplot(GEOM.VER,GEOM.ITRI,'v1.5');
        qtriplot(reg,'v1.5');
        qtriplot('panel v1.5= 2 1');
    end
end