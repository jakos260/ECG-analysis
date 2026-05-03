clearvars
% global measrd measwrd
measwrd=[];
measrd=[];

close all
% S.geomdir = fullfile('C:\Users\Damp2\Documents\ECG_simulation\STW\Data\geometries\');
S.geomdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/';

S.refLeads=1; % 1: zeromean BSM and A-matrix, -2: 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured. 2: apply WCT to A-matrix and BSM.
casedir		='';
S.dirout = fullfile('./results/');
S.lpass=5%;   % # samples in lowpassma used to filter fw results and BSM!!!
S.clusterDist=30; % used in multifocisacn

% which leads should be used in the initial esitimate
leadset='all';
group=[];
sub={9};%{'b01', 'b03'}%{ 9 10 7 8}; % can be an array now
saveCase =1%1; % 2 = overschrijven ecgspecs, 1 read ecgspecs when available
saveresult=1%1; % save results of initial estimate and optimization
S.skipexisting=0; % skip if result files already present in target folder ( recover from crash)
S.doinit=true; % if not true read initial estimate from file
if ~S.doinit
    S.skipexisting=false;
end
S.doinv=true % perform inverse after intial estimate
S.interpol=0; % if 1 removed leads are interpolated, no interpolation for bucket

showplots=true;
plotegm=true; % plot egm on catherer positions
S.sinkScan= 0;
S.clusters = 1%6;%16; % 4
wrd=0%[0,1]; % use weighed rd in inverse
S.useScaling =1;
S.blmode=1%; % oostep1: 1 blcorrection qrsonset-T, 2: (Vstim-10)-T.
%initialtype
S.initmode=0 % methods for intial estimate
% 0: use normal fociscan
% 1: test with fixed initial estimate focus at stimulation site
% n ('number'): test nearest n points at same surface
% negative: test with focus on oposite surface
S.issinus=1 % % will change beat folder also

layfile='pig_adam.mla';%'pig64.mla';
% layfile='Bucket01Removed.lay';
type='ventricles';

S.dostampl=false; % fit AP amplitude on ST segmet.
% initialActTime=0; % override. 0 means scaling with maxvelocity and ANIS. Make [] for automatic setting in multifociscan

% mudep=0
% mudep=mudep/33

S.DATE=datestr(now,'yyyymmdd');
% S.DATE='20141031' % overrule when continuing earlier analysis

 for sub=sub 
% if 1 % for only works for number arrays


% S.INITIALVELOCITY=[0.2,0.4,0.6,0.8,1.0];
S.INITIALVELOCITY=[0,0.4,1.0];
S.INITIALVELOCITY=[0]%0.4;

S.maxvelocity=[2.5];%1.6;
velo = 2.0;

S.ANIS = [0.01,1,2.5];
S.ANIS=[0.01 2.0]%,2,2.5]%,2];
S.ANIS=1.0;
% S.ANIS=[0.01]

S.INITANIS=[0.01,1,2.5];
S.INITANIS=1.0;%2.5;
% S.INITANIS=2.5


% if sub{1}==2
%     subject='ARVCp02';
%     %     bsmfile=fullfile('C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File01_beat7.ecg');beat='01';
%     bsmdir='/Users/peteroosterhoff/Documents/Werk/ARVC/DATA/ARVCp02/ECG/export/beats/';
%     S.dirout= ['/Users/peteroosterhoff/Documents/Werk/ARVC/DATA/results/' subject '/'];
%     S.refLeads=-2;
% else
if sub{1} == 7
    subject='Pig07';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/AVG/beats/';
    %     S.dirout= ['.\results\' subject '\'];
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
elseif sub{1} == 8
    subject='Pig08';
%     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/Pacing/AVG/beats/';
bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/AVG/beats/'; % dummy prepare for SR
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
elseif sub{1} == 9
    %         subject='Pig09r';
    subject='Pig09';
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats3/';BAH,
    %     NOOIT MEER DOEN!!
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
    %         bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/SR/AVG/beats/';
    
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/SRST/beats/';
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    %     bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\avobeats\';
    %     S.dirout= ['.\results\' subject '\avobeats\'];
elseif sub{1} == 10
    subject='Pig10';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/Biosemi/export/AVG/beats/';
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    
elseif strcmp(sub{1},'b01')
    subject='Bucket01';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/Biosemi/export/beats/';
%     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/Biosemi/export20141021/beats/';
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    layfile='bucket.lay';
    S.blmode=1;
elseif strcmp(sub{1},'b01s')
    subject='Bucketsmall01';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/Biosemi/export/beats/';
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    layfile='bucket.lay';
    S.blmode=1;
elseif strcmp(sub{1},'b03')
    subject='Bucket03';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket03/export/beats/';
    S.dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    layfile='bucket.lay';
    S.blmode=1;
    
else
    return
end

if S.issinus
    bsmdir=regexprep(bsmdir,'/export/','/export/SR/');
    % './biosemi/export/SR/beats' of './biosemi/export/SR/AVG/beats'
end
bsmfiles = dir(fullfile(bsmdir, '*.selecg'));


RESULTS =[];


S.diroutOrg = S.dirout;




tic
for ibeat = 1 : length(bsmfiles)
    %     for S.initmode=[1,-1]
    for wrd=0%wrd=wrd
        if wrd
            MUDEP=5e-5;%2.5e-6;
            murep=5e-5;%2.5e-6;
        else
            MUDEP=[ 5e-6,]%5e-6 1e-6 5e-7];%[5e-5,1e-5,5e-6];%5e-6;
            murep=5e-5;
        end
        
        clear meas;
        clear measinit;
        bsmfile = fullfile(bsmdir, bsmfiles(ibeat).name);
        beat = bsmfiles(ibeat).name(1:end-7);
        orgname = beat(1:strfind(beat,'beat')-2);
        
        
        
        
        if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file')
            T=load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat']);
            remove = T.DATA.remove(1:64);
            clear T;
        else
            error('Export file not available');
            %         A=loadmat(bsmfile);
            %         remove = zeros(64,1);
        end
        
        
        if ~ischar(sub{1})
%             remove([31,63])=1;
            warning('NOT skipping front legs');
        else
            remove([32 63 64 ])=1% 13 14 19 20 25])=1; % bucket03 ECG with P-wave. WHY?
        end
        
        
        
        
        % Load all data from selected subject
        disp('===============================================================')
        disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
        %
        GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',S.ANIS(1),'group',group,'basedir',S.geomdir);
        
        
        S.dirout = S.diroutOrg;
        heartpath=fullfile(S.geomdir, group,  subject,  subject );
        
        if ~exist(S.dirout,'dir')
            mkdir(S.dirout)
        end
        
        if saveresult
            save(fullfile(S.dirout, [subject '_' S.DATE '_Settings.mad']),'S')
        end
        
        
        GEOM.beat=beat;
        %     useLeadsold=1:size(GEOM.LAY,1) - 1;
        %     useLeadsold(remove ==1)=[];
        
        %         remove(21)=1 % test
        useLeads=find(remove==0)'; % original BSM numbering, including removed channels
        remLeads=find(remove==1);
        if ischar(sub{1})
            remLeads(remLeads==32 | remLeads==63 |remLeads==64)=[]; % Bucket unused channels
        end
 
        
        if S.interpol
            warning('S.interpolating missing leads!');
            T=intripol( GEOM.TVER,GEOM.TITRI,useLeads);
            Trem=T(remLeads,:);
            GEOM.BSM(remLeads,:)=Trem*GEOM.BSM(useLeads,:);
            remove(remLeads)=0;
            useLeads=find(remove==0)'; % original BSM numbering, including removed channels
            remLeads=find(remove==1);
        end
        
        
        L= GEOM.LAY(2:end,:);
        L(ismember(L(:,3),find(remove)),:)=[];
        for i=1:size(L,1)
            %         L(i,3) = L(i,3) - sum(remove (1:i));
            L(i,3) = L(i,3) - sum(remove (1:L(i,3))); % BSM numbering after removing channels
        end
        
        GEOM.LAY = [GEOM.LAY(1,:); L];
        
        
        %     useLeads= [40 41 42 46 47 48 52 53 54];
        %         useLeads= [38 40 42 43 45 47 49 51 53];
        % %     useLeads = round(rand(9,1) * (length(useLeads)-1));
        %     GEOM.LAY= loadmat('9lds.mla');
        GEOM = selectLeads(GEOM,useLeads,S.refLeads);
        GEOM.SPECS = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim',1);
%         GEOM.SPECS.onsetqrs=GEOM.SPECS.onsetqrs-10
%         GEOM.SPECS.qrsduration=GEOM.SPECS.qrsduration+10
%         GEOM.SPECS.time_Jpoint=GEOM.SPECS.time_Jpoint+10
        %                 continue % when just re-annotating ecg beats
        %         GEOM = selectLeads(GEOM,useLeads,1);
        
        
        %     GEOM = prepare_geom(GEOM,[bsmfile(1:end-7) '.spe'],saveCase,GEOM.ECGextra );
        
        %                 GEOM.BSM(:,GEOM.specs(2):end) = killhum(GEOM.BSM(:,GEOM.specs(2):end),50,1000,0.05);
        %                 GEOM.BSM = lowpassma(GEOM.BSM,20);
        
        %     GEOM.specs(6)=0.005513;
        %     GEOM.specs(7)=0.055512;
        %     GEOM.specs(8)=1.317503;
        %     GEOM.pS = GEOM.specs(6:8);
        
        
        %     eval(['BSM' num2str(ibeat) '=baselinecor(GEOM.BSM(:,GEOM.SPECS.q(2):end),1,GEOM.specs(5)-GEOM.specs(2));']);
        
        if S.dostampl
            % Amplitude
            rmsbsm=rms(GEOM.BSM);
            [~,ist]=min(rmsbsm((25:100)+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint));
            ist=ist+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint;
            figure;plot(rms(GEOM.BSM))
            hold on
            plot(ist,rmsbsm(ist),'r*');
            %         ist=350; % defaultposition for ST amplitude estimation.
            Ampl=sqrAmpl(GEOM,ist);
            if showplots
                qtriplot('reset');
                qtriplot(GEOM.VER,GEOM.ITRI,'Ventricles');
                qtriplot(Ampl);
                qtriplot([-16.8	27.3	888.3],[],'LVLatENdoII');
            end
        end
        
        
        
        close all;
        
        for initialvelocity=S.INITIALVELOCITY;
            for initanis=S.INITANIS;
                for anis=S.ANIS
                    for mudep=MUDEP
                        clear meas;
                        clear measinit;
                        if 1%S.clusters == 1 % 25 mm
                            %                     S.dirout = fullfile(S.dirout, 'single', sprintf('%s_im%d_S.wrd%d_iV%0.1f_iS.ANIS%.01fS.ANIS%0.1f',S.DATEstr(now,'yyyymmdd'),S.initmode,wrd,initialvelocity,initanis,anis));
                            S.dirout = fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), sprintf('%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2fmudep%0.2e',S.DATE,S.initmode,wrd,initialvelocity,initanis,anis,mudep));
                        end
                        if ~exist(S.dirout,'dir')
                            mkdir(S.dirout)
                        end
                        
                        %                 fnameout=[subject '_' beat '_' type '_' num2str(velo,'%0.1f') '_' num2str(anis,'%0.2f')];
                        fnameout=sprintf('%s_%s_%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f',subject,beat,type,S.initmode,wrd,initialvelocity,initanis,anis);
                        
                        if ~S.skipexisting || length(dir(fullfile(S.dirout, [fnameout '*.mat'])))<2
                            
                            if initanis==1
                                GEOM.DIST25=GEOM.DIST;
                            else
                                GEOM.DIST25=loadmat(fullfile(S.geomdir,subject,[subject num2str(initanis) '_' type '.dist']));
                                if isempty(GEOM.DIST25)
                                    error(' requested distance file not found');
                                end
                            end
                            if anis==1
                                GEOM.DIST2W=GEOM.DIST;
                            else
                                GEOM.DIST2W=loadmat(fullfile(S.geomdir,subject,[subject num2str(anis) '_' type '.dist']));
                            end
                            GEOM.anisotropyRatio=anis;
                            
                            %         drawnow
                            
                            %         S.clusters =1;
                            
                            % initial scan
                            disp(['anisotropyRatio: ' num2str(GEOM.anisotropyRatio)])
                            surfnames={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % Simplified
 
                            S.diroutnew=S.dirout;
                            ddir=dir(fullfile(S.dirout, [fnameout '*init.mat']));
                            if isempty(ddir)
                                trydir=dir(fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters), sprintf('%s_im%d_wrd%d_iV%0.1f_iANIS%.01fANIS%0.2f*',S.DATE,S.initmode,wrd,initialvelocity,initanis,anis)));
                                trydir=trydir(cell2mat({trydir.isdir}));
                                for td=trydir'
                                    ddir=dir(fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters),td.name, [fnameout '*init.mat']));
                                    if ~isempty(ddir)
                                        S.diroutnew=fullfile(S.diroutOrg, sprintf('cluster%d',S.clusters),td.name);
                                        break
                                    end;
                                end
                            end
                            if isempty(ddir)
                                if wrd
                                    S.diroutnew=regexprep(S.dirout,'wrd1','wrd0');
                                    ddir=dir(fullfile(S.diroutnew, [regexprep(fnameout,'wrd1','wrd0') '*init.mat']));
                                else
                                    S.diroutnew=regexprep(S.dirout,'wrd0','wrd1');
                                    ddir=dir(fullfile(S.diroutnew, [regexprep(fnameout,'wrd0','wrd1') '*init.mat']));
                                end
                            end
                            
                            if ~S.doinit && ~isempty(ddir)
                                clear measinit
                                %                             fnamefilt=regexprep(fnameout,'wrd?','wrd*');
                                load(fullfile(S.diroutnew,ddir(1).name),'measinit');
                            else
                                %                         for iinit=1:max(abs(S.initmode),1)
                                if S.initmode==0
                                    %     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_wall(GEOM,'clusters',S.clusters);
                                    %     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_v12(GEOM,'clusters',S.clusters,'usecor',1,'usetime',60,'velocity',0.5);
                                    %     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_pvd(GEOM,'clusters',S.clusters,'usecor',1,'usetime',60,'velocity',0.4);
                                    %     [measinit.foci,measinit.dep,measinit.outp] = multifociscan_pvd(GEOM,'S.clusters',6,'usecor',1,'issinus',1,'initialvelocity',0.4,'maxvelocity',0.9);
                                    %                 [measinit.foci,measinit.dep,measinit.outp] = multifociscan_esc2013(GEOM,'clusters',S.clusters,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'initialacttime',initialActTime);
                                    [measinit.foci,measinit.dep,measinit.outp] = multifociscan_publicationPigsPO(GEOM,'clusters',S.clusters,'usecor',1,'issinus',S.issinus,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'showplots',showplots,'blmode',S.blmode);
%                                     [measinit.foci,measinit.dep,measinit.outp] = multifociscan_publicationPigs(GEOM,'clusters',S.clusters,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'showplots',showplots,'blmode',S.blmode);
                                else
                                    [stimpos,cathsurf]=GetCatheterPosition2(sub{1},bsmfiles(ibeat).name);
                                    if S.initmode>0 % scan initmode vertices nearest to stimulation site
                                        [pnearest, distp, trinearest, focus, distver]=findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,S.initmode);
%                                         focus=1160 % generate data for paper, Bucket01 ?
%                                         warning('Bucket01 fixed focus');
                                    else % scan initmode vertices nearest to vertex opposite to stimulation site
                                        oppsurf=unique(surfnames(~strcmp(cathsurf,surfnames)));
                                        [pnearest, distp, trinearest, focus, distver]=findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,-S.initmode);
                                    end
                                    [measinit.foci,measinit.dep,measinit.outp] = ...
                                        multifociscan_publicationPigsPO(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'focus',focus,'showplots',showplots,'blmode',S.blmode);
                                    %                     multifociscan_esc2013(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'initdep',focus,'initialacttime',initialActTime);
                                    
                                    %                 measinit.outp=NaN(1,3);
                                end
                                
                                disp(['min deptime ' num2str(min(measinit.dep))])
                                
                                measinit.anisotropyRatio = GEOM.anisotropyRatio;
                                measinit.cor = measinit.outp(end,1);
                                measinit.rd  = measinit.outp(end,2);
                                measinit.rep = initRep(GEOM,measinit.dep);
                                %%
                            end
                            if ~S.doinv
                                meas.depfinal=measinit.dep;
                                meas.repfinal=measinit.rep;
                                meas.corfinal=measinit.cor;
                                meas.rdfinal=measinit.rd;
                                meas.log='';
                                meas.iterfinal=NaN;
                            else
                                
                                pp=30;
                                %                         meas=inverse_pvd(GEOM,measinit.dep ,measinit.rep,'estimateampl',0,...
                                %                             'casedir',S.dirout,...
                                %                             'repopt','apd',... %'rep'
                                %                             'maxiter',40,...
                                %                             'mudep',mudep,...
                                %                             'murep',murep,...
                                %                             'weighed',0,...
                                %                             'minrd',0.15,...
                                %                             'mode',4,... % 1 or 4
                                %                             'weighedrd',wrd);
                                
                                meas=inverse_po(...
                                    ...meas=inverseSimplex_po(...
                                    GEOM,measinit.dep ,measinit.rep,'estimateampl',0,...
                                    'casedir',S.dirout,...
                                    'repopt','apd',... %'rep'
                                    'maxiter',400,...
                                    'mudep',mudep,...
                                    'murep',murep,...
                                    'minrd',0.15,...
                                    'mode',1,... % 1 or 4
                                    'weighedrd',wrd,...
                                    'blmode',S.blmode,...
                                    'lpass',S.lpass);
                                if showplots
                                    figure('Name','Final activation');showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',[find(meas.depfinal==min(meas.depfinal(GEOM.typ==1))) find(meas.depfinal==min(meas.depfinal(GEOM.typ==2)))]);view(0,0);
                                    figure('Name','Final repolarization');showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal)
                                    figure('Name','Final ARI');showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal)
                                end
                                
                                %     writeresults;
                                RESULTS = [RESULTS;   GEOM.VER(meas.depfinal==min(meas.depfinal),:) ];
                                
                                %     figure(pp+10);showpatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',measinit.foci);view(0,0);
                                %     figure(pp+10);showpatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',find(meas.depfinal==min(meas.depfinal)));view(0,0);
                                %%
                                %     figure(pp+20);showpatch(GEOM.VER,GEOM.ITRI,meas.repfinal,'nodes',measinit.foci);view(0,0);
                                %     figure(pp+30);showpatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal,'nodes',measinit.foci);view(0,0);
                                %     figure(pp+31);clf;plot(meas.depfinal(GEOM.Rfreewallver==0),meas.repfinal(GEOM.Rfreewallver==0)-meas.depfinal(GEOM.Rfreewallver==0),'.r');
                                %     hold on;plot(meas.depfinal(GEOM.Rfreewallver==1),meas.repfinal(GEOM.Rfreewallver==1)-meas.depfinal(GEOM.Rfreewallver==1),'.b')
                                %             t=0:GEOM.SPECS.endtwave;
                                t=0:GEOM.SPECS.qrstduration-1;
                                T=ones(length(GEOM.VER),1)*t;
                                if showplots
                                    initdep=measinit.dep;
                                    %find init opposite
                                    [~,verfocus]=min(measinit.dep);
                                    [~, ~, ~, altfocus]=findnearest(GEOM.VER(verfocus,:),GEOM.VER,GEOM.ITRI,GEOM.typ,surfnames{GEOM.typ(verfocus)},-1);
                                    [~,altdep] = multifociscan_publicationPigsPO(GEOM,'clusters',S.clusters,'usecor',1,'issinus',0,'focus',altfocus,'initialvelocity',initialvelocity,'maxvelocity',S.maxvelocity,'showplots',0,'blmode',S.blmode);
                                    altdep(altdep<min(altdep(GEOM.typ==1)))=min(altdep(GEOM.typ==1)); % nothing before epi
                                    altrep=initRep(GEOM,altdep,0);
                                    AltPSIAinitnohakkel =lowpassma(GEOM.AMA*getSmode(T,altdep,altrep,GEOM.SPECS,1, GEOM),S.lpass);

                                    
                                    PSIA =lowpassma(GEOM.AMA*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,1,GEOM),S.lpass);
                                    
%                                     % hack nothing before epi breakthrough
%                                     initdep(initdep<min(initdep(GEOM.typ==1)))=min(initdep(GEOM.typ==1)); % nothing before epi

                                    initrep=initRep(GEOM,initdep,0);
                                    PSIAinitnohakkel =lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,1, GEOM),S.lpass);
                                    PSIAinit = lowpassma(GEOM.AMA*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,1),S.lpass);
%                                     PSIAinitendowait =lowpassma(GEOM.AMA*getSmode(T,initdep,initrep,GEOM.SPECS,4, GEOM),S.lpass);
                                    %                             load('/Users/peteroosterhoff/Documents/MATLAB/Emmer/actint.mat','ActInt');
                                    %                             PSIANeedle =lowpassma(GEOM.AMA*getSmode(T,ActInt,ActInt+meas.repfinal-meas.depfinal,GEOM.SPECS,4,GEOM),S.lpass);
                                    
                                    figure('Name','Final QRST');clf
                                    sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),'',GEOM.LAY,1.5,'b',1,0);
                                    hold on
                                    sigplot(PSIA,'',GEOM.LAY,1.5,'r',1,0);
                                    %                             sigplot(PSIANeedle,'',GEOM.LAY,1.5,'k',1,0);
                                    
                                    figure('Name','Final QRS');clf
                                    sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:(GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint-1)),'',GEOM.LAY,1.5,'b',1,0);
                                    hold on
                                    sigplot(PSIA(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'r',1,0);
                                    %                             sigplot(PSIANeedle(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'k',1,0);
                                    sigplot(PSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'m',1,0);
                                    hold on
                                    sigplot(AltPSIAinitnohakkel(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'k',1,0);
%                                     sigplot(PSIAinitendowait(:,1:GEOM.SPECS.time_Jpoint),'',GEOM.LAY,1.5,'m:',1,0);
                                    


                                    
                                    %                             if exist('PSIAprev','var')
                                    %                                 sigplot(PSIAprev,'',GEOM.LAY,1.5,'g',1,0);
                                    %                             end
                                    PSIprev=PSIA;
                                    hold on
                                    
                                    if 0
                                        figure
                                        plot(rms(GEOM.BSM(:,GEOM.SPECS.onsetqrs:end)))
                                        hold all
                                        plot(rms(PSIAinit));
                                        figure('Name','Inverse/bsm');
                                        bsmb=GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave-1);
                                        ratio=mean(PSIA./bsmb,1);
                                        plot(ratio,'-r');
                                        hold on
                                        ratioinit=mean(PSIAinit./bsmb,1);
                                        plot(ratioinit,'-g');
                                        fprintf('average ratio init:%0.1f, Inv:%0.1f\n',mean(ratioinit),mean(ratio));
                                        ylim([-3 3]);
                                    end
                                    
                                    if 0 % endo epi plot init, altinit, depfinal
                                        qtriplot({{'reset'};{'horizontal 2'};{'vertical 3'};{'funcolor tim'};{'mouse fun'}});
                                        EndoEpiPlot(GEOM,ActInt,'needles',[1,1;2,1],'BOTH');
                                        EndoEpiPlot(GEOM,measinit.dep,'init',[1,1;2,1],'BOTH');
                                        EndoEpiPlot(GEOM,altdep,'altinit',[1,2;2,2],'BOTH');
                                        EndoEpiPlot(GEOM,meas.depfinal,'final',[1,3;2,3],'BOTH');
                                    end
                                    [mindep,imindep]=min(meas.depfinal);
                                    idepearly=find(meas.depfinal<=(mindep+2) & meas.depfinal>(mindep));
                                    qtriplot({{'reset'};{'horizontal 2'};{GEOM.VER,GEOM.ITRI,'init'};{measinit.dep};...
                                        {GEOM.VER(measinit.foci,:),[],'fociinit'};{'color white'};...
                                        {GEOM.VER,GEOM.ITRI,'final'};{meas.depfinal};{GEOM.VER(imindep,:),[],'focusfinal'};{'color green'};...
                                        {GEOM.VER(idepearly,:),[],'focifinal'};{'color white'};...
                                        {'panel final=2 1'};{'panel focusfinal=2 1'};{'panel focifinal=2 1'};...
                                        {'funcolor tim'}});
                                    pause
                                    
                                    if plotegm
                                        plotegm
                                    end
                                    
                                end
                                %     PSIA = PSIA(:,1:GEOM.qrsduration);
                                %     BSM=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(2)+size(PSIA,2)-1);
                                %     meas.rdfinalqrs =  norm(BSM - PSIA,'fro')/norm(BSM,'fro');
                                
                            end
                            
                            meas.velocity = velo;
                            meas.anis=anis;
                            meas.mudep=mudep;
                            meas.murep=murep;
                            meas.specs=GEOM.SPECS;
                            measinit.specs=GEOM.SPECS;
                            %     meas.vt=vt;
                            measinit.velocity = velo;
                            %     measinit.vt = vt;
                            measinit.anis=anis;
                            if saveresult
                                save(fullfile(S.dirout, [fnameout '_' S.DATE 'init.mat']),'measinit')
                                save(fullfile(S.dirout ,[fnameout '_' S.DATE '.mat']),'meas')
                                if S.dostampl
                                    save(fullfile(S.dirout, [GEOM.subject '_' beat '_' type '_IST' num2str(ist) '_' num2str(anis) '_' S.DATE 'stampl.mat']),'Ampl')
                                end
                            end
                            
                            %                 A=GEOM.VER(find(meas.depfinal == min(meas.depfinal)),:);
                            %                 xlswrite(GEOM.subject,A,['velocity' num2str(velo) 'anis' num2str(anis)],['A' num2str(ibeat)]);
                            
                            %                 geomheart = 'C:\Inge\invedl\dataVarken08\model\Pig08\Pig08_ventricles.tri';
                            %                 [VER, ITRI] = loadtri(geomheart);
                            %                 [p,q] =min(meas.depfinal);
                            %                 xlswrite(num2str(GEOM.beat),[num2cell(GEOM.beat),' ',num2str(measinit.rd),' ',num2str(measinit.cor),' ',num2str(measinit.foci),' ',num2str(GEOM.VER(measinit.foci,:)),' ',num2str(meas.rdfinal),' ',num2str(meas.corfinal),' ',q,' ', num2str(GEOM.VER(q,:))    ]);
                            %
                            %             pause(15);
                            if wrd
                                measwrd=meas;
                            else
                                measrd=meas;
                            end
                            if ~isempty(measwrd) && ~isempty(measrd)
                                Comparefun(measrd.depfinal,measwrd.depfinal,GEOM.VER,GEOM.ITRI);
                                %                             pause
                            end
                            %                         end
                        end
                    end
                end
            end
        end
    end
end
toc

if 0
    qtriplot({{'reset'}
        {'horizontal 2'}
        {'vertical 1'}
        });
    
    reg=zeros(length(GEOM.VER),1);
    regop= surflapl(GEOM.VER/1000,GEOM.ITRI,1);
    for i=1:length(GEOM.VER)
        reg(i)=rms(regop*GEOM.DIST25(:,i) / 0.8);
    end
    qtriplot(GEOM.VER,GEOM.ITRI,'v0.8');
    qtriplot(reg,'v0.8');
    qtriplot('panel v0.8= 1 1');
    for i=1:length(GEOM.VER)
        reg(i)=rms(regop*GEOM.DIST25(:,i) / 1.5);
    end
    qtriplot(GEOM.VER,GEOM.ITRI,'v1.5');
    qtriplot(reg,'v1.5');
    qtriplot('panel v1.5= 2 1');
end
 end