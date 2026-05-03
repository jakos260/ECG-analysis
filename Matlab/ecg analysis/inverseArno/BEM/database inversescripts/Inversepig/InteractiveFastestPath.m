if ~exist('debug','var') || ~debug
    clearvars
    close all
    % geomdir = fullfile('C:\Users\Damp2\Documents\ECG_simulation\STW\Data\geometries\');
    geomdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/';
    
    casedir		='';
    dirout = fullfile('./results/');
    lpass=5;   % # samples in lowpassma used to filter reults
    clusterDist=30; % used in multifocisacn
    
    % which leads should be used in the initial esitimate
    leadset='all';
    group=[];
    sub=9;
    saveCase =1; % 2 = overschrijven
    skipexisting=0; % skip if result files already present in target folder ( recover from crash)
    sinkScan= 0;
    layfile='pig64.mla';
    type='ventricles';
    clusters = 1;
    wrd=0; % use weighed rd in inverse
    useScaling =1;
    if wrd
        mudep=2.5e-6;
        murep=2.5e-6;
    else
        mudep=1e-5;
        murep=1e-5;
    end
    %
    %     initialvelocity=0.4;
    %     maxvelocity=1.2;
    % %     velo = 1.0;
    %     anis = 0.01;
    
    
    %initialtype
    
    initmode=1; % methods for intial estimate
    % 0: use normal fociscan
    % 1: test with fixed initial estimate focus at stimulation site
    % n ('number'): test nearest n points at same surface
    % negative: test with focus on oposite surface
    
    if sub==1
        subject='v3';
        %     bsmfile=fullfile('C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File01_beat7.ecg');beat='01';
        bsmfile='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\sinus_closedThorax_BSMonly_beat3.mat'; beat='09';
        dirout= ['.\results\' subject '\' ];
    elseif sub == 8
        subject='Pig08';
        bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig08\ECG\beats\';
        dirout= ['.\results\' subject '\'];
    elseif sub == 9
        %     subject='Pig09_refined';
        subject='Pig09';
        %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
        bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
        dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
        %     bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\avobeats\';
        %     dirout= ['.\results\' subject '\avobeats\'];
    else
        return
    end
    
    %debug
    PathName ='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats3/';
    %     FileName={    'pig09_005_155_184_RVApexPostEpiThrx1SyncVoff_20130211T133347_beat000.selecg'
    %         'pig09_005_283_312_RVApexPostEndoThrx1SyncVoff_20130207T090348_beat000.selecg'};
    %
    %     FileName={'pig09_008_75_104_LVApexEpiThrx1AsyncVoff_20130613T194832_beat000.selecg'
    %         'pig09_008_231_260_LVApexEndoThrx1_AsyncVoff_20130613T160301_beat000.selecg'};
    %
    %         FileName={'pig09_007_303_332_LVLatEndoThrx1SyncVoff_20130211T141135_beat000.selecg'
    %             'pig09_007_580_609_LVLatEpiThrx1SyncVoff_20130211T143819_beat000.selecg'};
    
    
    if ~exist('FileName') || isempty(FileName)
        % bsmfiles = dir([bsmdir '*.selecg']);
        [FileName,PathName]=uigetfile('*.selecg','Select endo and epi files for one position',bsmdir,'MultiSelect','on');
    end
    if ~iscell(FileName)
        FileName={FileName};
    end
    
    RESULTS =[];
    
    
    diroutOrg = dirout;
    anis=1; % dummy value
    for ifile = 1 : length(FileName)
        bsmfile = fullfile(bsmdir, FileName{ifile});
        beat = FileName{ifile}(1:end-7);
        orgname = beat(1:strfind(beat,'beat')-2);
        if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file')
            load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'])
            remove = DATA.remove(1:64);
            clear DATA;
        else
            error('Export file not available');
        end
        GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',anis,'group',group,'basedir',geomdir);
        GEOM.beat=beat;
        useLeads=1:size(GEOM.LAY,1) - 1;
        useLeads(remove ==1)=[];
        L= GEOM.LAY(2:end,:);
        for i=1:size(L,1)
            L(i,3) = L(i,3) - sum(remove (1:i));
        end
        L(remove ==1,:)=[];
        GEOM.LAY = [GEOM.LAY(1,:); L];
        GEOM = selectLeads(GEOM,useLeads,1);
        
        % effect of sternum
        %         A=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/pig09_ventricles_edl.mat');
        A=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/pig09_ventricles_NoSternum_edl.mat');
        GEOM.AMAnosternum=zeromean( 40*A(useLeads,:));
        A=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/pig09_ventricles_Sternum_edl.mat');
        GEOM.AMAsternum=zeromean( 40*A(useLeads,:));
        A=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/pig09_ventricles_nolungs_edl.mat');
        GEOM.AMAnolungs=zeromean( 40*A(useLeads,:));
        
        GEOM.SPECS = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',1,'dovstim',1);
        
        PSIREF = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),5));
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
        GEOM.initialActTime=initialActTime;
        
        
        surfnames={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % Simplified
        [cathpos,cathsurf]=GetCatheterPosition(sub,FileName{ifile});
        [pnearest, distp, trinearest, focus(ifile), distver]=findnearest(cathpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,initmode);
        oppsurf=unique(surfnames(~strcmp(cathsurf,surfnames)));
        [pnearest, distp, trinearest, focusopp(ifile), distver]=findnearest(cathpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,initmode);
        GEOMA(ifile)=GEOM;
    end
end % debug skip
% aniss=[0.01,0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0];
aniss=[0.1,1,2.0];
%%
close all
clear DIST2W
for ianis=1:length(aniss)
    anisdistpath=fullfile(geomdir , subject,  [subject num2str(aniss(ianis)) '_' type ]);
    if aniss(ianis)==1
        ADJ2W(ianis,:,:)=GEOMA(1).ADJ;
        DIST2W(ianis,:,:)=GEOMA(1).DIST;
    else
        if exist([anisdistpath '.adj2w'],'file')
            ADJ2W(ianis,:,:)=loadmat([anisdistpath '.adj2w']);
            if length(ADJ2W(ianis,:,:)) ~= length(GEOMA(1).ADJ)
                ADJ2W(ianis,:,:)=calcAnisADJ(GEOMA(1),aniss(ianis));
                savemat([anisdistpath '.adj2w'],squeeze(ADJ2W(ianis,:,:)));
            end
        else
            ADJ2W(ianis,:,:)=calcAnisADJ(GEOMA(1),aniss(ianis));
            savemat([anisdistpath '.adj2w'],squeeze(ADJ2W(ianis,:,:)));
        end
        if exist([anisdistpath '.dist'],'file')
            DIST2W(ianis,:,:)=loadmat([anisdistpath '.dist']);
            if length(DIST2W(ianis,:,:)) ~= length(GEOMA(1).ADJ)
                DIST2W(ianis,:,:)=graphdist(squeeze(ADJ2W(ianis,:,:)));
                savemat([anisdistpath '.dist'],squeeze(DIST2W(ianis,:,:)));
            end
        else
            DIST2W(ianis,:,:)=graphdist(squeeze(ADJ2W(ianis,:,:)));
            savemat([anisdistpath '.dist'],squeeze(DIST2W(ianis,:,:)));
        end
    end
    
    % aniss(ianis)
    % temp=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',aniss(ianis),'group',group,'basedir',geomdir);
    % DIST2W(ianis,:,:)=temp.DIST2W;
    % size(DIST2W)
    % qtriplot(squeeze(DIST2W(ianis,:,:))');
end

anis=0.1;
maxVelocity=1.5;
initialVelocity=0.4;
ianis=find(aniss==anis);
ianisinit=find(aniss==1.0);
qrsstart=false;
funscal=2;
sigwidth=110;
usesternum=0;
initialActTime=27;
opposite=false;
depdelay=0;
nohakkel=true;
smode=1;
blflag=1;
blstring={'none','onsetQRS','start beat','Before Vstim'};


fh=figure(34);
clf
set(fh,'WindowStyle','modal');
updateqtriplot=false;
calcfw=true;
loadanis=true;
stopflag=false;
selgeom=0; % 0 all signals, >0 GEOMA(selgeom) with opposite
k=' '; % notempty so plot at first iteration
while ~stopflag
    if isempty(k)
        waitforbuttonpress;
        k=get(gcf,'currentcharacter');
        %     k=getkey;
        switch k
            case ''
                continue;
            case ' '
                stopflag=true;
                continue;
            case 'a'
                ianis=max(1,ianis-1);
                loadanis=true;
                calcfw=true;
            case 'A'
                ianis=min(length(aniss),ianis+1);
                loadanis=true;
                calcfw=true;
            case 'g'
                selgeom=max(0,selgeom-1);
                calcfw=true;
            case 'G'
                selgeom=min(length(GEOMA),selgeom+1);
                calcfw=true;
            case 'q'
                ianisinit=max(1,ianisinit-1);
                loadanis=true;
                calcfw=true;
            case 'Q'
                ianisinit=min(length(aniss),ianisinit+1);
                loadanis=true;
                calcfw=true;
            case 'v'
                initialVelocity=max(0,initialVelocity-0.1);
                calcfw=true;
                %veloi=max(1,veloi-1);
            case 'V'
                initialVelocity=min(2.0,initialVelocity+0.1);
                calcfw=true;
                %             veloi=min(length(FWendo.velos),veloi+1);
            case 'f'
                maxVelocity=max(0.3,maxVelocity-0.1);
                calcfw=true;
                %veloi=max(1,veloi-1);
            case 'F'
                maxVelocity=min(2.0,maxVelocity+0.1);
                calcfw=true;
                %             veloi=min(length(FWendo.velos),veloi+1);
            case's'
                funscal=funscal/2;
            case 'S'
                funscal=funscal*2;
            case'd'
                depdelay=max(-16,depdelay-4);
                calcfw=true;
                
            case 'D'
                depdelay=min(60,depdelay+4);
                calcfw=true;
                
            case 't'
                sigwidth=max(10,sigwidth-5);
            case 'T'
                sigwidth=min(sigwidth+5,300);
            case 'p'
                qrsstart=~qrsstart;
                calcfw=true;
            case 'P'
                %             qrsstart=false;
                %             calcfw=true;
                pause
            case 'i'
                initialActTime=max(-5, initialActTime-5);
                calcfw=true;
            case 'I'
                initialActTime=min(60, initialActTime+5);
                calcfw=true;
            case 'm'
                %                 usesternum=max(0, usesternum-1);
                smode=1;
                calcfw=true;
            case 'M'
                smode=4;
                calcfw=true;
                %                 usesternum=min(2, usesternum+1);
            case 'u'
                updateqtriplot=1;
            case 'U'
                updateqtriplot=2
            case {'o','O'}
                opposite=~opposite;
                calcfw=true;
            case {'h','H'}
                nohakkel=~nohakkel;
                calcfw=true;
            case 'b'
                blflag=max(1,blflag-1);
                calcfw=true;
            case 'B'
                blflag=min(length(blstring),blflag+1);
                calcfw=true;
        end
    end
    k='';
    anis=aniss(ianis);
    %     velo=FWendo.velos(veloi);
    
    if qrsstart
        syncstr='Onset QRS';
    else
        syncstr='Pace';
    end
    
    
    
    if loadanis
        for igeom=1:length(GEOMA)
            GEOMA(igeom).DIST2W2=DIST2W(ianis,:,:);
            GEOMA(igeom).DIST2W1=DIST2W(ianisinit,:,:);
        end
        loadanis=false;
    end
    
    
    if selgeom==0
        igeomsel=1:length(GEOMA);
    else
        igeomsel=selgeom;
    end
    
    fh=figure(34);
    clf
    hold on
    if selgeom==0
        set(fh,'Name',sprintf('BSM Endo(red) and Epi(blue) and forward signals (-.) starting at %s for %dms.Delay: %d. INIT: Velo%.01f, Anis %0.2f for 1st %dms.  Anis:%0.2f. maxVelo:%0.1f. NoHakkel:%d Smode:%d BL:%s',...
            syncstr,sigwidth,depdelay,initialVelocity,aniss(ianisinit),initialActTime, anis,maxVelocity,nohakkel,smode,blstring{blflag}));
    else
        set(fh,'Name',sprintf('BSM Endo(red) and Epi(blue) and forward signals (-.) starting at %s for %dms.Delay %d. INIT: Velo%.01f, Anis %0.2f for 1st %dms.  Anis:%0.2f. maxVelo:%0.1f. NoHakkel:%d Smode:%d BL:%s',...
            syncstr,sigwidth,depdelay,initialVelocity,aniss(ianisinit),initialActTime, anis,maxVelocity,nohakkel,smode,blstring{blflag}));
    end
    
    
    for igeom=igeomsel
        if qrsstart
            qrsduration=GEOMA(igeom).SPECS.time_Jpoint;
            sigstart=GEOMA(igeom).SPECS.onsetqrs;
            initActTimeCor=0;
        else
            qrsduration=GEOMA(igeom).SPECS.time_Jpoint+GEOMA(igeom).SPECS.onsetqrs-GEOMA(igeom).SPECS.time_Vstim;
            sigstart=GEOMA(igeom).SPECS.time_Vstim;
            initActTimeCor=GEOMA(igeom).SPECS.onsetqrs-GEOMA(igeom).SPECS.time_Vstim;
        end
        if calcfw
            % calc signals
            %             if igeom==0
            %                 foc=focusopp(selgeom);
            %             else
            foco=focusopp(igeom);
            foc=focus(igeom);
            %             end
            if initialActTime<0
                usedInitialActTime=GEOMA(igeom).initialActTime+initActTimeCor;
                sprintf('InitialActTime(used): %d(%d)\n',GEOMA(igeom).initialActTime,usedInitialActTime);
            else
                usedInitialActTime=initialActTime;%-depdelay;
            end

            if usedInitialActTime>depdelay && initialVelocity>0
                dep=DIST2W(ianisinit,:,foc)/initialVelocity+depdelay;
                %                 depw=DIST2W(ianis,:,foc)/maxVelocity+depdelay;
                depw=DIST2W(ianis,:,foc)/initialVelocity+depdelay;
                
                
                %         dep=min(GEOM.DIST2W(:,1791),GEOM.DIST2W(:,1556))/initialVelocity; % Endo midden in de wand
                depAbove=find(dep>usedInitialActTime);
                %             dep(depAbove) = initialActTime + (dep(depAbove) - initialActTime)*initialVelocity / velo;
                %             dt0=min(20,initialActTime);
                dt0=usedInitialActTime;
                dep(depAbove)=depdelay+dt0 + (depw(depAbove) - dt0) * initialVelocity / min(maxVelocity,(max(depw) - dt0) * initialVelocity/ (qrsduration - dt0-depdelay) );
                %             dep(depAbove)=initialActTime + (dep(depAbove)-initialActTime) * initialVelocity / min(maxVelocity,(max(dep) - initialActTime) * initialVelocity/ (qrsduration - initialActTime) );
                %             DEP(igeom,:)=dep;
                if selgeom>0
                    depo=DIST2W(ianisinit,:,foco)/initialVelocity+depdelay;
                    depw=DIST2W(ianis,:,foco)/initialVelocity+depdelay;
                    
                    %         dep=min(GEOM.DIST2W(:,1791),GEOM.DIST2W(:,1556))/initialVelocity; % Endo midden in de wand
                    depAbove=find(depo>usedInitialActTime);
                    %             dep(depAbove) = initialActTime + (dep(depAbove) - initialActTime)*initialVelocity / velo;
                    %             dt0=min(20,initialActTime);
                    dt0=usedInitialActTime;
                    %                     fixen!!!
                    depo(depAbove)=depdelay + dt0 + (depw(depAbove) - dt0) * initialVelocity / min(maxVelocity,(max(depw) - dt0) * initialVelocity/ (qrsduration - dt0-depdelay) );
                    %             dep(depAbove)=initialActTime + (dep(depAbove)-initialActTime) * initialVelocity / min(maxVelocity,(max(dep) - initialActTime) * initialVelocity/ (qrsduration - initialActTime) );
                    %             DEP(igeom,:)=dep;
                end
            else
                dep1=DIST2W(ianis,:,foc); % if initialActTime~=0, apparently initVelocity==0
                dep=depdelay + dep1/min(maxVelocity,max(dep1)*(qrsduration-depdelay));
                if selgeom>0 
                    
                    dep1=DIST2W(ianis,:,foco); % if initialActTime~=0, apparently initVelocity==0
                    depo=depdelay + dep1/min(maxVelocity,max(dep1)*(qrsduration-depdelay));
                end
                
                
            end
            DEP(igeom,:)=dep;
            if smode==1
                rep        = 300 * ones(size(GEOMA(igeom).VER,1),1) + GEOMA(igeom).SPECS.time_apexT;
            else
                rho=zeros(size(dep));
                for i=1:length(dep)
                    a=find(GEOMA(igeom).DIST(i,:) < 50 & GEOMA(igeom).DIST(i,:) > 0);
                    rho(i) = sum((dep(a) - dep(i)) ./ (GEOMA(igeom).DIST(i,a).^2));
                end
                %qrsduration = GEOM.SPECS.qrsduration;
                ari =  GEOMA(igeom).SPECS.time_apexT - mean(dep) + rho*7 + GEOMA(igeom).SPECS.repCorrection;
                
                rep =ari + dep;
            end
            %             switch usesternum
            %                 case 0
            
            if nohakkel
                FWSIG{igeom} =GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode,GEOMA(igeom));
                if selgeom>0, fwsigo =GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),depo,rep, GEOMA(igeom).SPECS,smode,GEOMA(igeom));  end
                %                 case 1
                fwnosternum=GEOMA(igeom).AMAnosternum*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode,GEOMA(igeom));
                %                 case 2
                fwsternum =GEOMA(igeom).AMAsternum*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode,GEOMA(igeom));
                fwnolungs=GEOMA(igeom).AMAnolungs*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode,GEOMA(igeom));
                
            else
                FWSIG{igeom} =lowpassma(GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode),0);
                if selgeom>0, fwsigo =lowpassma(GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),depo,rep, GEOMA(igeom).SPECS,smode),0);  end
                %                 case 1
                fwnosternum=lowpassma(GEOMA(igeom).AMAnosternum*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode),0);
                %                 case 2
                fwsternum =lowpassma(GEOMA(igeom).AMAsternum*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode),0);
                fwnolungs=lowpassma(GEOMA(igeom).AMAnolungs*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:400),dep,rep, GEOMA(igeom).SPECS,smode),0);
                %             end
                %             GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:sigwidth),dep,300*ones(size(GEOMA(igeom).VER,1),1),pS,[],1);
                
            end
            % Calc wrd, rd, corr QRS/QRST, same on 9 opposite leads
            tfw_start=GEOMA(igeom).SPECS.onsetqrs-GEOMA(igeom).SPECS.time_Vstim+1; % onset QRS assuming activation started at stim
            tfw_end=tfw_start+GEOMA(igeom).SPECS.time_Jpoint;
            tfwfull=tfw_start:tfw_end;
            tfwhalf=tfw_start:tfw_start+ceil((tfw_end-tfw_start)/2);
            tbsm_start=GEOMA(igeom).SPECS.onsetqrs;
            tbsm_end=tbsm_start+GEOMA(igeom).SPECS.time_Jpoint;
            tbsmfull=tbsm_start:tbsm_end;
            tbsmhalf=tbsm_start:tbsm_start+ceil((tbsm_end-tbsm_start)/2);
            
            
            
            %             RDfull(igeom)=norm(FWSIG{igeom}(:,tfwfull)-GEOMA(igeom).BSM(:,tbsmfull),'fro')/norm(GEOMA(igeom).BSM(:,tbsmfull),'fro');
            %             RDhalf(igeom)=norm(FWSIG{igeom}(:,tfwhalf)-GEOMA(igeom).BSM(:,tbsmhalf),'fro')/norm(GEOMA(igeom).BSM(:,tbsmhalf),'fro');
            %             RDWfull(igeom)=sum(rms(FWSIG{igeom}(:,tfwfull)-GEOMA(igeom).BSM(:,tbsmfull))./(0.001+rms(GEOMA(igeom).BSM(:,tbsmfull))))/length(tfwfull);
            %             RDWhalf(igeom)=sum(rms(FWSIG{igeom}(:,tfwhalf)-GEOMA(igeom).BSM(:,tbsmhalf))./(0.001+rms(GEOMA(igeom).BSM(:,tbsmhalf))))/length(tfwfull);
            %             CORfull(igeom)=compCor(FWSIG{igeom}(:,tfwfull),GEOMA(igeom).BSM(:,tbsmfull),1);
            %             CORhalf(igeom)=compCor(FWSIG{igeom}(:,tfwhalf),GEOMA(igeom).BSM(:,tbsmhalf),1);
            RDfull(igeom)=norm(FWSIG{igeom}(:,tfwfull)-GEOMA(igeom).BSM(:,tbsmfull),'fro')/norm(GEOMA(igeom).BSM(:,tbsmfull),'fro');
            RDhalf(igeom)=norm(FWSIG{igeom}(:,tfwhalf)-GEOMA(igeom).BSM(:,tbsmhalf),'fro')/norm(GEOMA(igeom).BSM(:,tbsmhalf),'fro');
            RDWfull(igeom)=sum(rms(FWSIG{igeom}(:,tfwfull)-GEOMA(igeom).BSM(:,tbsmfull))./(0.001+rms(GEOMA(igeom).BSM(:,tbsmfull))))/length(tfwfull);
            RDWhalf(igeom)=sum(rms(FWSIG{igeom}(:,tfwhalf)-GEOMA(igeom).BSM(:,tbsmhalf))./(0.001+rms(GEOMA(igeom).BSM(:,tbsmhalf))))/length(tfwfull);
            CORfull(igeom)=compCor(FWSIG{igeom}(:,tfwfull),GEOMA(igeom).BSM(:,tbsmfull),1);
            CORhalf(igeom)=compCor(FWSIG{igeom}(:,tfwhalf),GEOMA(igeom).BSM(:,tbsmhalf),1);
            fprintf('%s:; RDfulld(half):%0.3f(%0.3f), RDWfull(half):%0.3f(%0.3f), CORful(half):%0.3f(%0.3f)\n',...
                cell2mat(regexpi(FileName{igeom},'(ENDO|EPI)','match')),RDfull(igeom),RDhalf(igeom),RDWfull(igeom),RDWhalf(igeom),CORfull(igeom),CORhalf(igeom));
            
            %             TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
            %             TST.wrd = sum(rms(TST.RES) ./ (0.0010 + rms(INV.PHIREF))); % weighted rd
            %                 rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
            %     cors(inode) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
            % %     wrds(inode) = sum(rms(SCAN.PSIREF(:,1:maxt) - PSIA) ./ (0.0010 + rms(SCAN.PSIREF(:,1:maxt)))); % weighted rd
            
            
            
        end
        
        % d=0;% 10 pace sync, 0 qrs sync
        if  ~isempty(strfind(GEOMA(igeom).beat,'Endo'))
            linespec='r';
            linespeco='b';
        else
            linespec='b';
            linespeco='r';
        end
        
        
        %         figure(21);
        %         if igeom==1, clf; end
        %         sigplot(FWSIG{igeom}(:,tfwfull),'',GEOMA(igeom).LAY,funscal,['-.' linespec]);
        %         hold on
        %         sigplot(GEOMA(igeom).BSM(:,tbsmfull),'',GEOMA(igeom).LAY,funscal,linespec);
        %         hold on
        %         figure(34)
        
         switch blflag
            case 1
                bsm=GEOMA(igeom).BSM(:,sigstart:sigstart+sigwidth-1);
             case 2
                 bsmt=baselinecor(GEOMA(igeom).BSM,GEOMA(igeom).SPECS.onsetqrs,GEOMA(igeom).SPECS.endtwave);
                 bsm=bsmt(:,sigstart:sigstart+sigwidth-1);
             case 3
                 bsmt=baselinecor(GEOMA(igeom).BSM,1,GEOMA(igeom).SPECS.endtwave);
                 bsm=bsmt(:,sigstart:sigstart+sigwidth-1);
             case 4
                 bsmt=baselinecor(GEOMA(igeom).BSM,GEOMA(igeom).SPECS.time_Vstim-10,GEOMA(igeom).SPECS.endtwave);
                 bsm=bsmt(:,sigstart:sigstart+sigwidth-1);
                 
         end
         
        
        sigplot(max(-2/funscal,min(2/funscal,bsm)),'',GEOMA(igeom).LAY,funscal,linespec);
        hold on
       
        
        plot(0:1/(sigwidth-1):1,rms(GEOMA(igeom).BSM(:,sigstart:sigstart+sigwidth-1)),linespec);
        %     sigplot(ENDO.DATA.AVERAGEPWC(:,d+143:d+143+sigwidth),'old',ENDO.mlas,funscal,'red');
        
        
        % figure
        % d2=0; % 0 pace sync, -10 qrs sync
        sigplot(FWSIG{igeom}(:,1:sigwidth),'',GEOMA(igeom).LAY,funscal,['-.' linespec]);
        hold on
        plot(0:1/(sigwidth-1):1,rms(FWSIG{igeom}(:,1:sigwidth)),['-.' linespec]);
        if selgeom>0
            sigplot(fwsigo(:,1:sigwidth),'',GEOMA(igeom).LAY,funscal,['-.' linespeco]);
            hold on
        end
        
        %         sigplot(fwnosternum(:,1:sigwidth),'',GEOMA(igeom).LAY,funscal,['-' linespec]);
        %         hold on
        %                 sigplot(fwnolungs(:,1:sigwidth),'',GEOMA(igeom).LAY,funscal,['-.' linespec]);
        %         hold on
        
        
        
        sigplot(0.1*ones(size(GEOMA(igeom).BSM,1),sigwidth),'',GEOMA(igeom).LAY,funscal,':k');
        hold on
        sigplot(-0.1*ones(size(GEOMA(igeom).BSM,1),sigwidth),FileName{igeom},GEOMA(igeom).LAY,funscal,':k');
        hold on
        
        %         hold on
        %         sigplot(squeeze(FWendo.fwsig(anisi,veloi,:,1:sigwidth)),'epi',GEOMA(igeom).LAY,funscal,':r');
        
        
    end
    if updateqtriplot
        if selgeom==0
            for igeom=1:length(GEOMA)
                if igeom==1
                    qtriplot(sprintf('horizontal %d',length(GEOMA)));
                    qtriplot('vertical 1');
                    qtriplot('delete *');
                    qtriplot('funcolor tim');
                    qtriplot('funscale 1 110');
                    qtriplot('step 10');
                end
                qtriplot(GEOMA(igeom).VER,GEOMA(igeom).ITRI,FileName{igeom});
                qtriplot(DEP(igeom,:)');
                %             qtriplot(CalcTransmuralDiff(GEOMA(igeom).VER, GEOMA(igeom).ITRI, GEOMA(igeom).typ, DEP(igeom,:),GEOMA(igeom).DIST)');
                qtriplot(['marker vertex ',num2str(focus(igeom))]);
                qtriplot(['panel ' num2str(igeom) ',1']);
            end
        elseif updateqtriplot==1
            qtriplot('horizontal 2');
            qtriplot('vertical 1');
            qtriplot('delete *');
            qtriplot(GEOMA(igeomsel).VER,GEOMA(igeomsel).ITRI,'Diff');
            qtriplot('fun Diff=1');
            qtriplot('funcolor tim');
            qtriplot('funscale auto');
            qtriplot('step 2');
            %             qtriplot(DEP(igeom,:)');
            qtriplot(CalcTransmuralDiff(GEOMA(igeom).VER, GEOMA(igeom).ITRI, GEOMA(igeom).typ, DEP(igeom,:),GEOMA(igeom).DIST)');
            qtriplot(['marker vertex ',num2str(focus(igeom))]);
            qtriplot('panel 2 1');
            qtriplot(GEOMA(igeomsel).VER,GEOMA(igeomsel).ITRI,'Dep');
            qtriplot('fun Dep=2');
            qtriplot('funcolor2 tim');
            qtriplot('funscale2 1 110');
            qtriplot('step2 10');
            
            qtriplot(DEP(igeom,:)');
            %             qtriplot('funcolor tim');
            %             qtriplot('funscale 1 110');
            
            qtriplot(['marker vertex ',num2str(focus(igeom))]);
            qtriplot('panel 1 1');
        else
%             sliceheart(GEOMA(igeomsel).VER,GEOMA(igeomsel).ITRI,GEOMA(geomsel).type,fun,markers,funscalemax,rot
        end
        
        updateqtriplot=false;
    end
    calcfw=false;
end
set(fh,'WindowStyle','normal');
