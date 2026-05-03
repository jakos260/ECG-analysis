clearvars
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
    wrd=1; % use weighed rd in inverse
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
        bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats3/';
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
    %     FileName={'pig09_007_303_332_LVLatEndoThrx1SyncVoff_20130211T141135_beat000.selecg'
    %         'pig09_007_580_609_LVLatEpiThrx1SyncVoff_20130211T143819_beat000.selecg'};
    
    
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
        
        
        GEOM.SPECS = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',1);
        
        surfnames={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % Simplified
        [cathpos,cathsurf]=GetCatheterPosition(sub,FileName{ifile});
        [pnearest, distp, trinearest, focus(ifile), distver]=findnearest(cathpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,initmode);
        oppsurf=unique(surfnames(~strcmp(cathsurf,surfnames)));
        [pnearest, distp, trinearest, focusopp(ifile), distver]=findnearest(cathpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,-initmode);
        GEOMA(ifile)=GEOM;
    end
end % debug skip
aniss=[0.01,0.1,0.3,0.5,1.0,1.5,2.0,2.5];
velocities=[0.1,0.3,0.5,0.7,0.8,1.0,1.2,1.4,1.6];
inittimes=[0,10,20,30,40,50];
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

% anis=0.01
% maxVelocity=1.5;
% initialVelocity=0.4;
% ianis=find(aniss==anis);
% ianisinit=find(aniss==1.0);
qrsstart=false;
% pacesync=true;
funscal=2;
sigwidth=110;
usesternum=0;
initialActTime=20;
opposite=false;

doplot=false;
doplot2=false; % plot signals for RD etc calculation
updateqtriplot=false;


fh=figure(34);
clf
% set(fh,'WindowStyle','modal');


display('ianisinit vast op 1');
for ianisinit=1%:length(aniss)
    anisinit=aniss(ianisinit);
    for iinitialVelocity=1:length(velocities)
        initialVelocity=velocities(iinitialVelocity);
        for iinitialActTime=1:length(inittimes)
            initialActTime=inittimes(iinitialActTime);
            %             for imaxVelocity=1:length(velocities)
            maxVelocity=2.0;%velocities(imaxvelocity);
            for ianis=1:length(aniss)
                anis=aniss(ianis);
                
                for igeom=1:length(GEOMA)
                    GEOMA(igeom).DIST2W=DIST2W(ianis,:,:);
                end
                for iopposite=1:2
                    opposite=iopposite-1;
                    if qrsstart
                        syncstr='Onset QRS';
                    else
                        syncstr='Pace';
                    end
                    if doplot
                        
                        fh=figure(34);
                        clf
                        hold on
                        set(fh,'Name',sprintf('Endo(red) and Epi(blue) Opposite:%d, from %s for %dms.INIT: Velo%.01f, Anis %0.2f for 1st %dms.  Anis:%0.2f. maxVelo:%0.1f.',...
                            opposite,syncstr,sigwidth,initialVelocity,aniss(ianisinit),initialActTime, anis,maxVelocity));
                        
                    end
                    
                    
                    for igeom=1:length(GEOMA)
                        if qrsstart
                            qrsduration=GEOMA(igeom).SPECS.time_Jpoint;
                            sigstart=GEOMA(igeom).SPECS.onsetqrs;
                        else
                            qrsduration=GEOMA(igeom).SPECS.time_Jpoint+GEOMA(igeom).SPECS.onsetqrs-GEOMA(igeom).SPECS.time_Vstim;
                            sigstart=GEOMA(igeom).SPECS.time_Vstim;
                        end
                        % calc signals
                        if opposite
                            foc=focusopp(igeom);
                        else
                            foc=focus(igeom);
                        end
                        if initialActTime>0 && initialVelocity>0
                            dep=DIST2W(ianisinit,:,foc)/initialVelocity;
                            depw=DIST2W(ianis,:,foc)/initialVelocity;
                            
                            %         dep=min(GEOM.DIST2W(:,1791),GEOM.DIST2W(:,1556))/initialVelocity; % Endo midden in de wand
                            depAbove=find(dep>initialActTime);
                            %             dep(depAbove) = initialActTime + (dep(depAbove) - initialActTime)*initialVelocity / velo;
                            %             dt0=min(20,initialActTime);
                            dt0=initialActTime;
                            dep(depAbove)=dt0 + (depw(depAbove) - dt0) * initialVelocity / min(maxVelocity,(max(depw) - dt0) * initialVelocity/ (qrsduration - dt0) );
                            %             dep(depAbove)=initialActTime + (dep(depAbove)-initialActTime) * initialVelocity / min(maxVelocity,(max(dep) - initialActTime) * initialVelocity/ (qrsduration - initialActTime) );
                            %             DEP(igeom,:)=dep;
                               
                        else
                            dep1=DIST2W(ianis,:,foc); % if initialActTime~=0, apparently initVelocity==0
                            dep=initialActTime + dep1/min(maxVelocity,max(dep1)*(qrsduration-initialActTime));
                            
                            
                            
                        end
                        
                        if updateqtriplot
                            DEP(igeom,:)=dep;
                        end
                        
                        rep        = 300 * ones(size(GEOMA(igeom).VER,1),1) + GEOMA(igeom).SPECS.time_apexT;
                        %             switch usesternum
                        %                 case 0
                        FWSIG{igeom} =lowpassma(GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:200),dep,rep, GEOMA(igeom).SPECS,1,GEOMA(igeom)),5);
                        %                 case 1
                        %             fwnosternum=lowpassma(GEOMA(igeom).AMAnosternum*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:200),dep,rep, GEOMA(igeom).SPECS,1,GEOMA(igeom)),5);
                        %             %                 case 2
                        %             fwsternum =lowpassma(GEOMA(igeom).AMAsternum*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:200),dep,rep, GEOMA(igeom).SPECS,1,GEOMA(igeom)),5);
                        %             fwnolungs=lowpassma(GEOMA(igeom).AMAnolungs*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:200),dep,rep, GEOMA(igeom).SPECS,1,GEOMA(igeom)),5);
                        %             end
                        %             GEOMA(igeom).AMA*getSmode(ones(length(GEOMA(igeom).VER),1) * (1:sigwidth),dep,300*ones(size(GEOMA(igeom).VER,1),1),pS,[],1);
                        
                        
                        % Calc wrd, rd, corr QRS/QRST, same on 9 opposite leads
                        if qrsstart
                            tfw_start=1;
                        else
                            tfw_start=GEOMA(igeom).SPECS.onsetqrs-GEOMA(igeom).SPECS.time_Vstim+1;
                        end
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
                        
                        
                        
                        
                        RDfull(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)=norm(FWSIG{igeom}(:,tfwfull)-GEOMA(igeom).BSM(:,tbsmfull),'fro')/norm(GEOMA(igeom).BSM(:,tbsmfull),'fro');
                        RDhalf(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)=norm(FWSIG{igeom}(:,tfwhalf)-GEOMA(igeom).BSM(:,tbsmhalf),'fro')/norm(GEOMA(igeom).BSM(:,tbsmhalf),'fro');
                        RDWfull(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)=sum(rms(FWSIG{igeom}(:,tfwfull)-GEOMA(igeom).BSM(:,tbsmfull))./(0.001+rms(GEOMA(igeom).BSM(:,tbsmfull))))/length(tfwfull);
                        RDWhalf(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)=sum(rms(FWSIG{igeom}(:,tfwhalf)-GEOMA(igeom).BSM(:,tbsmhalf))./(0.001+rms(GEOMA(igeom).BSM(:,tbsmhalf))))/length(tfwfull);
                        CORfull(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)=compCor(FWSIG{igeom}(:,tfwfull),GEOMA(igeom).BSM(:,tbsmfull),1);
                        CORhalf(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)=compCor(FWSIG{igeom}(:,tfwhalf),GEOMA(igeom).BSM(:,tbsmhalf),1);
                     
                        if doplot2 && igeom==1 && iopposite==1
                            fh=figure(35);
                            clf
                            hold on
                            set(fh,'Name',sprintf('Endo(red) and Epi(blue) Opposite:%d, from %s for %dms.INIT: Velo%.01f, Anis %0.2f for 1st %dms.  Anis:%0.2f. maxVelo:%0.1f.RDWhalf%0.1f',...
                                opposite,syncstr,sigwidth,initialVelocity,aniss(ianisinit),initialActTime, anis,maxVelocity,...
                                RDWhalf(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom) ));
                            
                            sigplot(GEOMA(igeom).BSM(:,tbsmfull),'',GEOMA(igeom).LAY,funscal,'b');
                            hold on
                            sigplot(FWSIG{igeom}(:,tfwfull),'',GEOMA(igeom).LAY,funscal,'r');
                            hold on
                            fh2=figure(36);
                            set(fh2,'Name','RDWhalf');
                            sz=[length(aniss),length(velocities),length(inittimes),length(aniss),2,1];
                            index=sub2ind(sz,ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,1);
                            plot(index,RDWhalf(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom),'xr');
                            hold on
                            sig=[ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite];
                            %                             plot(repmat(index,1,length(sig)),sig,'+');
                            pause
                            plot(index,RDWhalf(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom),'xk');
                            
                        end
                        
                        
                        
                        
                        
                        %             fprintf('%s:; RDfulld(half):%0.3f(%0.3f), RDWfull(half):%0.3f(%0.3f), CORful(half):%0.3f(%0.3f)\n',...
                        %                 cell2mat(regexpi(FileName{igeom},'(ENDO|EPI)','match')),RDfull(igeom),RDhalf(igeom),RDWfull(igeom),RDWhalf(igeom),CORfull(igeom),CORhalf(igeom));
                        
                        %             TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
                        %             TST.wrd = sum(rms(TST.RES) ./ (0.0010 + rms(INV.PHIREF))); % weighted rd
                        %                 rds(inode) = norm(SCAN.PSIREF(:,1:maxt) - PSIA,'fro')/norm(SCAN.PSIREF(:,1:maxt),'fro');
                        %     cors(inode) =  compCor(SCAN.PSIREF(:,1:size(PSIA,2)),PSIA,corMode);
                        % %     wrds(inode) = sum(rms(SCAN.PSIREF(:,1:maxt) - PSIA) ./ (0.0010 + rms(SCAN.PSIREF(:,1:maxt)))); % weighted rd
                        
                        
                        
                        if doplot
                            
                            % d=0;% 10 pace sync, 0 qrs sync
                            if ~isempty(strfind(GEOMA(igeom).beat,'Endo'))
                                linespec='r';
                            else
                                linespec='b';
                            end
                            
                            
                            %         figure(21);
                            %         if igeom==1, clf; end
                            %         sigplot(FWSIG{igeom}(:,tfwfull),'',GEOMA(igeom).LAY,funscal,['-.' linespec]);
                            %         hold on
                            %         sigplot(GEOMA(igeom).BSM(:,tbsmfull),'',GEOMA(igeom).LAY,funscal,linespec);
                            %         hold on
                            %         figure(34)
                            
                            
                            
                            sigplot(max(-2/funscal,min(2/funscal,GEOMA(igeom).BSM(:,sigstart:sigstart+sigwidth-1))),'',GEOMA(igeom).LAY,funscal,linespec);
                            hold on
                            %     sigplot(ENDO.DATA.AVERAGEPWC(:,d+143:d+143+sigwidth),'old',ENDO.mlas,funscal,'red');
                            
                            
                            % figure
                            % d2=0; % 0 pace sync, -10 qrs sync
                            sigplot(FWSIG{igeom}(:,1:sigwidth),'',GEOMA(igeom).LAY,funscal,['-.' linespec]);
                            hold on
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
                        
                    end
                    if updateqtriplot
                        for igeom=1:length(GEOMA)
                            if igeom==1
                                qtriplot(sprintf('horizontal %d',length(GEOMA)));
                                qtriplot('vertical 1');
                                qtriplot('delete *');
                            end
                            qtriplot(GEOMA(igeom).VER,GEOMA(igeom).ITRI,['Ventricle' num2str(igeom)]);
                            qtriplot(DEP(igeom,:)');
                            qtriplot(['marker vertex ',num2str(focus(igeom))]);
                            qtriplot(['panel ' num2str(igeom) ',1']);
                        end
                        %         updateqtriplot=false;
                    end
                    calcfw=false;
                end
            end
        end
    end
end
fp=fopen([mfilename '.m'],'r');
codestr=fread(fp,inf,'char');
fclose(fp);
save(['BatchFastestPathOutput' datestr(now,'yyyymmddHHMMSS') '.mat'],'RD*','COR*','aniss','velocities','inittimes','codestr','PathName','FileName','qrsstart');
                
