% function showPigsPO()
disp(['starting ' mfilename]);
clearvars
global lpass
global clusterDist
clear qtriscriptl

makescript=true;%true; % create qtriplot script. One for all, one for each line/orgname.

qtriscriptall={'delete *','horizontal 6','vertical 9',};
ncol=3; % number of columns in per line view. number of row is hardcode in the next line:
% qtriscriptlinit={'delete *',['horizontal ' num2str(ncol+1)],'vertical 2',}; %init for per line script

geomdir = '/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries';


casedir		='';
dirout='./results/';
lpass=5;   % # samples in lowpassma used to filter reults
clusterDist=30; % used in multifociscan

% which leads should be used in the initial esitimate
leadset='all';
group=[];
sub='b01';
saveCase =1;
sinkScan= 0;
layfile='pig64.mla';
type='ventricles';
thoraxname='_thoraxpigs_adam.tri';
thoraxscale=0.1; % scale thorax to display ventricles and thorax simultaneously
wctorg=[31 63 64]; % default channels for pigs
% wctorg=[63 64 65]; % default channels for ams



if sub==1
    subject='v3';
    bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File01_beat7.ecg';beat='01';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File02_beat6.ecg';beat='02';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File03_beat6.ecg';beat='03';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File04_beat6.ecg';beat='04';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File05_beat10.ecg';beat='05';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File06_beat8.ecg';beat='06';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File07_beat9.ecg';beat='07';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File08_beat9.ecg';beat='08';
    bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\sinus_closedThorax_BSMonly_beat3.mat'; beat='09';
    dirout= ['.\results\' subject '\' ];
elseif sub == 5
    subject='v5';
    bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig05\ECG\101930_sr2_10_15_SR_2_ClosedThorax_20111212T143054.bsm'; beat='sinus';
    %     bsmfile='C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig05\ECG\101930_rvpace_0_5_RV Pace, Closed Thorax_20120209T164124.bsm'; beat='rpaced';
    dirout= ['.\results\' subject '\' ];
elseif sub == 7
    subject='Pig07';
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/beats/';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/blind/beatsPvD/';
    %     dirout= ['./results/' subject '/V0_BSMCT/']; % results files. Here a
    %     sub-folder script is created
    %     dirout='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/homogeneous/';
    %     dirout='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/PO_Ratio1.0QRSScaleShiftInit10';
    dirout='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/NotInWallmudep1e-4murep1e-4';
elseif sub == 71
    geomdir = '/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/geometries/v1';
    subject='Pig07';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/beats/';
    dirout= ['./results/' subject '/v1/'];
    
elseif sub == 72
    geomdir = '/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/geometries/v2';
    subject='Pig07';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/beats/';
    dirout= ['./results/' subject '/v2/'];
    
elseif sub == 73
    geomdir = '/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/geometries/v3';
    subject='Pig07';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/beats/';
    dirout= ['./results/' subject '/v3/'];
elseif sub == 8
    subject='Pig08';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/Pacing/blind/beats/';
    dirout='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig08/single/';
elseif sub == 9
    subject='Pig09';
    %     subjectbsm='Pig09'; % alternative subject id for bsm, used with multiple geometries with same bsm data
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/beats/';
    dirout='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/';
elseif sub == 10
    subject='Pig10';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/Biosemi/export/beats/';
    dirout='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig10/single/';
elseif sub=='b01'
    subject='Bucket01';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/Biosemi/export/beats/';
    %     dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/single'];
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Bucket01/cluster4'];
    thoraxname='_thorax.tri';
    thoraxscale=1;
elseif sub == 48
    subject='ajm048';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/CARTOData/ajm048/export/beats/';
    dirout=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/single/'];
    %     dirout=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    wctorg=[63 64 35]; % 35 should be 65, but this was not exported. default channels for ams
    UsePrevAsInit=true;
else
    return
end

% [orgnames,solutionflag,electrodes]=getelectrodemarkers(sub); % get electrode positions etc.
if ~exist('subjectbsm','var')
    subjectbsm=subject;
end


dirout=fullfile(dirout);
scriptdir=fullfile(dirout,'script');
if makescript
    if ~exist(scriptdir,'dir'), mkdir(scriptdir); end
    delete(fullfile(scriptdir,'*'));
    fnameall=sprintf('ID%d',sub);
end

bsmfiles = dir([bsmdir '*.selecg']);
initdatafiles = dir([dirout '*init.mat']);
findatafiles = dir([dirout '*2.mat']);
qtriplot('delete *');
qtriplot('horizontal 4')
qtriplot('vertical 9')
x=1;
y=1;
orgname =[];
pline=0; % y for use in script

fid=fopen([subject '_' num2str(sub) '.results'],'wt');
for ibeat=1:length(bsmfiles)
    bsmfile = [bsmdir bsmfiles(ibeat).name];
    %     load([dirout initdatafiles(ibeat).name ])
    %     load([dirout findatafiles(ibeat).name ])
    beat = bsmfiles(ibeat).name(1:end-7); % remove '.selecg'
    electrodes=GetCatheterPosition(sub,beat); % get electrode positions
    dirinit=dir(fullfile(dirout,[subject '_' beat '*init.mat']));
    try
        load(fullfile(dirout,dirinit.name)); % oostep1: error when 0 or more than 1. Catch error??
        load(fullfile(dirout,regexprep(dirinit.name,'init.mat','.mat')));
    catch
        warning(['skipped file:' bsmfiles(ibeat).name]);
        continue
    end
    
    prevname = orgname;
    orgname = beat(1:strfind(beat,'beat')-2);
    load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'])
    remove = DATA.remove;
    remove=DATA.remove(1:64); % correction for including EGM in signalss but not layout
    disp('limited remove on 64');
    
    
    clear DATA;
    % Load all data from selected subject
    disp('===============================================================')
    disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
    GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',1,'group',group,'basedir',geomdir);
    %     GEOM.BSM=lowpassma(GEOM.BSM,10);
    
    GEOM=prepare_geom(GEOM,[bsmfile(1:end-7) '.ecgspecs'],1 );
%     GEOM=prepare_geom(GEOM,[bsmfile(1:end-7) '.spe'],1 );
    
    if  GEOM.specs(5) > 1000
        error('hallo')
    end
    GEOM.beat=beat;
    useLeads = 1:size(GEOM.LAY,1) - 1;
    useLeads(remove ==1) = [];
    L = GEOM.LAY(2:end,:);
    for i=1:size(L,1)
        L(i,3) = L(i,3) - sum(remove (1:i));
        
    end
    L(remove ==1,:)=[];
    GEOM.LAY = [GEOM.LAY(1,:); L];
    GEOM = selectLeads(GEOM,useLeads,1);
    
    %WCT
    if ~isempty(wctorg)
        wctnew=[];
        for i=1:length(wctorg)
            wctnew(i)=find(L(:,3)==(wctorg(i)- sum(remove (1:wctorg(i)))));
        end
        GEOM.BSM=bsxfun(@minus,GEOM.BSM,mean(GEOM.BSM));
    end
    
    % GEOM.BSM = GEOM.BSM * 7;
    %     close all;
    
    
    qtriplot(GEOM.VER,GEOM.ITRI);
    qtriplot(meas.depfinal-min(meas.depfinal));
    %     qtriplot(measinit.dep - min(measinit.dep));
    qtriplot('step 5');
    
    
    
    if (x==1 && y==1) || strcmp(orgname,prevname)
        qtriplot(['panel ' num2str([x,y])]);
        pause(0.5);
        qtriplot(['panel ' num2str([x,y])]);
        x=x+1;
    else
        y=y+1;
        x=1;
        qtriplot(['panel ' num2str([x,y])]);
        pause(0.5);
        qtriplot(['panel ' num2str([x,y])]);
        x=2;
    end
    
    pline=y; % possibility to disconnect later
    % y is current line, x is next column
    if makescript
        fnamebeat=sprintf('ID%d_Beat%d',sub,ibeat);
        savetri(fullfile(scriptdir,[fnamebeat '.tri']),GEOM.VER,GEOM.ITRI);
        fun=meas.depfinal-min(meas.depfinal);
        funb=measinit.dep-min(measinit.dep);
        saveasci(fullfile(scriptdir,[fnamebeat '_opt.fun']),fun);
        saveasci(fullfile(scriptdir,[fnamebeat '_init.fun']),funb);
        if ~strcmp(orgname,prevname)
            % init
            %             pline=pline+1; not needed when using y
            qtriscriptl(pline).script={'delete *',['horizontal ' num2str(ncol+1)],'vertical 2',}; %init for per line script
            [dummy1, oname, dummy2] = fileparts(orgname);
            qtriscriptl(pline).fname=oname;
        end
        panelnames={};
        
        dmat=0.5*[0 0 0;1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];% make some noise!
        [dummy focioptn]=min(meas.depfinal);
        [dummy focibeginn]=min(measinit.dep);
        fociopt=GEOM.VER(focioptn,:);
        focibegin=GEOM.VER(focibeginn,:);
        savetri(fullfile(scriptdir,[fnamebeat '_focibegin.tri']),bsxfun(@plus,dmat,focibegin),[]);
        savetri(fullfile(scriptdir,[fnamebeat '_fociopt.tri']),bsxfun(@plus,dmat,fociopt),[]);
        try
            copyfile(fullfile(geomdir,subject,[subject '_lad.tri']),fullfile(scriptdir,[fnamebeat '_lad.tri']));
        end
        if focibegin==fociopt
            fbcolor='purple';
            focolor='red'; % was purple, with tim maybe orange?
        else
            fbcolor='green';
            focolor='red';% was red using heat.
        end
        resultscript={...
            ...sprintf('file %s=%s',[fnamebeat '_focibegin'],[fnamebeat '_focibegin.tri']),...
            ...sprintf('color %s=%s',[fnamebeat '_focibegin'],fbcolor),...
            sprintf('file %s=%s',[fnamebeat '_fociopt'],[fnamebeat '_fociopt.tri']),...
            sprintf('color %s=%s',[fnamebeat '_fociopt'],focolor),...
            sprintf('file %s=%s',[fnamebeat '_lad'],[fnamebeat '_lad.tri']),...
            sprintf('color %s=%s',[fnamebeat '_lad'],'black'),...
            };
        %             panelnames=[panelnames;{[fnamebeat '_focibegin'],[fnamebeat '_fociopt'],[fnamebeat '_electrodes']}];
        panelnames=[panelnames;{[fnamebeat '_fociopt'],[fnamebeat '_lad']}];
        if ~isempty(electrodes) %
            savetri(fullfile(scriptdir,[fnamebeat '_electrodes.tri']),bsxfun(@plus,dmat,electrodes),[]);
            panelnames=[panelnames {[fnamebeat '_electrodes']}];
            resultscript=[resultscript {...
                sprintf('file %s=%s',[fnamebeat '_electrodes'],[fnamebeat '_electrodes.tri']),...
                sprintf('color %s=blue',[fnamebeat '_electrodes'])...
                }];
            
        end
        geoscript={sprintf('file %s=%s',sprintf('beat%d',ibeat),[fnamebeat,'.tri']),...% plot geoemtries with activation
            sprintf('fun %s=1',sprintf('beat%d',ibeat)),...
            'step 10',...% 5 for pigs.
            'funscale 0,150'}; % not set (0,80?) for pigs
        geoscript{end+1}=sprintf('funfile %s=%s',sprintf('beat%d',ibeat),[fnamebeat '_opt.fun']);
        
        panelnames=[panelnames {sprintf('beat%d',ibeat)}];
        
        % scripts for moving to the right panel
        panellscript={};
        panelallscript={};
        for i=1:length(panelnames)
            panelallscript=[panelallscript {sprintf('panel %s=%d %d',panelnames{i},1+x-1 ,y)}]; % pig+signal first column
            panellscript=[panellscript {sprintf('panel %s=%d %d',panelnames{i},1+1+mod(x-2,ncol),ceil((x-1)/ncol))}]; % pig+signal first column
        end
        if x==2 % x==2, first beat from file
            % copy thorax, script thorax, read and write beat signal, movie
            if y==1 % only once
                geoscript{end+1}='funcolor tim';
                [VERt,ITRIt]=loadtri(fullfile(geomdir,subject,[subject thoraxname]));
                shift=mean(GEOM.VER,1)-mean(VERt*thoraxscale,1);
                savetri(fullfile(scriptdir,'thorax.tri'),bsxfun(@plus,VERt*thoraxscale,shift),ITRIt);
                savetri(fullfile(scriptdir,'origin.tri'),bsxfun(@plus,dmat,shift),[1 2 3]); % dummy geometry to hold rms
                copyfile('64.el',fullfile(scriptdir,'64.el'));
                
            end
            
            T=intripol(VERt,ITRIt,useLeads);
            VALS= T*GEOM.BSM;
            
            bsmpath=fullfile(scriptdir,[qtriscriptl(pline).fname '.asc']);
            VALSscale=(VALS+1.5)*80.0/3;
            %             saveasci(bsmpath,VALSscale);
            saveasci(bsmpath,VALS);
            
            bsmrms=rms(GEOM.BSM);
            rmspath=fullfile(scriptdir,[qtriscriptl(pline).fname '_rms.asc']);
            bsmrmsscale=(bsmrms+1.5)*80.0/3;
            %             saveasci(rmspath,bsmrmsscale);
            saveasci(rmspath,bsmrms);
            
            bsmscript={'file thorax=thorax.tri','panel thorax=1 1',['funfile thorax=' qtriscriptl(pline).fname '.asc'],...
                'step2 .2','fun thorax=2','funcolor2 heat','funscale2 -1.5,1.5','mouse thorax=vertex',...
                'elfile thorax=64.el',...'marker thorax=vertex 0','signal thorax=0,0.05,0.2,0.15,0.2'
                'file dummy=origin.tri','panel dummy=1 1','fun dummy=2',['funfile dummy=' qtriscriptl(pline).fname '_rms.asc'],...
                'marker dummy=red','signal dummy=1,0.05,0.2001,0.15,0.2'};
            
            % signal, init speed position,mouse
        else
            bsmscript={};
        end
        
        %build scripts
        qtriscriptall=[qtriscriptall  geoscript resultscript panelallscript];
        qtriscriptl(pline).script=[qtriscriptl(pline).script geoscript resultscript panellscript bsmscript];
    end
    
    
    
    
    %     ab=find(measinit.dep == min(measinit.dep)); % oostep1 crashes when
    %     % more points at same act time
    %     ae=find(meas.depfinal == min(meas.depfinal));
    %     A=[ab GEOM.VER(ab,:) measinit.rd measinit.cor ae GEOM.VER(ae,:) meas.rdfinal meas.corfinal];
    %     fprintf(fid,'%s\n',[beat '  ' num2str(A)]);
    
    %      figure(pp+10);showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',measinit.foci);view(0,0);
    %             fh=figure(pp+10);set(fh,'Name','Final: Activation time');showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',find(meas.depfinal==min(meas.depfinal)));view(0,0);
    %         figure(pp+20);showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal,'nodes',measinit.foci);view(0,0);
    %         figure(pp+30);showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal,'nodes',measinit.foci);view(0,0);
    %         figure(pp+31);clf;plot(meas.depfinal(GEOM.Rfreewallver==0),meas.repfinal(GEOM.Rfreewallver==0)-meas.depfinal(GEOM.Rfreewallver==0),'.r');
    %                   hold on;plot(meas.depfinal(GEOM.Rfreewallver==1),meas.repfinal(GEOM.Rfreewallver==1)-meas.depfinal(GEOM.Rfreewallver==1),'.b')
    %         t=0:GEOM.specs(5)-GEOM.specs(2);
    %         T=ones(length(GEOM.VER),1)*t;
    %         PSIA =lowpassma(GEOM.AMA*getSmode(T,meas.depfinal,meas.repfinal,GEOM.pS,[],4),lpass);
    %             figure(32);clf
    %             sigplot(GEOM.BSM(1:size(GEOM.LAY,1)-1,GEOM.specs(2):GEOM.specs(5)),'',GEOM.LAY,1,'b',1,0);
    %             leadv16(GEOM.BSM,'leadsystem','ams' );
    %             title(beat)
    %             saveas(gcf,['./figs/pig07/' GEOM.subject beat '.png' ]);
    %         hold on
    %         sigplot(PSIA(1:size(GEOM.LAY,1)-1,:),'',GEOM.LAY,1,'r',1,0);
    %         hold on
    
end
fclose(fid);
% pause
% % qtriplot('panel 4 9')
% qtriplot('bgdcolor white')
% qtriplot('funscale 0 110')
% qtriplot('funcolor hsv');
% qtriplot('step 10')
% qtriplot(['png ./figs/pig07/' GEOM.subject num2str(sub) '.png'])
% pause
% qtriplot('funscale 0 20')
% qtriplot('step 2.5')
% qtriplot(['png ./figs/pig07/' GEOM.subject num2str(sub) 'focus.png'])

% save scripts
fid=fopen(fullfile(scriptdir,[fnameall '.trp']),'w+');
for i=1:length(qtriscriptall)
    fprintf(fid,'%s\n',qtriscriptall{i});
end
fclose(fid);



for j=1:length(qtriscriptl)
    if ~isempty(wctnew)
        qtriscriptl(j).fname=[qtriscriptl(j).fname '_WCT'];
    end
    fid=fopen(fullfile(scriptdir,[qtriscriptl(j).fname '_opt.trp']),'w+');
    for i=1:length(qtriscriptl(j).script);
        fprintf(fid,'%s\n',qtriscriptl(j).script{i});
    end
    fclose(fid);
end

for j=1:length(qtriscriptl)
    fid=fopen(fullfile(scriptdir,[qtriscriptl(j).fname '_init.trp']),'w+');
    for i=1:length(qtriscriptl(j).script);
        fprintf(fid,'%s\n',regexprep(qtriscriptl(j).script{i},'_opt.fun','_init.fun'));
    end
    fclose(fid);
end

disp(['finished ' mfilename]);

