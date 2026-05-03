% function showPigsPO()
global lpass
global clusterDist
clear qtriscriptl
makescript=true; % create qtriplot script. One for all, one for each line/orgname.

qtriscriptall={'delete *','horizontal 6','vertical 9',};
ncol=3; % number of columns in per line view. number of row is hardcode in the next line:
qtriscriptlinit={'delete *',['horizontal ' num2str(ncol)],'vertical 2',}; %init for per line script

geomdir = '/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/geometries';


casedir		='';
dirout='./results/';
lpass=5;   % # samples in lowpassma used to filter reults
clusterDist=30; % used in multifociscan

% which leads should be used in the initial esitimate
leadset='all';
group=[];
sub=71;
saveCase =1;
sinkScan= 0;
layfile='pig64.mla';
type='ventricles';

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
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/beats/';
    dirout= ['./results/' subject '/'];
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
else
    return
end

[orgnames,focibegin,fociopt,electrodes]=getmarkers(sub); % get electrode positions etc.


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
    load([dirout initdatafiles(ibeat).name ])
    load([dirout findatafiles(ibeat).name ])
    
    beat = bsmfiles(ibeat).name(1:end-7);
    prevname = orgname;
    orgname = beat(1:strfind(beat,'beat')-2);
    load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'])
    remove = DATA.remove;
    clear DATA;
    % Load all data from selected subject
    disp('===============================================================')
    disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
    GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',1.5,'group',group,'basedir',geomdir);
    %     GEOM.BSM=lowpassma(GEOM.BSM,10);
    
    GEOM=prepare_geom(GEOM,[bsmfile(1:end-7) '.spe'],1 );
    
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
        saveasci(fullfile(scriptdir,[fnamebeat '.fun']),fun);
        if ~strcmp(orgname,prevname)
            % init
            %             pline=pline+1; not needed when using y
            qtriscriptl(pline).script=qtriscriptlinit;
            [dummy1, oname, dummy2] = fileparts(orgname);
            qtriscriptl(pline).fname=oname;
        end
        panelnames={};
        if length(orgnames) % empty is no results present
            savetri(fullfile(scriptdir,[fnamebeat '_focibegin.tri']),focibegin(ibeat,:),[]);
            savetri(fullfile(scriptdir,[fnamebeat '_fociopt.tri']),fociopt(ibeat,:),[]);            
            savetri(fullfile(scriptdir,[fnamebeat '_electrodes.tri']),electrodes(ibeat,:),[]);
            if focibegin(ibeat,:)==electrodes(ibeat,:)
                fbcolor='green';
            else fbcolor='lightgreen'; end
             if fociopt(ibeat,:)==electrodes(ibeat,:)
                focolor='purple';
            else focolor='red'; end
            resultscript={sprintf('file %s=%s',[fnamebeat '_focibegin'],[fnamebeat '_focibegin.tri']),...
                sprintf('color %s=%s',[fnamebeat '_focibegin'],fbcolor),...
                sprintf('file %s=%s',[fnamebeat '_fociopt'],[fnamebeat '_fociopt.tri']),...
                sprintf('color %s=%s',[fnamebeat '_fociopt'],focolor),...
                sprintf('file %s=%s',[fnamebeat '_electrodes'],[fnamebeat '_electrodes.tri']),...
                sprintf('color %s=black',[fnamebeat '_electrodes'])};
            panelnames=[panelnames;{[fnamebeat '_focibegin'],[fnamebeat '_fociopt'],[fnamebeat '_electrodes']}];
        else
            resultscript=[];
        end
        geoscript={sprintf('file %s=%s',sprintf('beat%d',ibeat),[fnamebeat,'.tri']),...
            sprintf('funfile %s=%s',sprintf('beat%d',ibeat),[fnamebeat '.fun']),...
            'step 5'}; % plot geoemtries
        panelnames=[panelnames {sprintf('beat%d',ibeat)}];
        
        % scripts for moving to the right panel
        panellscript={};
        panelallscript={};
        for i=1:length(panelnames)
            panelallscript=[panelallscript {sprintf('panel %s=%d %d',panelnames{i},x-1 ,y)}];
            panellscript=[panellscript {sprintf('panel %s=%d %d',panelnames{i},1+mod(x-2,ncol),ceil((x-1)/ncol))}];
        end
        
        %build scripts
        qtriscriptall=[qtriscriptall  geoscript resultscript panelallscript];
        qtriscriptl(pline).script=[qtriscriptl(pline).script geoscript resultscript panellscript];
    end
    
    
    
    
    ab=find(measinit.dep == min(measinit.dep));
    ae=find(meas.depfinal == min(meas.depfinal));
    A=[ab GEOM.VER(ab,:) measinit.rd measinit.cor ae GEOM.VER(ae,:) meas.rdfinal meas.corfinal];
    fprintf(fid,'%s\n',[beat '  ' num2str(A)]);
    
    %      figure(pp+10);showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',measinit.foci);view(0,0);
    %         figure(pp+10);showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',find(meas.depfinal==min(meas.depfinal)));view(0,0);
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
    fid=fopen(fullfile(scriptdir,[qtriscriptl(j).fname '.trp']),'w+');
    for i=1:length(qtriscriptl(j).script);
        fprintf(fid,'%s\n',qtriscriptl(j).script{i});
    end
    fclose(fid);
end




