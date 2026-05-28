
global lpass
global clusterDist

geomdir = '/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Results1/geometries';


casedir		='';
dirout='./results/';
lpass=5;   % # samples in lowpassma used to filter reults
clusterDist=30; % used in multifocisacn

% which leads should be used in the initial esitimate
leadset='all';
group=[];
sub=7;
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
dirout=fullfile(dirout);
bsmfiles = dir([bsmdir '*.selecg']);
initdatafiles = dir([dirout '*init.mat']);
findatafiles = dir([dirout '*2.mat']);
qtriplot('delete *')
qtriplot('horizontal 4')
qtriplot('vertical 9')
x=1;
y=1;
orgname =[];

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


    qtriplot(GEOM.VER,GEOM.ITRI)
        qtriplot(meas.depfinal-min(meas.depfinal));
%     qtriplot(measinit.dep - min(measinit.dep));
    qtriplot('step 5');


    if (x==1 && y==1) || strcmp(orgname,prevname)
        qtriplot(['panel ' num2str([x,y])])
        pause(0.5)
        qtriplot(['panel ' num2str([x,y])])
        x=x+1;       
    else
        y=y+1;
        x=1;
        qtriplot(['panel ' num2str([x,y])])
        pause(0.5)
        qtriplot(['panel ' num2str([x,y])])
        x=2;
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
pause
% qtriplot('panel 4 9')
qtriplot('bgdcolor white')
qtriplot('funscale 0 110')
qtriplot('funcolor hsv');
qtriplot('step 10')
qtriplot(['png ./figs/pig07/' GEOM.subject num2str(sub) '.png'])
% pause
% qtriplot('funscale 0 20')
% qtriplot('step 2.5')
% qtriplot(['png ./figs/pig07/' GEOM.subject num2str(sub) 'focus.png'])
