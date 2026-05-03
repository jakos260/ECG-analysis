
geomdir = fullfile('C:\Users\Damp2\Documents\ECG_simulation\STW\Data\geometries\');
casedir		='';
dirout = fullfile('.\results\');
lpass=5;   % # samples in lowpassma used to filter reults
clusterDist=30; % used in multifocisacn

% which leads should be used in the initial esitimate
leadset='all';
group=[];
sub=8;
saveCase =1; % 2 = overschrijven
sinkScan= 0;
layfile='pig64.mla';
type='ventricles';

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
    subject='Pig09_refined';
        bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\Pig09PacedBeats\beats\';
    bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\export\beats\';
    bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\newECG\beats1\';
    bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\Pig09Export\AVG\beats\';
    dirout= ['.\results\' subject '\'];
%     bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\avobeats\';
%     dirout= ['.\results\' subject '\avobeats\'];
else
    return
end

dirout= ['.\results\' subject '\ESC\'];

clusters = 1;
useScaling =1;
if clusters == 1 % 25 mm
    dirout = fullfile([dirout 'single\']);
end
if ~exist(dirout,'dir')
    mkdir(dirout)
end
bsmfiles = dir([bsmdir '*.selecg']);

RESULTS =[];


diroutOrg = dirout;
velo = 1.0;
anis = 2.5;

for ibeat = 1 : length(bsmfiles)
    
    clear meas;
    clear measinit;
    bsmfile = [bsmdir bsmfiles(ibeat).name];
    beat = bsmfiles(ibeat).name(1:end-7);
    orgname = beat(1:strfind(beat,'beat')-2);
    if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file')
        load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'])
        remove = DATA.remove(1:64);
        clear DATA;
    else
        A=loadmat(bsmfile);
        remove = zeros(64,1);
    end
    % Load all data from selected subject
    disp('===============================================================')
    disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
    %
    GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',anis,'group',group,'basedir',geomdir);
        
    dirout = diroutOrg;
    heartpath=fullfile(geomdir, group,  subject,  subject );
    
    if ~exist(dirout,'dir')
        mkdir(dirout)
    end
    
    GEOM.beat=beat;
    useLeads=1:size(GEOM.LAY,1) - 1;
    useLeads(remove ==1)=[];
    L= GEOM.LAY(2:end,:);
    for i=1:size(L,1)
        L(i,3) = L(i,3) - sum(remove (1:i));
    end
    L(remove ==1,:)=[];
    GEOM.LAY = [GEOM.LAY(1,:); L];
    
    
%     useLeads= [40 41 42 46 47 48 52 53 54];
%         useLeads= [38 40 42 43 45 47 49 51 53];
% %     useLeads = round(rand(9,1) * (length(useLeads)-1));
%     GEOM.LAY= loadmat('9lds.mla');
    
    GEOM.SPECS = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase);
    GEOM = selectLeads(GEOM,useLeads,1);
    
%     GEOM = prepare_geom(GEOM,[bsmfile(1:end-7) '.spe'],saveCase,GEOM.ECGextra );    

    %                 GEOM.BSM(:,GEOM.specs(2):end) = killhum(GEOM.BSM(:,GEOM.specs(2):end),50,1000,0.05);
    %                 GEOM.BSM = lowpassma(GEOM.BSM,20);
    
%     GEOM.specs(6)=0.005513;
%     GEOM.specs(7)=0.055512;
%     GEOM.specs(8)=1.317503;
%     GEOM.pS = GEOM.specs(6:8);
    
   
%     eval(['BSM' num2str(ibeat) '=baselinecor(GEOM.BSM(:,GEOM.SPECS.q(2):end),1,GEOM.specs(5)-GEOM.specs(2));']);
    
    close all;
    drawnow
    
    clusters =1;
    
    % initial scan
    disp(['anisotropyRatio: ' num2str(GEOM.anisotropyRatio)])

%     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_wall(GEOM,'clusters',clusters);
%     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_v12(GEOM,'clusters',clusters,'usecor',1,'usetime',60,'velocity',0.5);
%     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_pvd(GEOM,'clusters',clusters,'usecor',1,'usetime',60,'velocity',0.4);
%     [measinit.foci,measinit.dep,measinit.outp] = multifociscan_pvd(GEOM,'clusters',6,'usecor',1,'issinus',1,'initialvelocity',0.4,'maxvelocity',0.9);
    [measinit.foci,measinit.dep,measinit.outp] = multifociscan_esc2013(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',0.4,'maxvelocity',1.3);


    measinit.anisotropyRatio = GEOM.anisotropyRatio;
    measinit.cor = measinit.outp(end,1);
    measinit.rd  = measinit.outp(end,2);
    measinit.rep = initRep(GEOM,measinit.dep);
    %%
    pp=30;
    meas=inverse_pvd(GEOM,measinit.dep ,measinit.rep,'estimateampl',0,...
        'casedir',dirout,...
        'repopt','apd',... %'rep'
        'maxiter',40,...
        'mudep',1.5e-6,...
        'murep',1.5e-6,...
        'weighed',0,...
        'minrd',0.15,...
        'mode',4);
    figure(pp+10);showpatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',[find(meas.depfinal==min(meas.depfinal(GEOM.typ==1))) find(meas.depfinal==min(meas.depfinal(GEOM.typ==2)))]);view(0,0);
    figure(pp+11);showpatch(GEOM.VER,GEOM.ITRI,meas.repfinal)
    figure(pp+12);showpatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal)
    
%     writeresults;
    RESULTS = [RESULTS;   GEOM.VER(meas.depfinal==min(meas.depfinal),:) ];
    
%     figure(pp+10);showpatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',measinit.foci);view(0,0);
%     figure(pp+10);showpatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',find(meas.depfinal==min(meas.depfinal)));view(0,0);
    %%
%     figure(pp+20);showpatch(GEOM.VER,GEOM.ITRI,meas.repfinal,'nodes',measinit.foci);view(0,0);
%     figure(pp+30);showpatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal,'nodes',measinit.foci);view(0,0);
%     figure(pp+31);clf;plot(meas.depfinal(GEOM.Rfreewallver==0),meas.repfinal(GEOM.Rfreewallver==0)-meas.depfinal(GEOM.Rfreewallver==0),'.r');
%     hold on;plot(meas.depfinal(GEOM.Rfreewallver==1),meas.repfinal(GEOM.Rfreewallver==1)-meas.depfinal(GEOM.Rfreewallver==1),'.b')
    t=0:GEOM.SPECS.endtwave;
    T=ones(length(GEOM.VER),1)*t;
    PSIA =lowpassma(GEOM.AMA*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,[],4),lpass);
    figure(pp+32);clf
    sigplot(GEOM.BSM,'',GEOM.LAY,1.5,'b',1,0);
    hold on
    sigplot(PSIA,'',GEOM.LAY,1.5,'r',1,0);
    hold on
%     PSIA = PSIA(:,1:GEOM.qrsduration);
%     BSM=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(2)+size(PSIA,2)-1);
%     meas.rdfinalqrs =  norm(BSM - PSIA,'fro')/norm(BSM,'fro');
    
    
    
    meas.velocity = velo;
    meas.anis=anis;
%     meas.vt=vt;
    measinit.velocity = velo;
%     measinit.vt = vt;
    measinit.anis=anis;
    if saveCase
        save([dirout GEOM.subject '_' beat '_' type '_' num2str(velo) '_' num2str(anis) '_' date 'init.mat'],'measinit')
        save([dirout GEOM.subject '_' beat '_' type '_' num2str(velo) '_' num2str(anis) '_' date '.mat'],'meas')
    end
    
    %                 A=GEOM.VER(find(meas.depfinal == min(meas.depfinal)),:);
    %                 xlswrite(GEOM.subject,A,['velocity' num2str(velo) 'anis' num2str(anis)],['A' num2str(ibeat)]);
    
    %                 geomheart = 'C:\Inge\invedl\dataVarken08\model\Pig08\Pig08_ventricles.tri';
    %                 [VER, ITRI] = loadtri(geomheart);
    %                 [p,q] =min(meas.depfinal);
    %                 xlswrite(num2str(GEOM.beat),[num2cell(GEOM.beat),' ',num2str(measinit.rd),' ',num2str(measinit.cor),' ',num2str(measinit.foci),' ',num2str(GEOM.VER(measinit.foci,:)),' ',num2str(meas.rdfinal),' ',num2str(meas.corfinal),' ',q,' ', num2str(GEOM.VER(q,:))    ]);
    %
    
end


