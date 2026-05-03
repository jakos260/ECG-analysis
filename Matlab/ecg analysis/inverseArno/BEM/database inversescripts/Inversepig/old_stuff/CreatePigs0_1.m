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
skipexisting=1; % skip if result files already present in target folder ( recover from crash)
sinkScan= 0;
layfile='pigs_adam.mla';%'pig64.mla';
type='ventricles';
clusters = 1;
wrd=1; % use weighed rd in inverse
useScaling =1;
dostampl=false; % fit AP amplitude on ST segmet.
% initialActTime=0; % override. 0 means scaling with maxvelocity and anis. Make [] for automatic setting in multifociscan
if wrd
    mudep=2.5e-6;
    murep=2.5e-6;
else
    mudep=1e-5;
    murep=1e-5;
end

INITIALVELOCITY=[0.2,0.4,0.6,0.8,1.0];
maxvelocity=1.6;
velo = 1.0;
ANIS = [0.1,0.2,0.3,0.4,0.5,0.8,1,1.5];
INITANIS=[0.5 ,0.8, 1.0,1.5,2.5];


%initialtype

initmode=1; % methods for intial estimate
% 0: use normal fociscan
% 1: test with fixed initial estimate focus at stimulation site
% n ('number'): test nearest n points at same surface
% negative: test with focus on oposite surface

for initialvelocity=INITIALVELOCITY;
    for initanis=INITANIS;
        for anis=ANIS 





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
%     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/SRST/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    %     bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\avobeats\';
    %     dirout= ['.\results\' subject '\avobeats\'];
else
    return
end


if clusters == 1 % 25 mm
    dirout = fullfile(dirout, 'single', sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.1f',datestr(now,'yyyymmdd'),initmode,wrd,initialvelocity,initanis,anis));
end
if ~exist(dirout,'dir')
    mkdir(dirout)
end
bsmfiles = dir([bsmdir '*.selecg']);

RESULTS =[];


diroutOrg = dirout;

for ibeat =  1 : length(bsmfiles)
    
    clear meas;
    clear measinit;
    bsmfile = [bsmdir bsmfiles(ibeat).name];
    beat = bsmfiles(ibeat).name(1:end-7);
    orgname = beat(1:strfind(beat,'beat')-2);
    fnameout=[subject '_' beat '_' type '_' num2str(velo) '_' num2str(anis)];
    
    if ~skipexisting || length(dir(fullfile(dirout, [fnameout '*.mat'])))<2
        
        
        
        if exist([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'],'file')
            load([bsmdir(1:strfind(bsmdir,'beats')-1) orgname '.mat'])
            remove = DATA.remove(1:64);
            clear DATA;
        else
            error('Export file not available');
            %         A=loadmat(bsmfile);
            %         remove = zeros(64,1);
        end
        % Load all data from selected subject
        disp('===============================================================')
        disp(['Loading the ' type  ' of subject ' subject '  selected beat:' beat ])
        %
        GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',anis,'group',group,'basedir',geomdir);
        if initanis~=2.5 % override hardcoded 2.5 in multifociscan
            if initanis==1
                GEOM.DIST25=GEOM.DIST;
            else
                GEOM.DIST25=loadmat(fullfile(geomdir,subject,[subject num2str(initanis) '_' type '.dist']));
            end
        end
        
        
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
        GEOM = selectLeads(GEOM,useLeads,1);
        GEOM.SPECS = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim',1);
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
        
        if dostampl
            % Amplitude
            rmsbsm=rms(GEOM.BSM);
            [~,ist]=min(rmsbsm((25:100)+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint));
            ist=ist+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint;
            figure;plot(rms(GEOM.BSM))
            hold on
            plot(ist,rmsbsm(ist),'r*');
            %         ist=350; % defaultposition for ST amplitude estimation.
            Ampl=sqrAmpl(GEOM,ist);
            qtriplot('reset');
            qtriplot(GEOM.VER,GEOM.ITRI,'Ventricles');
            qtriplot(Ampl);
            qtriplot([-16.8	27.3	888.3],[],'LVLatENdoII');
        end
        
        
        
        close all;
        %         drawnow
        
        %         clusters =1;
        
        % initial scan
        disp(['anisotropyRatio: ' num2str(GEOM.anisotropyRatio)])
        surfnames={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % Simplified
        
        for iinit=1:max(abs(initmode),1)
            if initmode==0
                %     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_wall(GEOM,'clusters',clusters);
                %     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_v12(GEOM,'clusters',clusters,'usecor',1,'usetime',60,'velocity',0.5);
                %     [measinit.foci,measinit.dep,measinit.outp]=multifociscan_pvd(GEOM,'clusters',clusters,'usecor',1,'usetime',60,'velocity',0.4);
                %     [measinit.foci,measinit.dep,measinit.outp] = multifociscan_pvd(GEOM,'clusters',6,'usecor',1,'issinus',1,'initialvelocity',0.4,'maxvelocity',0.9);
                %                 [measinit.foci,measinit.dep,measinit.outp] = multifociscan_esc2013(GEOM,'clusters',clusters,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',maxvelocity,'initialacttime',initialActTime);
                [measinit.foci,measinit.dep,measinit.outp] = multifociscan_publicationPigs(GEOM,'clusters',clusters,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',maxvelocity);
            else
                [stimpos,cathsurf]=GetCatheterPosition2(sub,bsmfiles(ibeat).name);
                if initmode>0
                    [pnearest, distp, trinearest, focus, distver]=findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,initmode);
                    
                else
                    oppsurf=unique(surfnames(~strcmp(cathsurf,surfnames)));
                    [pnearest, distp, trinearest, focus, distver]=findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,oppsurf,-initmode);
                end
                [measinit.foci,measinit.dep,measinit.outp] = ...
                    multifociscan_publicationPigs(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',maxvelocity,'focus',focus);
                %                     multifociscan_esc2013(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',maxvelocity,'initdep',focus,'initialacttime',initialActTime);
                
%                 measinit.outp=NaN(1,3);
            end
            
            disp(['min deptime ' num2str(min(measinit.dep))])
            
            measinit.anisotropyRatio = GEOM.anisotropyRatio;
            measinit.cor = measinit.outp(end,1);
            measinit.rd  = measinit.outp(end,2);
            measinit.rep = initRep(GEOM,measinit.dep);
            %%

            pp=30;
            meas=inverse_pvd(GEOM,measinit.dep ,measinit.rep,'estimateampl',0,...
                'casedir',dirout,...
                'repopt','apd',... %'rep'
                'maxiter',1,...
                'mudep',mudep,...
                'murep',murep,...
                'weighed',0,...
                'minrd',0.15,...
                'mode',4,...
                'weighedrd',wrd);
            figure('Name','Final activation');showPatch(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',[find(meas.depfinal==min(meas.depfinal(GEOM.typ==1))) find(meas.depfinal==min(meas.depfinal(GEOM.typ==2)))]);view(0,0);
            figure('Name','Final repolarization');showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal)
            figure('Name','Final ARI');showPatch(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal)
            
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
            PSIA =lowpassma(GEOM.AMA*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,4),lpass);
            figure('Name','Final QRST');clf
            sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),'',GEOM.LAY,1.5,'b',1,0);
            hold on
            sigplot(PSIA,'',GEOM.LAY,1.5,'r',1,0);
            hold on
            %     PSIA = PSIA(:,1:GEOM.qrsduration);
            %     BSM=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(2)+size(PSIA,2)-1);
            %     meas.rdfinalqrs =  norm(BSM - PSIA,'fro')/norm(BSM,'fro');
            
            
            
            meas.velocity = velo;
            meas.anis=anis;
            meas.mudep=mudep;
            meas.murep=murep;
            %     meas.vt=vt;
            measinit.velocity = velo;
            %     measinit.vt = vt;
            measinit.anis=anis;
            if saveCase
                save(fullfile(dirout, [fnameout '_' date 'init.mat']),'measinit')
                save(fullfile(dirout ,[fnameout '_' date '.mat']),'meas')
                if dostampl
                    save(fullfile(dirout, [GEOM.subject '_' beat '_' type '_IST' num2str(ist) '_' num2str(anis) '_' date 'stampl.mat']),'Ampl')
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
        end
    end
end
        end
    end
end


