
global lpass
global clusterdist


if ispc
    geomdir = 'C:\Users\Peter.Damp2\Documents\Data\geometries\'; % Peter van Dam's laptop
else
    geomdir = '/Users/peteroosterhoff/Documents/Werk/Brugada/DATA/geometries'; % Peter O's MacBook
end;



casedir		='';

dirout=fullfile('.','resuls'); %'.\results\';

lpass=5;   % # samples in lowpassma used to filter reults
clusterdist=30; % used in multifocisacn

% which leads should be used in the initial esitimate
leadset='all';

sub=2;
savecase = 1;
sinkscan= 0;
doatria=0;
group = 'brugada';layfile='ams65.mla';type='ventricles';
if sub==1
    subject='ajm82';
    bsmdir='c:\users\peter.damp2\documents\data\measurements\brugada\ajm082\bspm\export\beats\';
elseif sub ==2
    subject='ajm059';
    bsmdir='/users/peteroosterhoff/documents/werk/brugada/data/measurements/ajm059/bspm/export/beats/';
    bsmdir='C:\Users\Peter.Damp2\Documents\Data\measurements\Brugada\Brug59\ECG\avgbeats\' 
elseif sub ==3
    subject='ajm119';
    bsmdir='c:\users\peter.damp2\documents\data\measurements\brugada\ajm119\ecg\bspm\export\beats\';
else
    return;
end
%%
bsmfiles = dir([bsmdir  '*.beatecg']);
bsmfiles = dir([bsmdir  '*.mes']);
for ibeat = 15:length(bsmfiles);
    dirout= fullfile('.','results',subject); %['./results/' subject '\'];
    if ~exist(dirout,'dir')
        mkdir(dirout)
    end
    
    beat = bsmfiles(ibeat).name(1:end-7);
    bsmfile = [bsmdir bsmfiles(ibeat).name];
    % load all data from selected subject
    disp('===============================================================')
    disp(['loading the ' type  ' of subject ' subject '  selected beat:' beat ])
    geom=invinit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyratio',2.5,'group',group,'basedir',geomdir,'usemean',1);
    %     [geom.bsm,geom.spikesamples]=removepacingspike(geom.bsm);
    %     geom = prepare_geom(geom,[bsmdir beat '.spe'],savecase );
    geom = prepare_geom(geom,[bsmfile '.spe'],savecase );
    %     [geom.bsm,geom.spikesamples]=removepacingspike(geom.bsm);
    if doatria
        geoma = invinit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type','atria','anisotropyratio',1,'group',group,'basedir',geomdir,'usemean',1);
        % geoma=prepare_geom(geoma,fullfile(geomdir, group, geom.subject , [geom.subject beat 'atria' '.spe']),savecase );
        geoma.specs = [0.25 1 110 geom.specs(2) geom.specs(5:end)']';
        geoma.bsm= geom.bsm;
        geoma.ps=geom.ps;
        geoma.beat=beat;
    end
    geom.beat=beat;
    
    % close all;
    drawnow;
    eval(['ECG' num2str(ibeat) '= geom.BSM;']);
    %%
    for k=2:-1:1
        
        if k ==0
            useleads = 1:length(geom.tver);
            leadsys  = 'intripol';
            mu       = 1.5e-4;
        elseif k == 1
            useleads = geom.v12;
            leadsys  = 'v12';
            mu       = 0.5e-4;
        else
            useleads = 1:size(geom.LAY,1) - 1; % oostep1: was lay
            leadsys  = 'bsmnim';
            mu       = 1.5e-4;
        end
        geom = selectleads(geom,useleads,1);
        if doatria
            geoma = selectleads(geoma,useleads,1);
            % initial scan
            disp(['atrial scan anisotropyratio: ' num2str(geoma.anisotropyratio)])
            [measinita.foci,measinita.dep,measinita.outp]=multifociscanatria(geoma,1,0);
            measinita.anisotropyratio = geoma.anisotropyratio;
            measinita.cor = measinita.outp(end,1);
            measinita.rd  = measinita.outp(end,2);
            measinita.rep=calcreptims(geoma,measinita.dep,0.00090);
            measinita.rep = measinita.rep - mean(measinita.rep) + 140;
        end
        disp(['ventricular scan anisotropyratio: ' num2str(geom.anisotropyRatio)]) % oostep1: was geom.anisotropyratio
        % [measinit.foci,measinit.dep,measinit.outp]=multifociscanpacing(geom,1,0);
        [measinit.foci,measinit.dep,measinit.outp]=multifociscan_wall(geom,'clusters',6,'usecor',1,'issinus',1,'velocity',0.8);

%         [measinit.foci,measinit.dep,measinit.outp]=multifociscanbrugada(geom,1,0);
        
        measinit.cor = measinit.outp(end,1);
        measinit.rd  = measinit.outp(end,2);
        measinit.rep = initrepppd(geom,measinit.dep,useleads);
        
        tmpmeas=inverse(geom,measinit.dep ,measinit.rep,'estimateampl',1,'casedir',dirout);
        measinit.amplitude = tmpmeas.amplfinal;
        
        %     measinit.rep=calcreptims(geom,measinit.dep,0.0079);
        
        %%
        pp=30;
        meas=inverse(geom,measinit.dep ,measinit.rep,'estimateampl',0,...
            'casedir',dirout,...
            'repopt','apd',...
            'maxiter',40,...
            'mudep',mu,...
            'murep',mu,...
            'weighed',0,...
            'minrd',0.15,...
            'leads',useleads,...
            'mode',4);
        %         figure(pp+10);showpatch(geom.ver,geom.itri,meas.depfinal,'nodes',measinit.foci);view(0,0);
        if length(useleads) ~= length(geom.v12)
            qtriplot('delete finalbsm')
            qtriplot(geom.ver,geom.itri,'finalbsm');
            qtriplot(meas.depfinal-min(meas.depfinal));
            qtriplot('step 10');
            qtriplot('panel 2,1');
            qtriplot('panel 2,1');
        else
            qtriplot('delete finalv12')
            qtriplot(geom.ver,geom.itri,'finalv12');
            qtriplot(meas.depfinal-min(meas.depfinal));
            qtriplot('step 10');
            qtriplot('panel 2,2');
            qtriplot('panel 2,2');
        end
        %         figure(pp+11);showpatch(geom.VER,geom.ITRI,meas.depfinal,'nodes',find(meas.depfinal==min(meas.depfinal)));view(0,0);
        %         figure(pp+20);showpatch(geom.ver,geom.itri,meas.repfinal,'nodes',measinit.foci);view(0,0);
        %         figure(pp+30);showpatch(geom.ver,geom.itri,meas.repfinal-meas.depfinal,'nodes',measinit.foci);view(0,0);
        %         figure(pp+31);clf;plot(meas.depfinal(geom.rfreewallver==0),meas.repfinal(geom.rfreewallver==0)-meas.depfinal(geom.rfreewallver==0),'.r');
        %         hold on;plot(meas.depfinal(geom.rfreewallver==1),meas.repfinal(geom.rfreewallver==1)-meas.depfinal(geom.rfreewallver==1),'.b')
        %
        t=0:geom.specs(5)-geom.specs(2);
        t=ones(length(geom.ver),1)*t;
        psia =lowpassma(geom.ama*getsmode(t,meas.depfinal,meas.repfinal,geom.ps,[],4),lpass);
        %         figure(pp+33);clf;plot(rms(geom.bsm(:,geom.specs(2):end)),'b');hold on;plot(rms(psia),'r')
        if length(useleads) ~= length(geom.v12)
            figure(pp+32);clf
            sigplot(geom.bsm(1:size(geom.lay,1)-1,geom.specs(2):geom.specs(5)),'',geom.lay,1,'b',0.25,0);
            hold on
            sigplot(psia(1:size(geom.lay,1)-1,:),'',geom.lay,1,'r',1,0);
            rrms=rms(psia);
            rrms=2*rrms/max(rrms)*geom.lay(1,2);
            t=(1:length(rrms))/length(rrms)*geom.lay(1,1);
            plot(t,rrms,'r')
            rrms=rms(baselinecor(geom.bsm(1:size(geom.lay,1)-1,geom.specs(2):geom.specs(5))));
            rrms=2*rrms/max(rrms)*geom.lay(1,2);
            t=(1:length(rrms))/length(rrms)*geom.lay(1,1);
            plot(t,rrms,'b')
        else
            figure(pp+33);clf
            l12 = [geom.v12(end-2:end) geom.v12(1:6)];
            ecg12 = baselinecor(geom.bsm([7:9 1:6],geom.specs(2):geom.specs(5)));
            sigplot(ecg12,'',loadmat('9lds.mla'),1,'b',1,0);
            hold on
            sigplot(psia([7:9 1:6],:),'',loadmat('9lds.mla'),1,'r',0,0);
            rrms=rms(psia);
            rrms=2*rrms/max(rrms)*3;
            t=(1:length(rrms))/length(rrms)*3;
            plot(t,rrms,'r')
            rrms=rms(ecg12);
            rrms=2*rrms/max(rrms)*3;
            t=(1:length(rrms))/length(rrms)*3;
            plot(t,rrms,'b')
        end
        qtriplot('funcolor tims');
        
        hold on
        if savecase
            meas.geom = geom;
            save([dirout geom.subject '_' leadsys '_' beat '_' type '_' date 'init.mat'],'measinit')
            save([dirout geom.subject '_' leadsys '_' beat '_' type '_' date '.mat'],'meas')
        end
        saveasci([dirout subject '_' leadsys '_' beat '_' type '_' date '.dep'],meas.depfinal);
        %%
        
        pp=300;
        meas2=inverse(geom,measinit.dep ,measinit.rep,'estimateampl',0,...
            'casedir',dirout,...
            'repopt','apd',...
            'estimatedamplitude', measinit.amplitude,...
            'maxiter',40,...
            'mudep',mu,...
            'murep',mu,...
            'weighed',0,...
            'minrd',0.15,...
            'leads',useleads,...
            'mode',4);
        %         figure(pp+10);showpatch(geom.ver,geom.itri,meas.depfinal,'nodes',measinit.foci);view(0,0);
        if length(useleads) ~= length(geom.v12)
            qtriplot('delete finalbsm2')
            qtriplot(geom.ver,geom.itri,'finalbsm2');
            qtriplot(meas2.depfinal-min(meas2.depfinal));
            qtriplot('step 10');
            qtriplot('horizontal 3');
            qtriplot('vertical 2');
            qtriplot('panel 3,1');
            qtriplot('panel 3,1');
        else
            qtriplot('delete finalv122')
            qtriplot(geom.ver,geom.itri,'finalv122');
            qtriplot(meas2.depfinal-min(meas2.depfinal));
            qtriplot('step 10');
            qtriplot('horizontal 3');
            qtriplot('vertical 2');
            qtriplot('panel 3,2');
            qtriplot('panel 3,2');
        end
        %         figure(pp+11);showpatch(geom.ver,geom.itri,meas2.depfinal,'nodes',find(meas2.depfinal==min(meas2.depfinal)));view(0,0);
        %         figure(pp+20);showpatch(geom.ver,geom.itri,meas2.repfinal,'nodes',meas2init.foci);view(0,0);
        %         figure(pp+30);showpatch(geom.ver,geom.itri,meas2.repfinal-meas2.depfinal,'nodes',meas2init.foci);view(0,0);
        %         figure(pp+31);clf;plot(meas2.depfinal(geom.rfreewallver==0),meas2.repfinal(geom.rfreewallver==0)-meas2.depfinal(geom.rfreewallver==0),'.r');
        %         hold on;plot(meas2.depfinal(geom.rfreewallver==1),meas2.repfinal(geom.rfreewallver==1)-meas2.depfinal(geom.rfreewallver==1),'.b')
        %
        t=0:geom.specs(5)-geom.specs(2);
        t=ones(length(geom.ver),1)*t;
        psia =lowpassma(geom.ama*getsmode(t,meas2.depfinal,meas2.repfinal,geom.ps,[],4),lpass);
        %         figure(pp+33);clf;plot(rms(geom.bsm(:,geom.specs(2):end)),'b');hold on;plot(rms(psia),'r')
        if length(useleads) ~= length(geom.v12)
            figure(pp+32);clf
            sigplot(geom.bsm(1:size(geom.lay,1)-1,geom.specs(2):geom.specs(5)),'',geom.lay,1,'b',0.25,0);
            hold on
            sigplot(psia(1:size(geom.lay,1)-1,:),'',geom.lay,1,'r',1,0);
            rrms=rms(psia);
            rrms=2*rrms/max(rrms)*geom.lay(1,2);
            t=(1:length(rrms))/length(rrms)*geom.lay(1,1);
            plot(t,rrms,'r')
            rrms=rms(baselinecor(geom.bsm(1:size(geom.lay,1)-1,geom.specs(2):geom.specs(5))));
            rrms=2*rrms/max(rrms)*geom.lay(1,2);
            t=(1:length(rrms))/length(rrms)*geom.lay(1,1);
            plot(t,rrms,'b')
        else
            figure(pp+33);clf
            l12 = [geom.v12(end-2:end) geom.v12(1:6)];
            ecg12 = baselinecor(geom.bsm([7:9 1:6],geom.specs(2):geom.specs(5)));
            sigplot(ecg12,'',loadmat('9lds.mla'),1,'b',1,0);
            hold on
            sigplot(psia([7:9 1:6],:),'',loadmat('9lds.mla'),1,'r',0,0);
            rrms=rms(psia);
            rrms=2*rrms/max(rrms)*3;
            t=(1:length(rrms))/length(rrms)*3;
            plot(t,rrms,'r')
            rrms=rms(ecg12);
            rrms=2*rrms/max(rrms)*3;
            t=(1:length(rrms))/length(rrms)*3;
            plot(t,rrms,'b')
        end
        qtriplot('funcolor tims');
        
        hold on
        if savecase
            meas2.geom = geom;
            save([dirout geom.subject '_2_' leadsys '_' beat '_' type '_' date '.mat'],'meas2')
        end
        saveasci([dirout subject '_2_' leadsys '_' beat '_' type '_' date '.dep'],meas.depfinal);
        
        %%
        if doatria
            pp=130;
            measa=inverse(geoma,measinita.dep ,measinita.rep,'estimateampl',0,...
                'casedir',dirout,...
                'repopt','apd',...
                'maxiter',40,...
                'mudep',1.5e-4,...
                'murep',1.5e-4,...
                'weighed',0,...
                'minrd',0.15,...
                'leads',useleads,...
                'mode',4);
            %             figure(pp+10);showpatch(geoma.ver,geoma.itri,measa.depfinal,'nodes',measinita.foci);view(0,0);
            
            %             figure(pp+11);showpatch(geoma.ver,geoma.itri,measa.depfinal,'nodes',find(measa.depfinal==min(measa.depfinal)));view(0,0);
            %             figure(pp+20);showpatch(geoma.ver,geoma.itri,measa.repfinal,'nodes',measinita.foci);view(0,0);
            %             figure(pp+30);showpatch(geoma.ver,geoma.itri,measa.repfinal-measa.depfinal,'nodes',measinita.foci);view(0,0);
            %             figure(pp+31);clf;plot(measa.depfinal,measa.repfinal-measa.depfinal,'.r');
            %             t=0:geoma.specs(3)-geoma.specs(2);
            %             t=ones(length(geoma.ver),1)*t;
            %             psia =lowpassma(geoma.ama*getsmode(t,measa.depfinal,measa.repfinal,geoma.ps,[],4),lpass);
            %             figure(pp+33);clf;plot(rms(geoma.bsm(:,geoma.specs(2):end)),'b');hold on;plot(rms(psia),'r')
            %             if k==0
            %                 figure(pp+32);clf
            %                 sigplot(geoma.bsm(1:size(geoma.lay,1)-1,geoma.specs(2):geoma.specs(3)),'',geoma.lay,1,'b',1,0);
            %                 hold on
            %                 sigplot(psia(1:size(geoma.lay,1)-1,:),'',geoma.lay,1,'r',1,0);
            %             else
            %
            %             end
            hold on
            if savecase
                meas.geom = geoma;
                save([dirout geoma.subject '_' leadsys '_' beat '_atria_' date 'init.mat'],'measinit')
                save([dirout geoma.subject '_' leadsys '_' beat '_atria_' date '.mat'],'meas')
            end
        end
    end
end