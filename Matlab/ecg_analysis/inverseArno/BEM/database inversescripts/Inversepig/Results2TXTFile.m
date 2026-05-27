display(['Starting ' mfilename]);
clearvars
ignoreavg=false; % true will calculate average of non-averaged beats only. Set to false for avg-only beat sets.
loadexisting=false%true; % When true will attempt to read the .mad files and only recreate the summary files.
plotact=true; % plot activation times vs distance to stim

% Will be set to true after processing with selstring==''!!!!
ellist={'LVEndo','RVEndo','Epi','Epi2','Epi3','Epi4','Epi5','Epi6','Epi7','Epi8','Epi9','Epi10'};
% selstring='' % If not empty, only select result files containing this string (Endo, Epi, LV, RV)
% selstring='Epi'
% SELSTRING={''};
SELSTRING={'','Endo','Epi'}; % make sure '' is first if loadexisting==false
% selstring='Endo'
% Rinit={};
% Rinv={};

sub =7;
for sub=[7 8 9 10];
ignorestring=''; % When filename contains this string file is excluded from selection

switch sub
    case 7
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/CORxcorrrmsshiftNewBSMNewGeom';
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/CORNoShiftNewBSMNewGeom';
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single';% controle PvD Peters Oude Beats
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/NewPeter20121022'; % opnieuw controle PvD. Eigen code useleads gefixed.
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/PO_OptimRatio1.0_20121022';
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/PvDHomogeenVarkenPig07_1.0_Vl1.5';
        % basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/homogeneous';
        % basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/NotInWallmudep1e-4murep1e-4';
        % basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/PO_Ratio1.0QRSScaleShiftInit10';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/Ingetest';
        basepath='';
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig07/Pig07_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig07/Pig07_ventricles.typ');
        bsmpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/AVG';
        subj='Pig07';
        ignorestring='LVAnteriorBaseEndo_20120627T181103_PNaT_beat000'; % double
        
    case 8
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig08/single/20130214';
%         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig08';
        
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig08/Pig08_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig08/Pig08_ventricles.typ');
        bsmpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/Pacing/AVG';
        subj='Pig08';
    case 9
        basepath='';%/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/Single_20130409_mudep1.5e-4_murep1.5e-9_FixedRepol_NotInWall';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/Single_20130425_mudep1e-5_murep1e-5_Fixedrepol_NotInWall';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ESC/Pig09_refined/ESC/single';
        
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.tri');
        %         [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09_refined/Pig09_refined_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.typ');
        %         surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09_refined/Pig09_refined_ventricles.typ');
        bsmpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG';
        ignorestring='pig09_009'; % When filename contains this string file is excluded from selection
        subj='Pig09';
    case 91
        basepath='';%/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single';
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09r/Pig09r_ventricles.typ');
        bsmpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG';
        sub=9
        subj='Pig09';
    case 10
        % basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/20130214';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/Single_20130409_mudep1.5e-4_murep1.5e-9_FixedRepol_NotInWall';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig10/single';
        
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig10/Pig10_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig10/Pig10_ventricles.typ');
        bsmpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/Biosemi/export/AVG';
        subj='Pig10';
end
beatpath=fullfile(bsmpath,'beats');

% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/test'%TEST
basepath=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results_b/' subj '/cluster1/minrd015_stuck_mode1_mudepscan'];
% basepath=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subj '/cluster1/HakkelGetS'];
% basepath=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subj '/cluster1/minrd015_mudep5e-5_blmode1'];
if ~exist('basepath','var') || isempty(basepath)
    basepath=uigetdir(['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subj '/cluster1'],sprintf('Select result folder to process for subject %d',sub));
end


basepath=regexprep(basepath,'\$|/$',''); % remove trailing / or \





D=dir(fullfile(basepath,'*.mat'));

if length(D)>=2
    BASEPATH={basepath};
else
    DD=dir(fullfile(basepath,'*'));
    isdir=cell2mat({DD.isdir});
    DD=DD(isdir);
    DD(1:2)=[]; % remove . and ..
    for i=1:length(DD)
        BASEPATH{i}=fullfile(basepath,DD(i).name);
    end
end
for selstringcell=SELSTRING
    selstring=cell2mat(selstringcell);
    % warning('inbp not starting at 1');
    for ibp=1:length(BASEPATH)
        clearvars -except ibp BASEPATH VERv ITRIv surftype bsmpath sub beatpath ignoreavg loadexisting ellist selstring ignorestring plotact SELSTRING
        fileinit={};
        fileinv={};
        cathpos=[];
        realsurf={};
        deltavstim=[];
        StimendQRS=[];
        QRSAmplitude=[];
        D=dir(fullfile(BASEPATH{ibp},'*.mat'));
        focusinv=[];
        for i=1:length(D)
            %     t=load(fullfile(basepath,D(i).name));
            fprintf('Pre-processing Folder %d of %d. File %d of %d\n',ibp,length(BASEPATH),i,length(D));
            if D(i).name(end-7:end)=='init.mat'
                t=load(fullfile(BASEPATH{ibp},D(i).name),'measinit');
                if exist('Rinit','var')
                    Rinit(end+1)=orderfields(t.measinit);
                else
                    Rinit=orderfields(t.measinit);
                end
                fileinit{end+1}=fullfile(BASEPATH{ibp},D(i).name);
                if isempty(cathpos)
                    [cathpos,realsurf, cathmap]=GetCatheterPosition2(sub,D(i).name); % first iteration of for loop
                else
                    [cathpos(end+1,:),realsurf(end+1), cathmap(end+1)]=GetCatheterPosition2(sub,D(i).name);
                end
                k=strfind(D(i).name,'_beat');
                beatfilename=fullfile(beatpath,[D(i).name(7:(k+7)) '.ecgspecs']);
                if ~exist(beatfilename,'file');
                    beatfilename=fullfile(beatpath,[D(i).name(7:(k+7)) '.spe']);
                    error('Spec file not found');%Old spec format, not checked');
                end
                ecgfilename=regexprep(beatfilename,'.ecgspecs$','.selecg');
                ecg=loadmat(ecgfilename);
                rrms=rms(zeromean(ecg));
                SPECS=loadmat(beatfilename);
                if isempty(SPECS) || length(SPECS)<13
                    deltavstim(end+1)=0;
                    StimendQRS(end+1)=0;
                    QRSAmplitude(end+1)=0;
                else
                    deltavstim(end+1)=SPECS(2)-SPECS(13); % Vstim to onsetQRS
                    StimendQRS(end+1)=SPECS(2)-SPECS(13)+SPECS(4);
                    QRSAmplitude(end+1)=max(max(abs(ecg(:,SPECS(2):SPECS(2)+SPECS(4)))));
                end
            else
                t=load(fullfile(BASEPATH{ibp},D(i).name),'meas');
                if ~isfield(t.meas,'log')
                    t.meas.log='';
                    t.meas.iterfinal=NaN;
                end
                if exist('Rinv','var')
                    Rinv(end+1)=orderfields(t.meas);
                else
                    Rinv=orderfields(t.meas);
                end
                fileinv{end+1}=fullfile(BASEPATH{ibp},D(i).name);
            end
            
        end
        % fileinv=fileinv';
        % fileinit=fileinit';
        
        if length(Rinv)~=length(Rinit)
            error('different number of initial and finale solutions');
        end
        
        for i=1:length(fileinit)
            if ~strncmp(fileinit{i},fileinv{i},length(fileinit{i})-8);
                error('Initial and final solution do not match');
            end
        end
        
        
        for i=1:length(Rinv)
            [~, focusinv(i,1)]=min(Rinv(i).depfinal);
            focusinvxyz(i,:)=VERv(focusinv(i),:);
            rdinv(i,1)=Rinv(i).rdfinal;
            corinv(i,1)=Rinv(i).corfinal;
        end
        
        for i=1:length(Rinit)
            focusinit(i,1)=Rinit(i).foci(1);
            focusinitxyz(i,:)=VERv(focusinit(i),:);
            rdinit(i,1)=Rinit(i).rd;
            corinit(i,1)=Rinit(i).cor;
        end
        
        
        
        %
        % M{1,1}=D(1).name;
        % M{2,1}=D(2).name;
        % M{:,2}={focusinit};
        % M{:,3}=focusinitxyz(:,1);
        % M{:,4}=focusinitxyz(:,2);
        % M{:,5}=focusinitxyz(:,3);
        % M{:,6}=focusinv;
        % M{:,7}=focusinvxyz(:,1);
        % M{:,8}=focusinvxyz(:,2);
        % M{:,9}=focusinvxyz(:,3);
        % dlmwrite('test.xls',M);
        
        % SUrtype 1-7: Epi, LVEndo,RVEndo,Mitral Valve, Tricuspid Valve,RVOT,Aorta
        surfnames={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % simplified
        
        %
        focusinitsurf=surfnames(surftype(focusinit))';
        focusinvsurf=surfnames(surftype(focusinv))';
        
        % oppsurfinit=unique(surfnames(~strcmp(focusinitsurf,surfnames)));
        % [pnearest, distp, trinearest, oppverinit, distver]=findnearest(oppsurfinit,VERv,ITRIv,surftype,oppsurfinit,1);
        % oppsurfinv=unique(surfnames(~strcmp(focusinvsurf,surfnames)));
        % [pnearest, distp, trinearest, oppverinv, distver]=findnearest(oppsurfinv,VERv,ITRIv,surftype,oppsurfinv,1);
        [~,foldername,ext]=fileparts(BASEPATH{ibp});
        if ~isempty(ext)
            foldername=[foldername ext]; % decimal dot mistaken for extension seperator.
        end
        
        if loadexisting && exist(fullfile(BASEPATH{ibp},[foldername '.mad']),'file')
            load(fullfile(BASEPATH{ibp},[foldername '.mad']),'-mat','C');
        else
            clear C
            l=length(Rinit);
            for i=1:l
                oppsurfinit=unique(surfnames(~strcmp(focusinitsurf(i),surfnames)));
                [pnearest, distp, trinearest, oppverinit, distver]=findnearest(VERv(focusinit(i),:),VERv,ITRIv,surftype,oppsurfinit,1);
                oppsurfinv=unique(surfnames(~strcmp(focusinvsurf(i),surfnames)));
                [pnearest, distp, trinearest, oppverinv, distver]=findnearest(VERv(focusinv(i),:),VERv,ITRIv,surftype,oppsurfinv,1);
                % TBD update cathpos.surfpos using .exactpos and VERv.
                
                
                
                % MEASUREMENT DATA
                C(i).filename=fileinv{i};
                C(i).Cathposx=cathpos(i).surfpos(1);
                C(i).Cathposy=cathpos(i).surfpos(2);
                C(i).Cathposz=cathpos(i).surfpos(3);
                C(i).Realsurf=realsurf{i};
                C(i).StimEndQRS=StimendQRS(i);
                C(i).StimOnsetQRS=deltavstim(i);
                C(i).OnsetEndQRS=StimendQRS(i)-deltavstim(i);
                C(i).QRSAmplitude=QRSAmplitude(i);
                
                % INITIAL ESTIMATE
                C(i).Initrd=Rinit(i).rd;
                C(i).Initcor=Rinit(i).cor;
                C(i).maxvelo=Rinit(i).outp(end,5);
                C(i).intialActTime=Rinit(i).outp(end,6);
                C(i).Initdep0=min(Rinit(i).dep);
                C(i).Initdepstim=deltavstim(i)+min(Rinit(i).dep);
                C(i).Initdispfocus=Rinit(i).dep(focusinit(i))-Rinit(i).dep(oppverinit);
                C(i).Initfocus=focusinit(i);
                C(i).Initx=focusinitxyz(i,1);
                C(i).Inity=focusinitxyz(i,2);
                C(i).Initz=focusinitxyz(i,3);
                C(i).Initsurf=focusinitsurf{i};
                C(i).Initdist=norm(focusinitxyz(i,:)-cathpos(i).surfpos);
                for j=1:length(ellist)
                    C(i).(['Initdact' ellist{j}])=NaN;
                end
                C(i).meanabsInitdact=NaN;
                C(i).dispInitdact=NaN;
                
                
                %OPTIMIZATION
                C(i).Invrd=Rinv(i).rdfinal;
                C(i).Invcor=Rinv(i).corfinal;
                C(i).Invdep0=min(Rinv(i).depfinal);
                C(i).Invdepstim=deltavstim(i)+min(Rinv(i).depfinal);
                C(i).Invdispfocus=Rinv(i).depfinal(focusinv(i))-Rinv(i).depfinal(oppverinv);
                C(i).Invfocus=focusinv(i);
                C(i).Invx=focusinvxyz(i,1);
                C(i).Invy=focusinvxyz(i,2);
                C(i).Invz=focusinvxyz(i,3);
                C(i).Invsurf=focusinvsurf{i};
                C(i).Invdist=norm(focusinvxyz(i,:)-cathpos(i).surfpos);
                
                %     ellist={'LVEndo','RVEndo','Epi','Epi2','Epi3','Epi4','Epi5','Epi6','Epi7','Epi8','Epi9','Epi10'};
                for j=1:length(ellist)
                    C(i).(['Invdact' ellist{j}])=NaN;
                end
                C(i).meanInvdact=NaN;
                C(i).dispInvdact=NaN;
                
                % ACTIVATION TIMES VS CATHETER
                DACTinit=[];
                DACTinv=[];
                invact=[];
                initact=[];
                dist=[];
                for iel=1:length(cathmap(i).dep)
                    [pnearest,~,trinearest,~,~,la,mu]=findnearest(cathmap(i).exactpos(iel,:),VERv,ITRIv,surftype,cathmap(i).surf{iel},1);
                    initact(iel)=trifun(Rinit(i).dep(ITRIv(trinearest,:)),la,mu)+deltavstim(i); % activation measured rel to stim
                    invact(iel)=trifun(Rinv(i).depfinal(ITRIv(trinearest,:)),la,mu)+deltavstim(i);
                    initdact=initact(iel)-cathmap(i).dep(iel);
                    invdact=invact(iel)-cathmap(i).dep(iel);
                    dist(iel)=norm(pnearest-cathpos(i).surfpos,'fro');
                    % vaste lijst met ellabels, scan op naam, zet waarde in veld,
                    % init op NaN???
                    k=find(strcmpi(cathmap(i).label{iel},ellist));
                    if ~isempty(k)
                        DACTinit(end+1)=initdact;
                        DACTinv(end+1)=invdact;
                        C(i).(['Initdact' ellist{k}])=initdact;
                        C(i).(['Invdact' ellist{k}])=invdact;
                    end
                end
                if plotact
                    labels={cathmap(i).label};
                    labels=labels{1};
                    iRVE=find(strncmpi('RVendo',labels,3));
                    iLVE=find(strncmpi('LVendo',labels,3));
                    iEndo=[iRVE iLVE];
                    iEpi=find(strncmpi('Epi',labels,3));
                    
                    hact=figure(31);
                    cla
                    hold on
                    legendstr={};
                    if ~isempty(iEndo)
                        plot(dist(iEndo),initact(iEndo),'+r');
                        plot(dist(iEndo),invact(iEndo),'*r');
                        plot(dist(iEndo),cathmap(i).dep(iEndo),'*b');
                        legendstr={'InitEndo','InvEndo','MeasEndo'};
                    end
                    if ~isempty(iEpi)
                        plot(dist(iEpi),initact(iEpi),'sr');
                        plot(dist(iEpi),invact(iEpi),'or');
                        plot(dist(iEpi),cathmap(i).dep(iEpi),'ob');
                        legendstr=[legendstr,{'InitEpi','InvEpi','MeasEpi'}];
                    end
                    [~,tstr]=fileparts(fileinv{i});
                    title(tstr,'Interpreter','none');
                    xlabel('Dist to stim[mm]');
                    ylabel('Act. time [ms]');
                    legend(legendstr,'Location','Best'); % ro
                    saveas(hact,regexprep(fileinv{i},'.mat$','.fig'),'fig');
                    saveas(hact,regexprep(fileinv{i},'.mat$','.jpg'),'jpg');
                end
                DACTinit(isnan(DACTinit))=[];
                DACTinv(isnan(DACTinv))=[];
                C(i).meanabsInitdact=mean(abs(DACTinit));
                if length(DACTinit)>1
                    C(i).dispInitdact=max(DACTinit)-min(DACTinit);
                else
                    C(i).dispInitdact=NaN;
                end
                C(i).meanabsInvdact=mean(abs(DACTinv));
                if length(DACTinv)>1
                    C(i).dispInvdact=max(DACTinv)-min(DACTinv);
                else
                    C(i).dispInvdact=NaN;
                end
                C(i).depfinal = Rinv(i).depfinal;
            end
            WriteStructsToText(fullfile(BASEPATH{ibp},[foldername '.txt']),C);
            save(fullfile(BASEPATH{ibp},[foldername '.mad']),'C');
        end
        % identify average beats, other beats per position
        
        fn={C.filename}';
        idrem=[];
        for i=1:length(fn)
            exportfn{i}=fn{i}(1:(strfind(fn{i},'_beat')-17));
            if (ignoreavg && ~isempty(strfind(fn{i},'_beat000'))) ||...
                    (~isempty(selstring) && isempty(strfind(fn{i},selstring)))||...
                    (~isempty(ignorestring) && ~isempty(strfind(fn{i},ignorestring)))
                idrem(end+1)=i;
            end
        end
        
        [un idx_last idx] = unique(exportfn(1,:));
        N = length(un);
        idx_unique = cell(1,N);
        for i = 1:N % for each unique element
            idx_unique{i} = find(strcmp(exportfn(1,:),un(i))==1);
            idx_unique{i}(ismember(idx_unique{i},idrem))=[]; % remove average beats
        end
        
        for i=1:length(idx_unique)
            if isempty(idx_unique{i})
                continue; % skip files for which we removed everything
            end
            avg(i).name=un{i};
            avg(i).Initrd=mean(cell2mat({Rinit(idx_unique{i}).rd}));
            %     avg(i).Initrdsd=std(cell2mat({Rinit(idx_unique{i}).rd}));
            avg(i).Initrdmin=min(cell2mat({Rinit(idx_unique{i}).rd}));
            avg(i).Initrdmax=max(cell2mat({Rinit(idx_unique{i}).rd}));
            
            avg(i).Initcor=mean(cell2mat({Rinit(idx_unique{i}).cor}));
            %     avg(i).Initcorsd=std(cell2mat({Rinit(idx_unique{i}).cor}));
            avg(i).Initcormin=min(cell2mat({Rinit(idx_unique{i}).cor}));
            avg(i).Initcormax=max(cell2mat({Rinit(idx_unique{i}).cor}));
            
            outpm=cell2mat({Rinit(idx_unique{i}).outp});
            avg(i).maxvelo=mean(outpm(5));
            avg(i).maxvelo=min(outpm(5));
            avg(i).maxvelo=max(outpm(5));
            
            avg(i).initialActTime=mean(outpm(6));
            avg(i).initialActTime=min(outpm(6));
            avg(i).initialActTime=max(outpm(6));
            
            
            avg(i).Initdep0=mean(cell2mat({C(idx_unique{i}).Initdep0}));
            avg(i).Initdep0sd=std(cell2mat({C(idx_unique{i}).Initdep0}));
            avg(i).Initdep0min=min(cell2mat({C(idx_unique{i}).Initdep0}));
            avg(i).Initdep0max=max(cell2mat({C(idx_unique{i}).Initdep0}));
            
            avg(i).Initdepstim=mean(cell2mat({C(idx_unique{i}).Initdepstim}));
            avg(i).Initdepstimabs=mean(abs(cell2mat({C(idx_unique{i}).Initdepstim})));
            
            avg(i).Initdispfocus=mean(cell2mat({C(idx_unique{i}).Initdispfocus}));
            avg(i).Initdispfocussd=std(cell2mat({C(idx_unique{i}).Initdispfocus}));
            avg(i).Initdispfocusmin=min(cell2mat({C(idx_unique{i}).Initdispfocus}));
            avg(i).Initdispfocusmax=max(cell2mat({C(idx_unique{i}).Initdispfocus}));
            
            
            avg(i).Realsurf=realsurf{idx_unique{i}(1)};
            avg(i).InitLREpi=sum(strcmp(realsurf(idx_unique{i}),focusinitsurf(idx_unique{i})'))/length(idx_unique{i});
            
            episcore=sum(strcmp('Epi',focusinitsurf(idx_unique{i})'))/length(idx_unique{i});
            if strfind(avg(i).Realsurf,'Epi')
                avg(i).InitEndoEpi=episcore;
            else
                avg(i).InitEndoEpi=1-episcore;
            end
            avg(i).Initdist=mean(cell2mat({C(idx_unique{i}).Initdist}));
            avg(i).Initdistsd=std(cell2mat({C(idx_unique{i}).Initdist}));
            avg(i).Initdistmin=min(cell2mat({C(idx_unique{i}).Initdist}));
            avg(i).Initdistmax=max(cell2mat({C(idx_unique{i}).Initdist}));
            
            for j=1:length(ellist)
                fn=['Initdact' ellist{j}];
                A=cell2mat({C(idx_unique{i}).(fn)});
                avg(i).(fn)=mean(A(~isnan(A)));
            end
            
            avg(i).meanmeanabsInitdact=nanmean(cell2mat({C(idx_unique{i}).meanabsInitdact}));
            avg(i).minmeanabsInitdact=max(cell2mat({C(idx_unique{i}).meanabsInitdact}));
            avg(i).maxmeanabsInitdact=min(cell2mat({C(idx_unique{i}).meanabsInitdact}));
            avg(i).meandispInitdact=nanmean(cell2mat({C(idx_unique{i}).dispInitdact}));
            avg(i).mindispInitdact=max(cell2mat({C(idx_unique{i}).dispInitdact}));
            avg(i).maxdispInitdact=min(cell2mat({C(idx_unique{i}).dispInitdact}));
            
            
            
            
            avg(i).Invrd=mean(cell2mat({Rinv(idx_unique{i}).rdfinal}));
            %     avg(i).Invrdsd=std(cell2mat({Rinv(idx_unique{i}).rd}));
            avg(i).Invrdmin=min(cell2mat({Rinv(idx_unique{i}).rdfinal}));
            avg(i).Invrdmax=max(cell2mat({Rinv(idx_unique{i}).rdfinal}));
            
            avg(i).Invcor=mean(cell2mat({Rinv(idx_unique{i}).corfinal}));
            %     avg(i).Invcorsd=std(cell2mat({Rinv(idx_unique{i}).cor}));
            avg(i).Invcormin=min(cell2mat({Rinv(idx_unique{i}).corfinal}));
            avg(i).Invcormax=max(cell2mat({Rinv(idx_unique{i}).corfinal}));
            
            avg(i).Invdep0=mean(cell2mat({C(idx_unique{i}).Invdep0}));
            avg(i).Invdep0sd=std(cell2mat({C(idx_unique{i}).Invdep0}));
            avg(i).Invdep0min=min(cell2mat({C(idx_unique{i}).Invdep0}));
            avg(i).Invdep0max=max(cell2mat({C(idx_unique{i}).Invdep0}));
            
            avg(i).Invdispfocus=mean(cell2mat({C(idx_unique{i}).Invdispfocus}));
            avg(i).Invdispfocussd=std(cell2mat({C(idx_unique{i}).Invdispfocus}));
            avg(i).Invdispfocusmin=min(cell2mat({C(idx_unique{i}).Invdispfocus}));
            avg(i).Invdispfocusmax=max(cell2mat({C(idx_unique{i}).Invdispfocus}));
            
            
            
            avg(i).InvLREpi=sum(strcmp(realsurf(idx_unique{i}),focusinvsurf(idx_unique{i})'))/length(idx_unique{i});
            
            episcore=sum(strcmp('Epi',focusinvsurf(idx_unique{i})'))/length(idx_unique{i});
            if strfind(avg(i).Realsurf,'Epi')
                avg(i).InvEndoEpi=episcore;
            else
                avg(i).InvEndoEpi=1-episcore;
            end
            avg(i).Invdist=mean(cell2mat({C(idx_unique{i}).Invdist}));
            avg(i).Invdistsd=std(cell2mat({C(idx_unique{i}).Invdist}));
            avg(i).Invdistmin=min(cell2mat({C(idx_unique{i}).Invdist}));
            avg(i).Invdistmax=max(cell2mat({C(idx_unique{i}).Invdist}));
            
            for j=1:length(ellist)
                fn=['Invdact' ellist{j}];
                A=cell2mat({C(idx_unique{i}).(fn)});
                avg(i).(fn)=mean(A(~isnan(A)));
            end
            
            avg(i).meanmeanabsInvdact=mean(cell2mat({C(idx_unique{i}).meanabsInvdact}));
            avg(i).maxmeanabsInvdact=min(cell2mat({C(idx_unique{i}).meanabsInvdact}));
            avg(i).minmeanabsInvdact=max(cell2mat({C(idx_unique{i}).meanabsInvdact}));
            avg(i).meandispInvdact=mean(cell2mat({C(idx_unique{i}).dispInvdact}));
            avg(i).maxdispInvdact=min(cell2mat({C(idx_unique{i}).dispInvdact}));
            avg(i).mindispInvdact=max(cell2mat({C(idx_unique{i}).dispInvdact}));
            
            
            
            A= cell2mat({C(idx_unique{i}).depfinal}) - mean(cell2mat({C(idx_unique{i}).depfinal}),2) * ones(1,length(idx_unique{i}));
            avg(i).stddep =  sqrt(mean(A(:).^2));
        end
        if ~isvarname('avg')
            continue % No results for this file, go to next.
        end
        
        nmeas=length(avg);
        imean=nmeas+1;
        avg(imean).name='MEAN';
        avg(imean).Initrd=nanmean(cell2mat({avg(1:nmeas).Initrd}));
        avg(imean).Initcor=nanmean(cell2mat({avg(1:nmeas).Initcor}));
        avg(imean).InitLREpi=nanmean(cell2mat({avg(1:nmeas).InitLREpi}));
        avg(imean).InitEndoEpi=nanmean(cell2mat({avg(1:nmeas).InitEndoEpi}));
        avg(imean).Initdist=nanmean(cell2mat({avg(1:nmeas).Initdist}));
        avg(imean).Initdep0=nanmean(cell2mat({avg(1:nmeas).Initdep0}));
        avg(imean).Initdepstim=nanmean(cell2mat({avg(1:nmeas).Initdepstim}));
        avg(imean).Initdepstimabs=nanmean(cell2mat({avg(1:nmeas).Initdepstimabs}));
        
        istd=imean+1;
        avg(istd).name='STD';
        avg(istd).Initrd=std(cell2mat({avg(1:nmeas).Initrd}));
        avg(istd).Initcor=std(cell2mat({avg(1:nmeas).Initcor}));
        avg(istd).InitLREpi=std(cell2mat({avg(1:nmeas).InitLREpi}));
        avg(istd).InitEndoEpi=std(cell2mat({avg(1:nmeas).InitEndoEpi}));
        avg(istd).Initdist=std(cell2mat({avg(1:nmeas).Initdist}));
        avg(istd).Initdep0=std(cell2mat({avg(1:nmeas).Initdep0}));
        avg(istd).Initdepstim=std(cell2mat({avg(1:nmeas).Initdepstim}));
        avg(istd).Initdepstimabs=std(cell2mat({avg(1:nmeas).Initdepstimabs}));
        
        avg(imean).Invrd=nanmean(cell2mat({avg(1:nmeas).Invrd}));
        avg(imean).Invcor=nanmean(cell2mat({avg(1:nmeas).Invcor}));
        avg(imean).InvLREpi=nanmean(cell2mat({avg(1:nmeas).InvLREpi}));
        avg(imean).InvEndoEpi=nanmean(cell2mat({avg(1:nmeas).InvEndoEpi}));
        avg(imean).Invdist=nanmean(cell2mat({avg(1:nmeas).Invdist}));
        avg(imean).meanmeanabsInitdact=nanmean(cell2mat({avg(1:nmeas).meanmeanabsInitdact}));
        avg(imean).minmeanabsInitdact=min(cell2mat({avg(1:nmeas).minmeanabsInitdact}));
        avg(imean).maxmeanabsInitdact=max(cell2mat({avg(1:nmeas).maxmeanabsInitdact}));
        avg(imean).meanmeanabsInvdact=nanmean(cell2mat({avg(1:nmeas).meanmeanabsInvdact}));
        avg(imean).minmeanInvitdact=min(cell2mat({avg(1:nmeas).minmeanabsInvdact}));
        avg(imean).maxmeanInvitdact=max(cell2mat({avg(1:nmeas).maxmeanabsInvdact}));
        avg(imean).meandispInitdact=nanmean(cell2mat({avg(1:nmeas).meandispInitdact}));
        avg(imean).mindispInitdact=min(cell2mat({avg(1:nmeas).mindispInitdact}));
        avg(imean).maxdispInitdact=max(cell2mat({avg(1:nmeas).maxdispInitdact}));
        avg(imean).meandispInvdact=nanmean(cell2mat({avg(1:nmeas).meandispInvdact}));
        avg(imean).mindispInvitdact=min(cell2mat({avg(1:nmeas).mindispInvdact}));
        avg(imean).maxdispInvitdact=max(cell2mat({avg(1:nmeas).maxdispInvdact}));
        
        avg(istd).Invrd=std(cell2mat({avg(1:nmeas).Invrd}));
        avg(istd).Invcor=std(cell2mat({avg(1:nmeas).Invcor}));
        avg(istd).InvLREpi=std(cell2mat({avg(1:nmeas).InvLREpi}));
        avg(istd).InvEndoEpi=std(cell2mat({avg(1:nmeas).InvEndoEpi}));
        avg(istd).Invdist=std(cell2mat({avg(1:nmeas).Invdist}));
        avg(istd).meanmeanabsInitdact=std1d(cell2mat({avg(1:nmeas).meanmeanabsInitdact}));
        avg(istd).minmeanabsInitdact=std1d(cell2mat({avg(1:nmeas).minmeanabsInitdact}));
        avg(istd).maxmeanabsInitdact=std1d(cell2mat({avg(1:nmeas).maxmeanabsInitdact}));
        avg(istd).meanmeanabsInvdact=std1d(cell2mat({avg(1:nmeas).meanmeanabsInvdact}));
        avg(istd).minmeanabsInvdact=std1d(cell2mat({avg(1:nmeas).minmeanabsInvdact}));
        avg(istd).maxmeanabsInvdact=std1d(cell2mat({avg(1:nmeas).maxmeanabsInvdact}));
        
        avg(istd).meandispInitdact=std1d(cell2mat({avg(1:nmeas).meandispInitdact}));
        avg(istd).mindispInitdact=std1d(cell2mat({avg(1:nmeas).mindispInitdact}));
        avg(istd).maxdispInitdact=std1d(cell2mat({avg(1:nmeas).maxdispInitdact}));
        avg(istd).meandispInvdact=std1d(cell2mat({avg(1:nmeas).meandispInvdact}));
        avg(istd).mindispInvitdact=std1d(cell2mat({avg(1:nmeas).mindispInvdact}));
        avg(istd).maxdispInvitdact=std1d(cell2mat({avg(1:nmeas).maxdispInvdact}));
        
        
        for j=1:length(ellist)
            fn=['Initdact' ellist{j}];
            avg(imean).(fn)=nanmean(abs(cell2mat({avg(1:nmeas).(fn)})));
            avg(istd).(fn)=std1d(abs(cell2mat({avg(1:nmeas).(fn)})));
            fn=['Invdact' ellist{j}];
            avg(imean).(fn)=nanmean(abs(cell2mat({avg(1:nmeas).(fn)})));
            avg(istd).(fn)=std1d(abs(cell2mat({avg(1:nmeas).(fn)})));
        end
        
        
        [~,foldername,ext]=fileparts(BASEPATH{ibp});
        if ~isempty(ext)
            foldername=[foldername ext]; % decimal dot mistaken for extension seperator.
        end
        WriteStructsToText(fullfile(BASEPATH{ibp},[foldername selstring '_Summary.txt']),avg);
        save(fullfile(BASEPATH{ibp},[foldername selstring '_Summary.mad']),'avg');
    end
    if isempty(selstring)
        loadexisting=true; % all non-summary files now up to date.
    end
end
end
