clearvars
savepng=false; 
sub=9;
switch sub
    case 7
        basepath='';
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig07/Pig07_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig07/Pig07_ventricles.typ');
    case 8
        basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ESC/Pig08';
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig08/Pig08_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig08/Pig08_ventricles.typ');
    case 9
        basepath='';
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.tri');
        %         [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09_refined/Pig09_refined_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.typ');
        %         surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09_refined/Pig09_refined_ventricles.typ');
        bsmpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG';
    case 10
        % basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/20130214';
        %         basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/Single_20130409_mudep1.5e-4_murep1.5e-9_FixedRepol_NotInWall';
        basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig10/single';
        
        [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig10/Pig10_ventricles.tri');
        surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig10/Pig10_ventricles.typ');
end
beatpath=fullfile(bsmpath,'beats');



% TEST
% PathName='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131111/20130918_im1_wrd1_iV0.4_mV1.5Anis0.1_Twospeed';
% PathName='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single.old/scaledPvD/multiple_initialAnis2.5/';
% FileName='Pig09_pig09_005_155_184_RVApexPostEpiThrx1SyncVoff_20130211T133347_beat000_ventricles_1_0.01_02-Oct-2013.mat';

% 20141208 Manuscript
PathName='/Users/peteroosterhoff/Documents/Werk/STW/Data/results_b/Pig09/cluster1/minrd015_stuck_mode1_mudepscan/20141117_im0_wrd0_iV0.4_iAnis1.0Anis0.01mudep5.00e-06';
FileName='Pig09_pig09_007_303_332_LVLatEndoThrx1SyncVoff_20130211T141135_beat000_ventricles_im0_wrd0_iV0.4_iAnis1.0Anis0.01_16-Nov-2014.mat';


% FileName='Pig09_pig09_005_155_184_RVApexPostEpiThrx1SyncVoff_20130211T133347_beat000_ventricles_1_0.1_18-Sep-2013.mat';
% %TEST
% D=dir(fullfile(PathName,'Pig09_pig09_005_155_184_RVApexPostEpiThrx1SyncVoff_20130211T133347_beat000_ventricles_1_0.1_11-Nov-2013.mat'));
% FileName={D.name};


if ~exist('FileName','var') || isempty(FileName)
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select  result files'...
        ,'/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single','MultiSelect','on');
end
if ~iscell(FileName)
    FileName={FileName};
end
for ifile=1:length(FileName);
    qtriplot('reset');
    FilePath=(fullfile(PathName,FileName{ifile}));
    ta=load(FilePath);
    if isfield(ta,'meas')
        dep=ta.meas.depfinal;
    elseif isfield(ta,'measinit')
        dep=ta.measinit.dep;
    else
        error('No depolarization times found');
    end
    [cathpos,realsurf, cathmap]=GetCatheterPosition2(sub,FileName{ifile});
    k=strfind(FileName{ifile},'_beat');
    beatfilename=fullfile(beatpath,[FileName{ifile}(7:(k+7)) '.ecgspecs']);
    if ~exist(beatfilename,'file');
        beatfilename=fullfile(beatpath,[FileName{ifile}(7:(k+7)) '.spe']);
        error('Old spec format, not jet implemented/checked');
    end
    ecgfilename=regexprep(beatfilename,'.ecgspecs$','.selecg');
    ecg=loadmat(ecgfilename);
    rrms=rms(zeromean(ecg));
    SPECS=loadmat(beatfilename);
    if isempty(SPECS) || length(SPECS)<13
        warning('No t=0 defined (sinus beat?)');
        deltavstim=0;
        StimendQRS=0;
        QRSAmplitude=0;
    else
        deltavstim=SPECS(2)-SPECS(13); % Vstim to onsetQRS
        StimendQRS=SPECS(2)-SPECS(13)+SPECS(4);
        QRSAmplitude=max(max(abs(ecg(:,SPECS(2):SPECS(2)+SPECS(4)))));
    end
    
    for icath=1:3
        switch icath
            case 1
                ellist={'Epi','Epi2','Epi3','Epi4','Epi5','Epi6','Epi7','Epi8','Epi9','Epi10'};
                cathname='EpiCath';
            case 2
                ellist={'RVEndo'};
                cathname='RVEndoCath';
            case 3
                ellist={'LVEndo'};
                cathname='LVEndoCath';
        end
        cathactt=NaN(1,length(ellist));
        actt=NaN(1,length(ellist));
        cathpos=NaN(length(ellist),3);
        
        for iel=1:length(cathmap.dep)
            [pnearest,~,trinearest,~,~,la,mu]=findnearest(cathmap.exactpos(iel,:),VERv,ITRIv,surftype,cathmap.surf{iel},1);
            ta=trifun(dep(ITRIv(trinearest,:)),la,mu)+deltavstim;
            dta=ta-cathmap.dep(iel);
            % vaste lijst met ellabels, scan op naam, zet waarde in veld,
            % init op NaN???
            k=find(strcmpi(cathmap.label{iel},ellist));
            if ~isempty(k)
                actt(k)=ta;
                cathpos(k,:)=pnearest;%cathmap.exactpos(iel,:);
                cathactt(k)=cathmap.dep(iel);
            end
        end
        idel=any(isnan(cathpos),2); % ellist without position
        cathpos(idel,:)=[];
        actt(idel)=[];
        cathactt(idel)=[];
        
        qtriplot(VERv,ITRIv,'Ventricles');
        qtriplot(dep,'Ventricles');
        qtriplot('funscale 0 125');
        qtriplot('step 5');
        qtriplot('funcolor tim');
        if size(cathpos,1)>1
            firstl=find(~isnan(cathact,1,'first'));
            lastl=find(~isnan(cathact,1,'last'));
            cathposf=cathpos(firstl:lastl,:);
            cathacf=cathact(firstl:lastl);
            id=[0,3,6,9,15,18,24,27,33,36]; % standard eldis
            idf=id(firstl:lastl); % where function valid
            idfi=[idf(1):1.5,idf(end)]; %
            
            % whole catheter
            [VERcc,ITRIcc]=make_snake(cathpos,id,12,0.9); % catheter contour
            qtriplot(VERcc,ITRIcc,[cathname '_c']);
            qtriplot('transparent 0.5');
            qtriplot('edge on');

            
            % part of catheter with activation times
            [VERca,ITRIca,intcathpos]=make_snake(cathposf,idfi,12,1); % interpolated catheter
            
            qtriplot(VERca,ITRIca,[cathname '_a']);
            %             if ~isempty(actti)
            %                 qtriplot(actti',[cathname '_a']);
            %             end
            qtriplot('back *=both');
        else
        end
        
    end
    
    if length(FileName)>1 || savepng
        fprintf('%s\n',FileName{ifile});
        pause
        if savepng
            pngdir=fullfile(PathName,'png');
            if ~exist(pngdir,'dir')
                mkdir(pngdir);
            end
            qtriplot(sprintf('png %s',fullfile(pngdir,FileName{ifile}(1:end-4))));
            display('Saved');
        end
    end
end




