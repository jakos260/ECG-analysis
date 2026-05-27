%% Script prepare.m 
%
% This script is used after ViewandExport to define the markers: 
% QRS onset, QRS end, Peak T-wave, etc...
% 
% The settings for the script are specified in the first part called
% 'INPUT'.
%
% The output is saved in a new folder 'beats',
% with every single beat as a seperate file with extensions '.selecg' &
% '.ecgspecs'.
%
% Version 1: 01-apr-2015

%% INPUT:

clear global lpass
maxbeatlength   = 1500;     % ?

subject         = 'NICE001';                            % name folder subject
namefile        = '';                                   % define filename. if namefile = '' --> all files with extension .mat are included.
pvcpath         = inversepath(1);                       % main folder study data defined in inversepath.m

bsmdir          = [pvcpath subject '/Biosemi/export/']; % define folder with matlab-file(s) produced with ViewandExport.m    
dirout          = [pvcpath subject '/results'];         % define folder to save results

% check if directory exists, else create
if exist(dirout, 'dir') == 7, else mkdir([pvcpath subject], 'results'); end

layfile         = 'ams65.mla';                          % default layout file
type            = 'ventricles';

beatlim         = 5;                                    % limit number of beats per file
skipexisting    = 0;                                    % when 1 will skip if file already exists
ibeat_start     = 1;                                    % if 0, also do average

pwavecorflag    = false;                                % true will use P-wave corrected version
averagebeatflag = false;                                % when true take average beat as beat 0.
removehum       = false;                                % use subtracthum signal, for now only if pwavecor and averagebeatflag os false

saveCase        = 2;

%%
% Check if average beat should be used: stored in subfolder AVG
if averagebeatflag || pwavecorflag, bsmdir = fullfile(bsmdir,'AVG'); end

if strcmp(namefile, ''), bsmfiles = dir(fullfile(bsmdir, '*.mat')); else bsmfiles = dir(fullfile(bsmdir, namefile)); end

for k = 1:length(bsmfiles),
    beatcount = 0;
    load(fullfile(bsmdir, bsmfiles(k).name));
    
    for ibeat = ibeat_start:length(DATA.SELBEATS),
        fn = fullfile(bsmdir ,'beats' ,[bsmfiles(k).name(1:end-4) '_beat' num2str(ibeat,'%03d') '.selecg']); % name .selecg
        
        if ((ibeat==0 && averagebeatflag) || (ibeat>0 && DATA.SELBEATS(ibeat))) && (skipexisting==0 || ~exist(regexprep(fn,'.selecg$','.ecgspecs'),'file'))
            if ~exist(dirout,'dir'),                    mkdir(dirout);                      end
            if ~exist(fullfile(bsmdir,'beats'),'dir'),  mkdir(fullfile(bsmdir,'beats'));    end

            if pwavecorflag && isfield(DATA,'BSMOUTPWC'),
                if ibeat, savemat(fn, DATA.BSMOUTPWC(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength)));
                    GEOM.BSM = DATA.BSMOUTPWC(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength));
                elseif averagebeatflag
                    savemat(fn,DATA.AVERAGEPWC);
                    GEOM.BSM = DATA.AVERAGEPWC;
                else
                    continue
                end
                
            else
                if pwavecorflag, warning('No P-wave corrected version available'); end
                if ibeat, 
                    if removehum,
                        if isfield(DATA,'BSMOUTSUB50Hz'),
                            savemat(fn, DATA.BSMOUTSUB50Hz(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength)));
                            GEOM.BSM = DATA.BSMOUTSUB50Hz(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength));
                        else
                            savemat(fn, subtracthum(DATA.BSMOUT(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength))));
                            GEOM.BSM    = subtracthum(DATA.BSMOUT(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength)));
                            fht         = figure('Name','Original signal (red) and after humm subtraction (black');
                            plot(rms(DATA.BSMOUT(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength))),'r');
                            hold on; 
                            plot(std(GEOM.BSM),'k');
                            pause
                            close(fht);
                        end
                        
                    else
                        savemat(fn, DATA.BSMOUT(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength)));
                        GEOM.BSM = DATA.BSMOUT(:,DATA.BEATS(ibeat):min(DATA.BEATS(ibeat+1),DATA.BEATS(ibeat)+maxbeatlength));
                    end
                elseif averagebeatflag
                    savemat(fn,DATA.AVERAGE);
                    GEOM.BSM = DATA.AVERAGE;
                else
                    continue
                end
            end
            
            bsmfile = fn;
            
            % Load all data from selected subject
            disp('===============================================================')
            disp(['Loading the ' type  ' of subject ' subject '  selected beat:' num2str(ibeat) ' from: ' bsmfiles(k).name])
            
            GEOM.LAY                            = loadmat(layfile);
            remove                              = DATA.remove;
            L                                   = GEOM.LAY(2:end,:);
            useLeads                            = find(remove == 0);
            L(ismember(L(:,3),find(remove)),:)  = [];
            for i = 1:size(L,1),    L(i,3)      = L(i,3) - sum(remove (1:L(i,3)));  end
            
            GEOM.LAY                            = [GEOM.LAY(1,:); L];
            GEOM.BSM                            = zeromean(GEOM.BSM(useLeads,:)); % Instead of selectLeads.m
            GEOM.subject                        = subject;
            GEOM.SPECS                          = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim','0');            
        end
        
        % count new and existing beat, but not the average beat
        if ibeat > 0 && DATA.SELBEATS(ibeat),
            beatcount = beatcount+1;
            if beatcount >=  beatlim && find(DATA.SELBEATS(ibeat:end),1,'last') > ibeat, fprintf('Number of beats limited on %d\n',beatlim); break; end
        end
    end
    
    display(beatcount);
    
end

close all;