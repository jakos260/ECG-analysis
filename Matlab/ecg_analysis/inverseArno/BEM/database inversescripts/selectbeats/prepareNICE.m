%% This is a version of preparePigs.m of PO [20150219], specific for PVC AMS STUDY
%
% function prepareNICE(sub, namefile)
%
% INPUT:
% sub       = NICE subject number
% namefile  = matlab file name [including extension .mat]
%
% Calls prepareECG.m
%
% Output files ...
%
% Arno Janssen 20150225: v 1.0

function prepareNICE(sub, namefile)

if nargin < 1, error('no subject number given'); end
if nargin < 2, namefile = ''; end

clear global lpass
maxbeatlength   = 1500;     % ?

% main path for PVC AMS STUDY:
pvcpath         = '/Users/arnojanssen/Documents/STW/PVCs/';

beatlim         = 5;            % limit number of beats per file
skipexisting    = 0;            % when 1 will skip if file already exists
ibeat_start     = 1;            % if 0, also do average

pwavecorflag    = false;        % true will use P-wave corrected version
averagebeatflag = false;        % when true take average beat as beat 0.
removehum       = false;        % use subtracthum signal, for now only if pwavecor and averagebeatflag os false

saveCase        = 2;
layfile         = 'ams65.mla';  % default layout file
type            = 'ventricles';

subject         = ['NICE00' num2str(sub)];
bsmdir          = [pvcpath subject '/Biosemi/export/'];
dirout          = [pvcpath subject '/results'];

% check if directory exists, else create
if exist(dirout, 'dir') == 7, else mkdir([pvcpath subject], 'results'); end

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
            
            GEOM.LAY            = [GEOM.LAY(1,:); L];
            GEOM.BSM            = zeromean(GEOM.BSM(useLeads,:)); % Instead of selectLeads.m
            GEOM.subject        = subject;
            GEOM.SPECS          = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim','0');            
        end
        
        % count new and existing beat, but not the average beat
        if ibeat > 0 && DATA.SELBEATS(ibeat),
            beatcount = beatcount+1;
            if beatcount>= beatlim && find(DATA.SELBEATS(ibeat:end),1,'last')>ibeat, fprintf('Number of beats limited on %d\n',beatlim); break; end
        end
    end
    
    display(beatcount);
    
end

close all;