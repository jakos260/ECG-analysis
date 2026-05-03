%% This is a version of ViewandExport.m of PO [20150219]
%
% function sel_beat_aj(filepath, ElecPos, defsel, defstart, defend)
%
% Beats can be selected and saved to export file as matlab structure.
%
% INPUT:
% filepath:     Name of path where ECG-files are stored [inlcuding file name and extension]
% ElecPos:      Name of electrode position file. For example 'amsbdf'
% defsel:       Vector with numbers of electrodes to evaluate
% defstart:     Define begin of time signal
% defend:       Define end of time signal as: + # ms compared to defstart
%
% Output files now constructed for 65 channels.
%
% Arno Janssen 20150219: v 1.0

function sel_beat_aj(filepath,ElecPos,defsel,defstart,defend)

%% Define begin parameters and ECG-filename and path
close all;

fl_sigplot      = 1; % plot signal
fl_selectbeats  = 1; % call Peters selectbeats
fl_onlyinmla    = false;
%subjcode        = 'NICE001';

if ~exist('ElecPos'), ElecPos = 'amsbdf'; end
if ~exist('filepath'), defaultpath      = ['/Users/arnojanssen/Documents/STW/']; % defaultpath      = ['/Users/arnojanssen/Documents/STW/PVCs/' subjcode '/Biosemi/'];
    [FileName, PathName, FilterIndex]   = uigetfile({'*.bdf;*.zip;*.XML;*.xml;*.csv;*.cipsecg'},'Select data file.',defaultpath);
    filepath                            = fullfile(PathName,FileName);
end

[folder, filename, ext] = fileparts(filepath);

%% Load raw data ECG-file
if strcmpi(ext,'.bdf'), 
    if exist('defsel'),     DefAns{1} = defsel;     else DefAns{1} = '[1:62 257 258 259]';    end % '[1:64 259 257 258]' based on 67 channels!
    if exist('defstart'),   DefAns{2} = defstart;   else DefAns{2} = '0';       end
    if exist('defend'),     DefAns{3} = defend;     else DefAns{3} = '+20';     end
   
    if nargin > 0, D = Readbdf(filepath,DefAns{1},DefAns{2},DefAns{3},'noui'); else D = Readbdf(filepath,DefAns{1},DefAns{2},DefAns{3}); end
    D.filepath = filepath;
else
    error('unknown file type');
end

%% Plot example signal
if fl_sigplot,
    XX = D.startpos+(1:size(D.signals,1))*D.Ts/1000;
    figure
    xlim([D.startpos,D.startpos+10]);
    plot(XX,D.signals(:,2)); % plot an example signal: channelnumber 2
    drawnow
    figure % dummy figure
end

%% Actual selection of beats and save to export folder.
if fl_selectbeats,
    % global Nchannels;
    % Nchannels = 60; % of toch 64?
    % psm; % maplab BSM layout
    
    [mla, standardc,leadsys] = elecpos2mlas(D.SelectChannels,ElecPos); % set ElecPos at start of this m-file.
    % mla seq L-R, T-D. org bdf chan numbers. standardc vr,vl,vf (Nijmegen), org bdf chan numbers.
    mlas        = mla; % mla for selected channels only. First line with grid size added later
    delchan     = [];
    standard    = zeros(1,9);
    
    % move grid so it begins in upper left corner [1,1]
    if mlas(1,3) == 0, mlas = mlas(2:end,:); end          % Remove grid size if present.
    mlas(:,1)               = mlas(:,1)-min(mlas(:,1))+1; % move grid to upper boundary
    mlas(:,2)               = mlas(:,2)-min(mlas(:,2))+1; % move grid to left boundary
    mlas                    = [max(mlas(:,1:2)),0; mlas]; % Now add first line with grid size
    
    sig1000 = interp1(1:size(D.signals,1),D.signals,1:1/D.Ts:size(D.signals,1),'pchip'); % D.Ts in ms, to Ts = 0.001 s RESAMPLE!
    
    waitfor(selectBSMbeat(sig1000','lay',mlas,'sampt',1/1000,'leadsys',leadsys));
    
    global DATA
    
    % save data to export folder:
    response = inputdlg({'Short description (empty aborts save)','Comment'});
    if ~isempty(response) && ~isempty(response{1}),
        D.comment   = response{2};
        folderexp   = [folder filesep 'export'];
        x           = exist(folderexp,'file');
        
        if x == 2, folderexp = folder; elseif x ~= 7, mkdir(folderexp); end % there is a folder named 'export' 
        save(fullfile(folderexp, [filename '_' num2str(D.startpos) '_' num2str(D.endpos) '_'  response{1} '_' datestr(now,30) '.mat']), 'DATA', 'D', 'mlas');
    end
end

close all