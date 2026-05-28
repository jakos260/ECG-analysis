% function [EGM, CARTO] = readmpt(fpath, subject)
%
% This function is used to convert the electrograms recorded with CARTO
% into a matlab format. All .mpt files within the subject directory are included.
% The data should be in the directory: ['/CARTO/Study/Patients/S1/PointsCommonData']
% Each .mpt file corresponds to a location on the heart.
%
% Also the data from SQL [save in a .xls format] are put in a matlab
% structure. The data should be in the directory: ['/exportSQL/' subject '.xls']
%
% INPUT:
% fpath     = Directory name of subject
% subject   = Subject name including quotes [ex. 'NICE001']
%
% OUTPUT:
% EGM       = Structure with all electrograms, filed per .mpt file.
% CARTO     = Structure with data from database file [SQL].
%
% NOTE: The signal names can be found in the SQL database tables: POINT_TABLE & CONFIG_CHANNEL_TABLE:
%
% ch_name = {'UNI-MAP-1'; 'UNI-MAP-2'; 'UNI-MAP-3'; 'UNI-MAP-4'; 'V1'; 'V2'; 'V3'; 'V4'; 'V5'; 'V6'; ...
%    'BI-M1-M2'; 'BI-M3-M4'; 'I'; 'II'; 'III'; 'aVL'; 'aVR'; 'aVF'};
% 
% The reference channel is defined in REF_ANNOTATION_MAP_SETUP_TABLE & CONFIG_CHANNEL_TABLE
%
% Version 1: 01-apr-2015

function [EGM, CARTO] = readmpt(fpath, subject)

%% Set parameters and names:

% If directory name is not defined:
if nargin < 1, subject = 'NICE001'; fpath = ['/Users/arnojanssen/Documents/STW/PVCs/' subject]; end

% Diretory path names extended:
mptpath             = [fpath '/CARTO/Study/Patients/S1/PointsCommonData'];
sqlpath             = [fpath '/exportSQL/' subject '.xls'];

% Load data from excell file:
CARTO.points        = xlsread(sqlpath, 'points');           % Point indices and index for .mpt files
CARTO.coordinates   = xlsread(sqlpath, 'coordinates');      % Includes 3D coordinates of points
CARTO.map_anno      = xlsread(sqlpath, 'map_annotation');   % Includes latencies [activation times]
CARTO.markers       = xlsread(sqlpath, 'multi_tag_table');  % Includes the markers that corrspond to the x-ray images

% for all .mpt files in the directory:
files               = dir(fullfile(mptpath,'*.mpt'));
numvec              = [];

% Get file names:
for k = 1:size(files,1), numvec(k) = str2num(files(k).name(2:end-4)); end

% Sort the files by number:
numvec = sort(numvec);

% Create structure for EGM data:
EGM                 = [];
EGM.sigs            = [];
EGM.filename        = [];

% Define standard settings:
nchline             = 33;
nsline              = 31;
headerlines         = 37;
endtrial            = 31;

%% Get data and file in EGM structure:
for k = 1:size(numvec,2),
    
    filename    = [mptpath '/P' num2str(numvec(k)) '.mpt'];
    fid         = fopen(filename);
    A           = fread(fid,'short');
    nch         = A(nchline);
    ns          = A(nsline);
    sigs        = zeros(nch,ns);
    
    % Forr all channels defined in header mpt-file:
    for m = 1:nch,
        vecsig(m,:) = headerlines + (1:ns) + (ns+endtrial)*(m-1);
        sigs(m,:)   = A(vecsig(m,:));
    end
    
    % File data:
    EGM(k).sigs     = sigs;
    EGM(k).filename = filename;
end

% save data to file:
if ~exist([fpath '/CATHDATA'],'dir'), mkdir([fpath '/CATHDATA']); end
save([fpath '/CATHDATA/CARTO.mat'], 'EGM', 'CARTO');