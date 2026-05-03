% function [stimpos, cathsurf, cathmap, bsm] = GetCatheterPosition2(sub,filename)
%
% This function is used to load the stimulation information from a text
% file, defined per subject in cathpos.m
%
% INPUT:
% sub       = Subject number.
% filebname = Name of beat-file.
%
% OUTPUT: 
% stimpos   = The stimulation position.
% cathsurf  = surface of stimulation?
% cathmap   = The catheter positions and activation times rel. to stim.
% bsm       = The body surface map.
%
% markers per orgname
% 20130821 oostep1 updated from Pig09_Catherpositions, added LVApex
% 20140319 oostep1 added Pig10
%
% Version 1: 01-apr-2015

function [stimpos, cathsurf, cathmap, bsm] = GetCatheterPosition2(sub,filename)

%% Set main variables:
surface                 = {};
bsm                     = [];
mfilepath               = mfilename('fullpath');
[datapath,dum1,dum2]    = fileparts(mfilepath);
stimpos                 = [];
cathsurf                = '';
cathmap                 = [];

[cathpos_file, solflag] = cathpos(sub);

% Get: filename parts for each data block [orgnames], [surface] of
% stimulation (EPI, LVEndo, RVEndo), position of stimulation electrode
% [stimelectrodes] and catheter positions [actmaps].
if solflag,
    [orgnames,surface,stimelectrodes,actmaps]   = readcatheterpos(cathpos_file);
    solutionflag                                = true;
else
    orgnames        = {};
    stimelectrodes  = [];
    solutionflag    = false;
end

%%
if solutionflag,
    % flatten orgnames, multiple files with same position
    flatorgnames    = {};
    flatelectrodes  = emptystruct(stimelectrodes);
    flatactmaps     = emptystruct(actmaps);
    flatsurface     = [];
    
    for i = 1:length(orgnames),
        if iscell(orgnames{i}),
            for j = 1:length(orgnames{i}),
                flatorgnames{end+1}     = orgnames{i}{j};
                flatelectrodes(end+1)   = stimelectrodes(i);
                
                if exist('surface','var'), flatsurface{end+1} = surface{i}; end
                
                flatactmaps(end+1)      = actmaps(i);
            end
        else
            flatorgnames{end+1}     = orgnames{i};
            flatelectrodes(end+1)   = stimelectrodes(i);
            if exist('surface','var'), flatsurface{end+1} = surface{i}; end
            
        end
    end
    
    for i = 1:length(flatorgnames),
        if strfind(filename,flatorgnames{i}),
            stimpos = flatelectrodes(i);
            
            if exist('surface','var'), cathsurf = flatsurface{i}; end
            
            cathmap     = flatactmaps(i);
            bsmfname    = fullfile(datapath,[flatorgnames{i} '_BSM.mat']);
            
            if exist(bsmfname,'file'), bsm = loadmat(bsmfname); end
            
            break
        end
    end
else
    stimpos = nan(1,3);
    endostr = regexp(filename,'Endo','match');
    
    if isempty(endostr),    endostr = '';   else endostr = endostr{1};  end
    
    epistr = regexp(filename,'Epi','match');
    
    if isempty(epistr),     epistr = '';    else epistr = epistr{1};    end
    
    lvstr = regexp(filename,'LV','match');
    
    if isempty(lvstr),      lvstr = '';     else lvstr = lvstr{1};      end
    
    rvstr = regexp(filename,'RV','match');
    
    if isempty(rvstr),      rvstr = '';     else rvstr = rvstr{1};      end
    
    cathsurf = [lvstr rvstr endostr epistr];
end

% Sub-function:
function S = emptystruct(obj)

% get fielnames:
fnlist  = fieldnames(obj);
evalstr = '';

% sum all field names in evalstr:
for i = 1:length(fnlist), evalstr = [evalstr sprintf(',''%s'',{}',fnlist{i})]; end

S       = eval(sprintf('struct(%s)',evalstr(2:end)));