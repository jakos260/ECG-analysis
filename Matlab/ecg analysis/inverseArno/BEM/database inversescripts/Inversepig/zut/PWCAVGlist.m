if ~exist('FileName','var') || isempty(FileName)
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select export files '...
        ,'/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export','MultiSelect','on');
end
if ~iscell(FileName)
    FileName={FileName};
end
for ifile=1:length(FileName);
    FilePath=(fullfile(PathName,FileName{ifile}));
S=load(FilePath);
fprintf('File %s: PWC: %d; AVG: %d PWCAVG: %d \n',FileName{ifile},isfield(S.DATA,'BSMOUTPWC'),isfield(S.DATA,'AVERAGE'),isfield(S.DATA,'AVERAGEPWC'));
end
