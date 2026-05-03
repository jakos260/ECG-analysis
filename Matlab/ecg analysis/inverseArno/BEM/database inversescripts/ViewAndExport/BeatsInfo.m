% BeatsInfo
% Plot timing and source file for beats in STW export file
clearvars
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select Biosemi export files',...
    '/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/BioSemi/export/AVG','MultiSelect','on');
if ~iscell(FileName)
    FileName={FileName};
end
for i=1:length(FileName)
    d=load(fullfile(PathName,FileName{i}));
    [~,bdfname,ext]=fileparts(d.D.filepath);
    bdfname=[bdfname ext];
    fprintf('File %s \n(extracted from %s)\n',FileName{i},bdfname);
    selbeats=find(d.DATA.SELBEATS);
    for j=1:length(selbeats)
        tbeat=d.DATA.BEATS(selbeats(j))/1000;
        fprintf('\t Beat %d at %0.1f (%0.1f)\n',selbeats(j),tbeat,d.D.startpos+tbeat);
    end
    
end