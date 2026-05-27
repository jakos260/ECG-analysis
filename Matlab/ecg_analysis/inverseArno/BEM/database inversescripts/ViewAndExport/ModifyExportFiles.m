clearvars
[FileNameIn,PathNameIn]=uigetfile({'*.mat';'*.*'}, 'Select export files',...
    '/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/','MultiSelect','on');
if ~iscell(FileNameIn)
    FileNameIn={FileNameIn};
end

for i=1:length(FileNameIn)
    clear DATA D mlas
    load(fullfile(PathNameIn,FileNameIn{i}));
    %     for r=[1 3 5 6 32 57]
    %         if ~find(DATA.remove==r)
    %             DATA.remove(end+1)=r;
    %         end
    %     end
    DATA.remove([31 32 63 64])=1;
% 
% D.signals=D.signals*4;
    save(fullfile(PathNameIn,FileNameIn{i}),'DATA','D','mlas');
end

