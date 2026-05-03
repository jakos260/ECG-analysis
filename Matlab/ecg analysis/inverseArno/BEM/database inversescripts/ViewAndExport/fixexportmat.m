clear all
clear global
close all
egmrange=65:76;
removerange=[1,5,25,33,39,45,57,63:68];

[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select STW export files to fix.','/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export','MultiSelect','on');
for i=1:length(FileName)
    matpath=fullfile(PathName,FileName{i});
    clear D DATA mlas
    load(matpath,'D','DATA','mlas');
    if ~isfield(D,'EGM') || size(D.signals,2)<max(egmrange)
        comment=D.comment;
        filepath=D.filepath;
        D=Readbdf(D.filepath,1:max(egmrange),D.startpos,D.endpos,'noui');
        D.comment=comment;
        D.filepath=filepath;
        DATA.EGM=interp1(1:size(D.signals,1),D.signals(:,egmrange),1:1/D.Ts:size(D.signals,1),'cubic')';
    end
    DATA.remove(removerange)=1;
    save(matpath,'D','DATA','mlas');
end

